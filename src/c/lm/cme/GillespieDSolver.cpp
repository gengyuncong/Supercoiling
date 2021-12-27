/*
 * Copyright 2012-2019 Johns Hopkins University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Developed by: Roberts Group
 *               Johns Hopkins University
 *               http://biophysics.jhu.edu/roberts/
 *
 * Author(s): Elijah Roberts, Max Klein
 */

#include <cmath>
#include <cstdio>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "lm/ClassFactory.h"
#include "lm/Math.h"
#include "lm/Print.h"
#include "lm/Tune.h"
#include "lm/Types.h"
#include "lm/cme/CMESolver.h"
#include "lm/cme/GillespieDSolver.h"
#include "lm/io/DegreeAdvancementTimeSeries.pb.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/io/SpeciesTimeSeries.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/message/WorkUnitOutput.pb.h"
#include "lm/rng/RandomGenerator.h"
#include "lm/rng/XORShift.h"
#include "lm/thread/Thread.h"
#include "lm/thread/Worker.h"
#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using std::list;
using std::map;
using std::string;
using std::vector;
using lm::rng::RandomGenerator;

namespace lm {
namespace cme {

bool GillespieDSolver::registered=GillespieDSolver::registerClass();

bool GillespieDSolver::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::me::MESolver","lm::cme::GillespieDSolver",&GillespieDSolver::allocateObject);
    return true;
}

void* GillespieDSolver::allocateObject()
{
    return new GillespieDSolver();
}

GillespieDSolver::GillespieDSolver()
:CMESolver((RandomGenerator::Distributions)(RandomGenerator::EXPONENTIAL|RandomGenerator::UNIFORM)),
 rngValues(NULL),expRngValues(NULL),nextRngValue(0),propensities(NULL)
{
}

GillespieDSolver::~GillespieDSolver()
{
    // Free any state.
    if (propensities != NULL) delete[] propensities; propensities = NULL;
    deallocateRngBuffers();
}

void GillespieDSolver::allocateRngBuffers()
{
    if (expRngValues == NULL || rngValues == NULL)
    {
        rngValues = new double[TUNE_LOCAL_RNG_CACHE_SIZE];
        expRngValues = new double[TUNE_LOCAL_RNG_CACHE_SIZE];
        nextRngValue = TUNE_LOCAL_RNG_CACHE_SIZE;
    }
}

void GillespieDSolver::deallocateRngBuffers()
{
    if (expRngValues != NULL) delete[] expRngValues; expRngValues = NULL;
    if (rngValues != NULL) delete[] rngValues; rngValues = NULL;
    nextRngValue = 0;
}

void GillespieDSolver::reset()
{
    CMESolver::reset();

    // Make sure we have allocated the RNG buffers.
    allocateRngBuffers();

    // Free any previous state.
    if (propensities != NULL) delete[] propensities; propensities = NULL;

    // Allocate reaction propensities table.
    propensities = new double[reactionModel->numberReactions];

    // Set the propensities to their initial values.
    for (uint i=0; i<reactionModel->numberReactions; i++)
    {
        propensities[i] = 0.0;
    }
}

void GillespieDSolver::getState(lm::io::TrajectoryState* state, uint trajectoryNumber)
{
    CMESolver::getState(state, trajectoryNumber);
}

void GillespieDSolver::setState(const lm::io::TrajectoryState& state, uint trajectoryNumber)
{
    CMESolver::setState(state, trajectoryNumber);

    // Set the propensities to their initial values.
    updateAllPropensities();
}

uint64_t GillespieDSolver::generateTrajectory(uint64_t maxSteps)
{
    if (reactionModel == NULL) throw Exception("GillespieDSolver did not have a reaction model.");
    if (propensities == NULL) throw Exception("GillespieDSolver state was not initialized.");

    // Make sure we have propensity functions for every reaction.
    for (uint i=0; i<reactionModel->numberReactions; i++)
        if (reactionModel->propensityFunctions[i] == NULL)
            throw Exception("A reaction did not have a valid propensity function",i);

    // Create local copies of the data for efficiency.
    const uint numberReactions = reactionModel->numberReactions;

    // Initialize the total propensity.
    double totalPropensity = 0.0;
    for (uint i=0; i<numberReactions; i++) totalPropensity += propensities[i];

    // Call the trajectory started method.
    trajectoryStarted();

    // Run the direct method.
    Print::printf(Print::VERBOSE_DEBUG, "Running Gillespie direct simulation for %d steps with %d species, %d reactions, %d limits, %d order parameters.", maxSteps, reactionModel->numberSpecies, reactionModel->numberReactions, timeLimits.limits_size()+stateLimits.limits_size(), numberOrderParameters);
    PROF_BEGIN(PROF_SIM_EXECUTE);
    uint64_t steps=0;
    while (true)
    {
        // See if we have finished the steps.
        if (steps >= maxSteps)
        {
            status = lm::message::WorkUnitStatus::STEPS_FINISHED;
            break;
        }

        // Increment the step counter.
        steps++;

        // See if we need to update our rng caches.
        if (nextRngValue >= TUNE_LOCAL_RNG_CACHE_SIZE)
        {
            rng->getRandomDoubles(rngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
            rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
            nextRngValue=0;
        }

        // Get the random values for this iteration though the loop.
        double randomValue = rngValues[nextRngValue];
        double expRandomValue = expRngValues[nextRngValue];

        // Go to the next rng pair.
        nextRngValue++;

        // Update the time for the next reaction. If we are outside the time limits, stop the trajectory.
        if (performTimeIncrement(expRandomValue/totalPropensity))
            break;

        // Calculate which reaction it was.
        double rngPropensity = randomValue*totalPropensity;
        uint r=0;
        for (; r<(numberReactions-1); r++)
        {
            if (rngPropensity < propensities[r])
                break;
            else
                rngPropensity -= propensities[r];
        }

        // Update species counts. If we are outside of the state limits, stop the trajectory.
        if (performReactionEvent(r))
            break;

        // Update the propensites given the reaction that occurred.
        updatePropensities(r);

        // Recalculate the total propensity.
        totalPropensity = 0.0;
        for (uint i=0; i<numberReactions; i++) totalPropensity += propensities[i];

        // If the total propensity is zero, stop the simulation.
        if (totalPropensity <= 0)
        {
            setTrajectoryToMaxTime();
            break;
        }

        //if (steps%1000000 == 0) Print::printf(Print::DEBUG, "Step %d: time=%e, count=%d, prop=%e, totprop=%e, op=%0.2f",steps,time,speciesCounts[0],propensities[0],totalPropensity, orderParameterValues[0]);
        //if ((speciesCountsPrevious[0]+speciesCountsPrevious[1]+speciesCountsPrevious[2]) != 1 || (speciesCounts[0]+speciesCounts[1]+speciesCounts[2]) != 1) Print::printf(Print::DEBUG, "CRASH Step %d: time=%e, prevcount=%d,%d,%d,%d,%d,%d,%d, count=%d,%d,%d,%d,%d,%d,%d, prop=%e, totprop=%e, op=%0.2f",steps,time,speciesCountsPrevious[0],speciesCountsPrevious[1],speciesCountsPrevious[2],speciesCountsPrevious[3],speciesCountsPrevious[4],speciesCountsPrevious[5],speciesCountsPrevious[6],speciesCounts[0],speciesCounts[1],speciesCounts[2],speciesCounts[3],speciesCounts[4],speciesCounts[5],speciesCounts[6],propensities[0],totalPropensity, orderParameterValues[0]);
        //else Print::printf(Print::DEBUG, "Step %d: time=%e, prevcount=%d,%d,%d,%d,%d,%d,%d, count=%d,%d,%d,%d,%d,%d,%d, prop=%e, totprop=%e, op=%0.2f",steps,time,speciesCountsPrevious[0],speciesCountsPrevious[1],speciesCountsPrevious[2],speciesCountsPrevious[3],speciesCountsPrevious[4],speciesCountsPrevious[5],speciesCountsPrevious[6],speciesCounts[0],speciesCounts[1],speciesCounts[2],speciesCounts[3],speciesCounts[4],speciesCounts[5],speciesCounts[6],propensities[0],totalPropensity, orderParameterValues[0]);
    }
    PROF_END(PROF_SIM_EXECUTE);

    Print::printf(Print::VERBOSE_DEBUG, "Generated trajectory through time %e.", time);

    // Call the trajectory finished method.
    trajectoryFinished();

    return steps;
}

void GillespieDSolver::updateAllPropensities()
{
    // Update the propensities.
    for (uint i=0; i<reactionModel->numberReactions; i++)
    {
        propensities[i] = reactionModel->propensityFunctions[i]->calculate(time, speciesCounts, reactionModel->numberSpecies);
    }
}

void GillespieDSolver::updatePropensities(uint sourceReaction)
{
    // Update the propensities of the dependent reactions.
    for (uint i=0; i<reactionModel->numberDependentReactions[sourceReaction]; i++)
    {
        uint r = reactionModel->dependentReactions[sourceReaction][i];
        propensities[r] = reactionModel->propensityFunctions[r]->calculate(time, speciesCounts, reactionModel->numberSpecies);
    }
}

}
}
