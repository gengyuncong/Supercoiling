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
 * Author(s): Elijah Roberts
 */

#include <limits>

#include "lm/ClassFactory.h"
#include "lm/Math.h"
#include "lm/Print.h"
#include "lm/Tune.h"
#include "lm/Types.h"
#include "lm/cme/CMESolver.h"
#include "lm/cme/NextReactionSolver.h"
#include "lm/io/DegreeAdvancementTimeSeries.pb.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/io/SpeciesTimeSeries.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/message/WorkUnitOutput.pb.h"
#include "lm/reaction/ReactionQueue.h"
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
using lm::reaction::ReactionQueue;
using lm::rng::RandomGenerator;

namespace lm {
namespace cme {

bool NextReactionSolver::registered=NextReactionSolver::registerClass();

bool NextReactionSolver::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::me::MESolver","lm::cme::NextReactionSolver",&NextReactionSolver::allocateObject);
    return true;
}

void* NextReactionSolver::allocateObject()
{
    return new NextReactionSolver();
}

NextReactionSolver::NextReactionSolver()
:CMESolver((RandomGenerator::Distributions)(RandomGenerator::EXPONENTIAL)),
expRngValues(NULL),reactionQueue(NULL)
{
    // Allocate the RNG buffers.
    allocateRngBuffers();
}

NextReactionSolver::NextReactionSolver(RandomGenerator::Distributions neededDists)
:CMESolver((RandomGenerator::Distributions)(RandomGenerator::EXPONENTIAL|neededDists)),
expRngValues(NULL),reactionQueue(NULL)
{
    // Allocate the RNG buffers.
    allocateRngBuffers();
}

NextReactionSolver::~NextReactionSolver()
{
    if (reactionQueue != NULL) delete reactionQueue; reactionQueue = NULL;
    deallocateRngBuffers();
}

void NextReactionSolver::allocateRngBuffers()
{
    if (expRngValues == NULL)
    {
        expRngValues = new double[TUNE_LOCAL_RNG_CACHE_SIZE];
        nextExpRngValue = TUNE_LOCAL_RNG_CACHE_SIZE;
    }
}

void NextReactionSolver::deallocateRngBuffers()
{
    if (expRngValues != NULL) delete[] expRngValues; expRngValues = NULL;
    nextExpRngValue = 0;
}

void NextReactionSolver::reset()
{
    CMESolver::reset();

    // Delete any previous reaction queue.
    if (reactionQueue != NULL) delete reactionQueue; reactionQueue = NULL;

    // Create the reaction queue.
    reactionQueue = new ReactionQueue(reactionModel->numberReactions);
}

void NextReactionSolver::getState(lm::io::TrajectoryState* state, uint trajectoryNumber)
{
    CMESolver::getState(state, trajectoryNumber);
}

void NextReactionSolver::setState(const lm::io::TrajectoryState& state, uint trajectoryNumber)
{
    CMESolver::setState(state, trajectoryNumber);

    // Initialize the reaction queue.
    updateAllReactionEventsInQueue();
}

uint64_t NextReactionSolver::generateTrajectory(uint64_t maxSteps)
{
    if (reactionModel == NULL) throw Exception("NextReactionSolver did not have a reaction model.");
    if (reactionQueue == NULL) throw Exception("NextReactionSolver state was not initialized.");

    // Make sure we have propensity functions for every reaction.
    for (uint i=0; i<reactionModel->numberReactions; i++)
        if (reactionModel->propensityFunctions[i] == NULL)
            throw Exception("A reaction did not have a valid propensity function",i);

    // Call the trajectory started method.
    trajectoryStarted();

    // Run the next reaction method.
    Print::printf(Print::DEBUG, "Running next reaction simulation with %d species, %d reactions, %d limits.", reactionModel->numberSpecies, reactionModel->numberReactions, timeLimits.limits_size()+stateLimits.limits_size());
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

        // Get the next reaction.
        uint r = reactionQueue->getNextReaction();
        double newTime = reactionQueue->getReactionEvent(r).time;

        // If the the time until the next reaction is infinite, stop the simulation.
        if (newTime == std::numeric_limits<double>::infinity())
        {
            setTrajectoryToMaxTime();
            break;
        }


        // Update the time. If we are outside the time limits, stop the trajectory.
        if (performTimeIncrement(newTime-time))
            break;

        // Update species counts. If we are outside of the state limits, stop the trajectory.
        if (performReactionEvent(r))
            break;

        // Update the reaction queue.
        updateReactionEventsInQueue(r);

    }
    PROF_END(PROF_SIM_EXECUTE);
    Print::printf(Print::DEBUG, "Generated trajectory through time %e.", time);

    // Call the trajectory finished method.
    trajectoryFinished();

    return steps;
}

void NextReactionSolver::updateAllReactionEventsInQueue()
{
    // Update all of the reactions.
    for (uint r=0; r<reactionModel->numberReactions; r++)
    {
        double propensity = reactionModel->propensityFunctions[r]->calculate(time, speciesCounts, reactionModel->numberSpecies);
        double newTime = std::numeric_limits<double>::infinity();
        if (propensity > 0.0)
        {
            if (nextExpRngValue >= TUNE_LOCAL_RNG_CACHE_SIZE)
            {
                rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
                nextExpRngValue=0;
            }
            newTime = time+expRngValues[nextExpRngValue++]/propensity;
        }
        reactionQueue->updateReactionEvent(r, newTime, propensity);
    }
}

void NextReactionSolver::updateReactionEventsInQueue(uint sourceReaction)
{
    // Update the source reaction.
    {
        double propensity = reactionModel->propensityFunctions[sourceReaction]->calculate(time, speciesCounts, reactionModel->numberSpecies);
        double newTime = std::numeric_limits<double>::infinity();
        if (propensity > 0.0)
        {
            if (nextExpRngValue >= TUNE_LOCAL_RNG_CACHE_SIZE)
            {
                rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
                nextExpRngValue=0;
            }
            newTime = time+expRngValues[nextExpRngValue++]/propensity;
        }
        reactionQueue->updateReactionEvent(sourceReaction, newTime, propensity);
    }

    // Update any dependent reactions.
    for (uint i=0; i<reactionModel->numberDependentReactions[sourceReaction]; i++)
    {
        uint r = reactionModel->dependentReactions[sourceReaction][i];

        // If the source reaction depends on itself, skip it.
        if (r == sourceReaction) continue;

        // Calculate the new propensity.
        double newPropensity = reactionModel->propensityFunctions[r]->calculate(time, speciesCounts, reactionModel->numberSpecies);
        double newTime = std::numeric_limits<double>::infinity();

        // If there is some propensity, figure out the new time.
        if (newPropensity > 0.0)
        {
            #if TUNE_NRM_REUSE_PROPENSITIES == 1
            double oldTime = reactionQueue->getReactionEvent(r).time;
            double oldPropensity = reactionQueue->getReactionEvent(r).propensity;

            // If this reaction had a previous time, reuse it.
            if (oldPropensity > 0.0)
            {
                newTime = time+((oldPropensity/newPropensity)*(oldTime-time));
            }

            // Otherwise choose a new time.
            else
            {
            #endif

                if (nextExpRngValue >= TUNE_LOCAL_RNG_CACHE_SIZE)
                {
                    rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
                    nextExpRngValue=0;
                }
                newTime = time+expRngValues[nextExpRngValue++]/newPropensity;
            #if TUNE_NRM_REUSE_PROPENSITIES == 1
            }
            #endif
        }
        reactionQueue->updateReactionEvent(r, newTime, newPropensity);
    }
}

}
}
