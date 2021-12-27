/*
 * Copyright 2019 Johns Hopkins University
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

#include <cmath>
#include "lm/ClassFactory.h"
#include "lm/Math.h"
#include "lm/Print.h"
#include "lm/Tune.h"
#include "lm/cme/FluctuatingNRSolver.h"
#include "lm/io/TrajectoryState.pb.h"
#include "lm/io/StochasticProcessTimeSeries.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/rng/RandomGenerator.h"
#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"
#include "robertslab/pbuf/NDArray.pb.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using lm::rng::RandomGenerator;
using robertslab::pbuf::NDArraySerializer;

namespace lm {
namespace cme {

bool FluctuatingNRSolver::registered=FluctuatingNRSolver::registerClass();

bool FluctuatingNRSolver::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::me::MESolver","lm::cme::FluctuatingNRSolver",&FluctuatingNRSolver::allocateObject);
    return true;
}

void* FluctuatingNRSolver::allocateObject()
{
    return new FluctuatingNRSolver();
}

FluctuatingNRSolver::FluctuatingNRSolver()
:NextReactionSolver(static_cast<RandomGenerator::Distributions>(RandomGenerator::NORMAL)),
normRngValues(NULL),
numberNoisyReactions(0),noisyReactions(NULL),
nextProcessUpdateTime(std::numeric_limits<double>::infinity()),processUpdateInterval(std::numeric_limits<double>::infinity()),numberProcesses(0),processes(NULL),writeProcessTimeSeries(false)
{
    // Allocate the RNG buffers.
    allocateRngBuffers();
}

FluctuatingNRSolver::~FluctuatingNRSolver()
{
    if (noisyReactions != NULL) delete[] noisyReactions; noisyReactions = NULL;
    for (size_t i=0; i<numberProcesses; i++) delete processes[i];
    if (processes != NULL) delete[] processes; processes = NULL;
    numberProcesses = 0;
    deallocateRngBuffers();
}

void FluctuatingNRSolver::allocateRngBuffers()
{
    NextReactionSolver::allocateRngBuffers();

    if (normRngValues == NULL)
    {
        normRngValues = new double[TUNE_LOCAL_RNG_CACHE_SIZE];
        nextNormRngValue = TUNE_LOCAL_RNG_CACHE_SIZE;
    }
}

void FluctuatingNRSolver::deallocateRngBuffers()
{
    NextReactionSolver::deallocateRngBuffers();
    if (normRngValues != NULL) delete[] normRngValues; normRngValues = NULL;
    nextNormRngValue = 0;
}

void FluctuatingNRSolver::setOutputOptions(const lm::input::OutputOptions& outputOptions)
{
    NextReactionSolver::setOutputOptions(outputOptions);
    writeProcessTimeSeries = outputOptions.has_stochastic_process_write_interval();
    processWriteInterval = writeProcessTimeSeries?outputOptions.stochastic_process_write_interval():std::numeric_limits<double>::infinity();
}

void FluctuatingNRSolver::setReactionModel(const lm::input::ReactionModel& rm)
{
    NextReactionSolver::setReactionModel(rm);

    if (!rm.has_noise_model()) THROW_EXCEPTION(lm::InputException, "FluctuatingNRSolver requires a noise model to be associated with the reaction model");

    // Free any memory used by the previous model.
    numberNoisyReactions = 0;
    if (noisyReactions != NULL) delete[] noisyReactions; noisyReactions = NULL;
    noisyReactionDependencies.clear();
    noisyReactionOriginalConstants.clear();
    for (size_t i=0; i<numberProcesses; i++) delete processes[i];
    if (processes != NULL) delete processes; processes = NULL;
    numberProcesses = 0;

    // Read the model.
    numberProcesses = rm.noise_model().number_processes();
    processUpdateInterval = rm.noise_model().process_update_interval();
    ndarray<uint32_t> T = NDArraySerializer::deserialize<uint32_t>(rm.noise_model().process_types());
    ndarray<double> K = NDArraySerializer::deserialize<double>(rm.noise_model().process_parameters());
    ndarray<uint32_t> D = NDArraySerializer::deserialize<uint32_t>(rm.noise_model().reaction_dependencies());
    if (T.shape.len != 1 || T.shape[0] != numberProcesses) THROW_EXCEPTION(lm::InputException, "Invalid shape for process_types");
    if (K.shape.len != 2 || K.shape[0] != numberProcesses || K.shape[1] != 10) THROW_EXCEPTION(lm::InputException, "Invalid shape for process_parameters");
    if (D.shape.len != 3 || D.shape[0] != reactionModel->numberReactions || D.shape[1] != 10 || D.shape[2] != numberProcesses) THROW_EXCEPTION(lm::InputException, "Invalid shape for reaction_dependencies");

    // Create the new processes.
    processes = new StochasticProcess*[numberProcesses];
    for (uint i=0; i<numberProcesses; i++)
    {
        switch (T[i])
        {
        case 0:
            processes[i] = new OrnsteinUhlenbeckProcess(K[utuple(i,0U)],K[utuple(i,1)]);
            break;
        case 1:
            processes[i] = new LogOrnsteinUhlenbeckProcess(K[utuple(i,0U)],K[utuple(i,1)]);
            break;
        default:
            THROW_EXCEPTION(lm::InputException, "Invalid stochastic process type: %d", T[i]);
        }
    }

    // Get a list of reactions that should be updated when the noise is updated.
    noisyReactions = new uint[rm.number_reactions()];
    for (uint i=0; i<D.shape[0]; i++)
    {
        // Go through every entry for this reaction and look for dependencies.
        bool anyDependencies = false;
        for (uint j=0; j<D.shape[1]; j++)
        {
            for (uint k=0; k<D.shape[2]; k++)
            {
                if (D[utuple(i,j,k)] > 0)
                {
                    noisyReactionDependencies.push_back(utuple(i,j,k));
                    noisyReactionOriginalConstants.push_back(reactionModel->propensityFunctions[i]->getConstants()[j]);
                    Print::printf(Print::DEBUG, "Saved initial value of rate constants for noisy reaction r=%d k=%d value=%0.3f",i,j,reactionModel->propensityFunctions[i]->getConstants()[j]);
                    anyDependencies = true;
                }
            }
        }

        // If there were any dependencies in this reaction, add it to the list.
        if (anyDependencies)
            noisyReactions[numberNoisyReactions++] = i;
    }
}

void FluctuatingNRSolver::reset()
{
    NextReactionSolver::reset();

    // Reset the processes.
    for (size_t i=0; i<numberProcesses; i++)
    {
        // Update the rng cache.
        if (nextNormRngValue >= TUNE_LOCAL_RNG_CACHE_SIZE)
        {
            rng->getNormRandomDoubles(normRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
            nextNormRngValue=0;
        }

        // Reset the process with a new random state.
        processes[i]->reset(normRngValues[nextNormRngValue++]);
    }

    // Reset all of the propensity function constants.
    updateReactionPropensities();

    // Clear the time series.
    processTimeSeries.clear();

    nextProcessUpdateTime = std::numeric_limits<double>::infinity();
}

void FluctuatingNRSolver::getState(lm::io::TrajectoryState* state, uint trajectoryNumber)
{
    NextReactionSolver::getState(state, trajectoryNumber);
}

void FluctuatingNRSolver::setState(const lm::io::TrajectoryState& state, uint trajectoryNumber)
{
    NextReactionSolver::setState(state, trajectoryNumber);

    nextProcessUpdateTime = roundNextMultiple(time, processUpdateInterval);
}

void FluctuatingNRSolver::trajectoryStarted()
{
    NextReactionSolver::trajectoryStarted();

    // Get the interval for process values.
    if (writeProcessTimeSeries)
    {
        // Initialize the time series, filling in the first set of values if necessary.
        if (!previouslyStarted && writeInitialTrajectoryState)
        {
            double* buffer = new double[numberProcesses];
            for (size_t i=0; i<numberProcesses; i++) buffer[i] = processes[i]->value;
            processTimeSeries.initialize(numberProcesses, processWriteInterval, time, buffer, numberProcesses);
            delete[] buffer;
        }
        else
        {
            processTimeSeries.initialize(numberProcesses, processWriteInterval, time);
        }
    }

}

void FluctuatingNRSolver::timeUpdated()
{
    NextReactionSolver::timeUpdated();

    // If we are writing process time steps, write out any steps before this event occurred.
    if (writeProcessTimeSeries)
    {
        // Write time steps until the next write time is past the current time.
        double* buffer = new double[numberProcesses];
        for (size_t i=0; i<numberProcesses; i++) buffer[i] = processes[i]->value;
        processTimeSeries.append(buffer, numberProcesses, time);
        delete[] buffer;
    }
}

void FluctuatingNRSolver::trajectoryFinished()
{
    NextReactionSolver::trajectoryFinished();

    // If we hit any limits, write out the final time if requested.
    if (status == lm::message::WorkUnitStatus::LIMIT_REACHED && writeFinalTrajectoryState)
    {
        // Write out the last time step, or ensure that it already has been.
        if (writeProcessTimeSeries)
        {
            double* buffer = new double[numberProcesses];
            for (size_t i=0; i<numberProcesses; i++) buffer[i] = processes[i]->value;
            processTimeSeries.appendFinal(buffer, numberProcesses, time);
            delete[] buffer;
        }
    }

    // If we have any processs time series data, add them to the output message.
    if (processTimeSeries.size() > 0)
    {
        // Mark that the message does contain some data.
        output->set_has_output(true);

        // Serialize the data into the message.
        lm::io::StochasticProcessTimeSeries* ds = output->mutable_stochastic_process_time_series();
        ds->set_trajectory_id(trajectoryId);
        processTimeSeries.serializeInto(ds->mutable_values(), ds->mutable_times());
    }
}

uint64_t FluctuatingNRSolver::generateTrajectory(uint64_t maxSteps)
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
    Print::printf(Print::DEBUG, "Running fluctuating next reaction simulation with %d species, %d reactions, %d limits.", reactionModel->numberSpecies, reactionModel->numberReactions, timeLimits.limits_size()+stateLimits.limits_size());
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
        double nexTime = reactionQueue->getReactionEvent(r).time;

        // If the next reaction time is before the process update time, perform the reaction.
        if (nexTime <= nextProcessUpdateTime)
        {
            // Update the time. If we are outside the time limits, stop the trajectory.
            if (performTimeIncrement(nexTime-time))
                break;

            // Update species counts. If we are outside of the state limits, stop the trajectory.
            if (performReactionEvent(r))
                break;

            // Update the reaction queue.
            updateReactionEventsInQueue(r);
        }

        // Otherwise, update the processes, update the queue, then try again.
        else
        {
            // Update the time until the next process update occurrs.
            if (performTimeIncrement(nextProcessUpdateTime-time))
                break;

            // Update the processes.
            for (size_t i=0; i<numberProcesses; i++)
            {
                // Update the rng cache.
                if (nextNormRngValue >= TUNE_LOCAL_RNG_CACHE_SIZE)
                {
                    rng->getNormRandomDoubles(normRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
                    nextNormRngValue=0;
                }

                // Get a normal rng.
                double normRng = normRngValues[nextNormRngValue++];

                // Update the process.
                processes[i]->update(processUpdateInterval, normRng);
            }

            // Increment the time to the next process update.
            nextProcessUpdateTime += processUpdateInterval;

            // Update the reaction propensities with the new values of the process.
            updateReactionPropensities();

            // Update the reaction queue using the new propensities.
            rescaleReactionEventsInQueue();
        }

    }
    PROF_END(PROF_SIM_EXECUTE);
    Print::printf(Print::DEBUG, "Generated trajectory through time %e.", time);

    // Call the trajectory finished method.
    trajectoryFinished();

    return steps;
}

void FluctuatingNRSolver::updateReactionPropensities()
{
    // Go through each dependency.
    for (size_t i=0; i<noisyReactionDependencies.size(); i++)
    {
        // Get the dependency.
        utuple d = noisyReactionDependencies[i];

        // Get the new process value.
        double v = processes[d[2]]->value;

        // Get the orignal rate constant.
        double k = noisyReactionOriginalConstants[i];

        // Calculate the new rate constant.
        double noisyk = v*k;

        // Update the rate constant.
        reactionModel->propensityFunctions[d[0]]->setConstant(d[1], noisyk);
    }
}

void FluctuatingNRSolver::rescaleReactionEventsInQueue()
{
    // Go through each reaction that has some noise.
    for (size_t i=0; i<numberNoisyReactions; i++)
    {
        // Get the reaction.
        uint r = noisyReactions[i];

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

StochasticProcess::StochasticProcess()
{
}

StochasticProcess::~StochasticProcess()
{
}

OrnsteinUhlenbeckProcess::OrnsteinUhlenbeckProcess(double variance, double tau)
:StochasticProcess(),variance(variance),tau(tau)
{
    _value = 0.0;
    value = 1.0;
}

OrnsteinUhlenbeckProcess::~OrnsteinUhlenbeckProcess()
{
}

void OrnsteinUhlenbeckProcess::reset(double normRng)
{
    _value = sqrt(variance)*normRng;
    value = ((_value+1.0)>=0.0)?(_value+1.0):(0.0);
}

void OrnsteinUhlenbeckProcess::update(double dt, double normRng)
{
    double mu = exp(-dt/tau);
    double sigma = sqrt(variance*(1-(mu*mu)));
    _value = _value*mu + sigma*normRng;
    value = ((_value+1.0)>=0.0)?(_value+1.0):(0.0);
}

LogOrnsteinUhlenbeckProcess::LogOrnsteinUhlenbeckProcess(double variance, double tau)
:variance(variance),tau(tau)
{
    _value = 0.0;
    value = log(_value);
}

LogOrnsteinUhlenbeckProcess::~LogOrnsteinUhlenbeckProcess()
{
}

void LogOrnsteinUhlenbeckProcess::reset(double normRng)
{
    _value = sqrt(variance)*normRng;
    value = exp(_value);
}

void LogOrnsteinUhlenbeckProcess::update(double dt, double normRng)
{
    double mu = exp(-dt/tau);
    double sigma = sqrt(variance*(1-(mu*mu)));
    _value = _value*mu + sigma*normRng;
    value = exp(_value);
}

}
}
