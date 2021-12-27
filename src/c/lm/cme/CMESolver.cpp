/*
 * University of Illinois Open Source License
 * Copyright 2008-2011 Luthey-Schulten Group,
 * Copyright 2012-2019 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
 *
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to
 * do so, subject to the following conditions:
 *
 * - Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimers.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimers in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Elijah Roberts, Max Klein
 */

#include <algorithm>
#include <cmath>
#include <limits>
#include <list>
#include <map>
#include <string>
#if defined(MACOSX)
#elif defined(LINUX)
#include <time.h>
#endif

#include "lm/Math.h"
#include "lm/Print.h"
#include "lm/Tune.h"
#include "lm/Types.h"
#include "lm/cme/CMESolver.h"
#include "lm/cme/ReactionModel.h"
#include "lm/types/OrderParameters.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/io/BarrierCrossingTimes.pb.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/io/SpeciesCounts.pb.h"
#include "lm/io/TrajectoryState.pb.h"
#include "lm/limit/LimitCheckMacros.h"
#include "lm/me/PropensityFunction.h"
#include "lm/me/TrajectoryHistogram.h"
#include "lm/message/WorkUnitStatus.pb.h"
#include "lm/oparam/OrderParameterFunction.h"
#include "lm/rng/RandomGenerator.h"
#include "lm/rng/XORShift.h"
#ifdef OPT_CUDA
#include "lm/rng/XORWow.h"
#endif
#include "lm/thread/Thread.h"
#include "lm/thread/Worker.h"
#include "lm/types/TrajectoryLimits.pb.h"
#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"
#include "robertslab/pbuf/NDArraySerializer.h"

#if defined(MACOSX)
#elif defined(LINUX)
#include <time.h>
#endif
#ifdef OPT_CUDA
#include "lm/rng/XORWow.h"
#endif

using std::list;
using std::map;
using std::string;
using robertslab::pbuf::NDArraySerializer;

namespace lm {
namespace cme {

CMESolver::CMESolver(RandomGenerator::Distributions neededDists)
:writeInitialTrajectoryState(false),writeFinalTrajectoryState(false),
 writeSpeciesTimeSeries(false),writeDegreeAdvancementTimeSeries(false),writeOrderParameterTimeSeries(false),
 degreeAdvancementWriteInterval(0.0),orderParameterWriteInterval(0.0),speciesWriteInterval(0.0),

 neededDists(neededDists),rng(NULL),reactionModel(NULL),
 hasUpdateReactionCountsListeners(false),hasUpdateSpeciesCountsListeners(false),numberOrderParameters(0),
 orderParameterFunctions(NULL),output(new lm::message::WorkUnitOutput()),status(lm::message::WorkUnitStatus::NONE),
 trajectoryId(std::numeric_limits<uint64_t>::max()),previouslyStarted(false),
 numberFptSpecies(0),numberFptOrderParameters(0),fptSpeciesValues(NULL),fptOrderParameterValues(NULL),
 totalSteps(0),speciesCounts(NULL),speciesCountsPrevious(NULL),trackDegreeAdvancements(false),degreeAdvancements(NULL),orderParameterValues(NULL),orderParameterValuesPrevious(NULL),
 time(0.0),timeStep(0.0),
 numberTrackingBarriers(0),trackingBarrierPriorCrossings(NULL),trackingBarrierTimesCounts(NULL),trackingBarrierTimesTimes(NULL),trackingBarrierTimesSteps(NULL),
 histogramBeginTime(0.0),histogramEndTime(std::numeric_limits<double>::infinity())
{
}

CMESolver::~CMESolver()
{
    // Free any model memory.
    if (reactionModel != NULL) delete reactionModel; reactionModel = NULL;

    // Free any output memory.
    if (output != NULL) delete output; output = NULL;

    // Free any memory associated with the state.
    if (degreeAdvancements != NULL) delete[] degreeAdvancements; degreeAdvancements = NULL;
    if (orderParameterFunctions != NULL)
    {
        for (int i=0; i<numberOrderParameters; i++) delete orderParameterFunctions[i];
        delete[] orderParameterFunctions; orderParameterFunctions = NULL;
    }
    if (orderParameterValues != NULL) delete orderParameterValues; orderParameterValues = NULL;
    if (orderParameterValuesPrevious != NULL) delete orderParameterValuesPrevious; orderParameterValuesPrevious = NULL;
    if (speciesCounts != NULL) delete[] speciesCounts; speciesCounts = NULL;
    if (speciesCountsPrevious != NULL) delete[] speciesCountsPrevious; speciesCountsPrevious = NULL;

    // Free any other memory.
    if (rng != NULL) delete rng; rng = NULL;
    if (fptSpeciesValues != NULL) delete[] fptSpeciesValues; fptSpeciesValues = NULL;
    if (fptOrderParameterValues != NULL) delete[] fptOrderParameterValues; fptOrderParameterValues = NULL;

    if (trackingBarrierPriorCrossings != NULL) delete[] trackingBarrierPriorCrossings; trackingBarrierPriorCrossings = NULL;
    if (trackingBarrierTimesCounts != NULL) delete[] trackingBarrierTimesCounts; trackingBarrierTimesCounts = NULL;
    if (trackingBarrierTimesTimes != NULL) delete[] trackingBarrierTimesTimes; trackingBarrierTimesTimes = NULL;
    if (trackingBarrierTimesSteps != NULL) delete[] trackingBarrierTimesSteps; trackingBarrierTimesSteps = NULL;

    // Free all of the histograms.
    for (size_t i=0; i<tilingHistograms.size(); i++)
        delete tilingHistograms[i];
    tilingHistograms.clear();
}

void CMESolver::setComputeResources(vector<int> cpus, vector<int> gpus)
{
    lm::me::MESolver::setComputeResources(cpus, gpus);

    // Create the appropriate RNG.
    if (neededDists != RandomGenerator::NONE)
    {
        if (gpus.size() > 0)
        {
            #ifdef OPT_CUDA
            // Create the cuda based rng.
            rng = new lm::rng::XORWow(gpus[0], 0, 0, neededDists);
            #endif
        }

        if (rng == NULL)
        {
            rng = new lm::rng::XORShift(0, 0);
        }
    }
}

void CMESolver::setLimits(const lm::types::TrajectoryLimits& newLimits)
{
    // Clear the limits.
    timeLimits.Clear();
    stateLimits.Clear();

    // Save the new limits into the correct lists.
    for (int i=0; i<newLimits.limits().size(); i++)
    {
        const lm::types::TrajectoryLimit& l = newLimits.limits(i);
        if (l.type() == lm::types::TrajectoryLimit::TIME)
            timeLimits.add_limits()->CopyFrom(l);
        else
            stateLimits.add_limits()->CopyFrom(l);
    }
}

void CMESolver::setBarriers(const lm::types::TrajectoryBarriers& newBarriers)
{
    // Clear the barriers.
    reflectingBarriers.Clear();
    trackingBarriers.Clear();

    // Delete the vectors for tracking barrier crossings.
    numberTrackingBarriers = 0;
    if (trackingBarrierPriorCrossings != NULL) delete[] trackingBarrierPriorCrossings; trackingBarrierPriorCrossings = NULL;
    if (trackingBarrierTimesCounts != NULL) delete[] trackingBarrierTimesCounts; trackingBarrierTimesCounts = NULL;
    if (trackingBarrierTimesTimes != NULL) delete[] trackingBarrierTimesTimes; trackingBarrierTimesTimes = NULL;
    if (trackingBarrierTimesSteps != NULL) delete[] trackingBarrierTimesSteps; trackingBarrierTimesSteps = NULL;

    // Save the new barriers.
    for (int i=0; i<newBarriers.barrier().size(); i++)
    {
        const lm::types::TrajectoryBarrier& b = newBarriers.barrier(i);
        if (b.behavior() == lm::types::TrajectoryBarrier::REFLECTING || b.behavior() == lm::types::TrajectoryBarrier::REFLECTING_INCREASING || b.behavior() == lm::types::TrajectoryBarrier::REFLECTING_DECREASING)
            reflectingBarriers.add_barrier()->CopyFrom(b);
        else if (b.behavior() == lm::types::TrajectoryBarrier::TRACKING || b.behavior() == lm::types::TrajectoryBarrier::TRACKING_INCREASING_INCLUSIVE || b.behavior() == lm::types::TrajectoryBarrier::TRACKING_DECREASING_INCLUSIVE || b.behavior() == lm::types::TrajectoryBarrier::TRACKING_INCREASING_EXCLUSIVE || b.behavior() == lm::types::TrajectoryBarrier::TRACKING_DECREASING_EXCLUSIVE)
            trackingBarriers.add_barrier()->CopyFrom(b);
    }

    // Create the vectors for tracking the barrier crossings.
    numberTrackingBarriers = static_cast<size_t>(trackingBarriers.barrier().size());
    if (numberTrackingBarriers > 0)
    {
        trackingBarrierPriorCrossings = new uint64_t[numberTrackingBarriers];
        trackingBarrierTimesCounts = new vector<int32_t>[numberTrackingBarriers];
        trackingBarrierTimesTimes = new vector<double>[numberTrackingBarriers];
        trackingBarrierTimesSteps = new vector<uint64_t>[numberTrackingBarriers];
    }

    // Mark that we have a listener to update whenever the speices counts changes.
    hasUpdateSpeciesCountsListeners = true;
}

void CMESolver::setOrderParameters(const lm::types::OrderParameters& ops)
{
    if (orderParameterFunctions != NULL)
    {
        for (int i=0;i<numberOrderParameters;i++) delete orderParameterFunctions[i];
        delete[] orderParameterFunctions; orderParameterFunctions = NULL;
    }
    if (orderParameterValues != NULL) delete orderParameterValues; orderParameterValues = NULL;
    if (orderParameterValuesPrevious != NULL) delete orderParameterValuesPrevious; orderParameterValuesPrevious = NULL;

    // Allocate space for the order parameters.
    numberOrderParameters = ops.order_parameter().size();
    orderParameterFunctions = new lm::oparam::OrderParameterFunction*[reactionModel->numberSpecies];
    orderParameterValues = new double[numberOrderParameters];
    orderParameterValuesPrevious = new double[numberOrderParameters];

    // Create the order parameter functions.
    lm::oparam::OrderParameterFunctionFactory fs;
    for (int i=0; i<numberOrderParameters; i++)
        orderParameterFunctions[i] = fs.createOrderParameterFunction(ops.order_parameter(i));

    // Mark that we have a listener to update whenever the speices counts changes.
    hasUpdateSpeciesCountsListeners = true;
}

void CMESolver::setTilings(const lm::types::Tilings& t)
{
    tilings.CopyFrom(t);
}

void CMESolver::setOutputOptions(const lm::input::OutputOptions& outputOptions)
{
//    workUnitCondenseOutput = outputOptions.condense_output();
    writeInitialTrajectoryState = outputOptions.write_initial_trajectory_state();
    writeFinalTrajectoryState = outputOptions.write_final_trajectory_state();

    //
    // Time series output options.
    //

    writeSpeciesTimeSeries = outputOptions.has_species_write_interval();
    speciesWriteInterval = writeSpeciesTimeSeries?outputOptions.species_write_interval():std::numeric_limits<double>::infinity();

    writeDegreeAdvancementTimeSeries = outputOptions.has_degree_advancement_write_interval();
    degreeAdvancementWriteInterval = writeDegreeAdvancementTimeSeries?outputOptions.degree_advancement_write_interval():std::numeric_limits<double>::infinity();

    writeOrderParameterTimeSeries = outputOptions.has_order_parameter_write_interval();
    orderParameterWriteInterval = writeOrderParameterTimeSeries?outputOptions.order_parameter_write_interval():std::numeric_limits<double>::infinity();

    //
    // Histogram output options.
    //

    // Delete any existing tiling histograms.
    for (size_t i=0; i<tilingHistograms.size(); i++)
        delete tilingHistograms[i];
    tilingHistograms.clear();

    // Create a histogram for each tiling we are tracking.
    for (int i=0; i<outputOptions.tiling_to_histogram().size(); i++)
    {
        size_t tilingIndex = outputOptions.tiling_to_histogram(i);
        if (tilingIndex >= static_cast<size_t>(tilings.tiling().size())) THROW_EXCEPTION(lm::RuntimeException, "tiling to for histogramming %d is greater than the number of tilings %d", tilingIndex, tilings.tiling().size());
        ndarray<int32_t> opIndices = NDArraySerializer::deserialize<int32_t>(tilings.tiling(i).order_parameter_indices());
        ndarray<double> opEdges = NDArraySerializer::deserialize<double>(tilings.tiling(i).edges());
        lm::me::TrajectoryHistogram<double>* histogram = new lm::me::TrajectoryHistogram<double>(opIndices, opEdges);
        tilingHistograms.push_back(histogram);
    }

    // Set the start and end times for histogramming.
    histogramBeginTime = outputOptions.has_histogram_begin_time()?outputOptions.histogram_begin_time():0.0;
    histogramEndTime = outputOptions.has_histogram_end_time()?outputOptions.histogram_end_time():std::numeric_limits<double>::infinity();
}

void CMESolver::setReactionModel(const lm::input::ReactionModel& rm)
{
    if (rm.number_reactions() != (uint)rm.reaction_size()) throw InvalidArgException("rm", "number of reaction does not agree with reaction list size");

    // Set the new reaction model.
    if (reactionModel != NULL) delete reactionModel;
    reactionModel = new ReactionModel(rm);

    // Allocate space for all of the possible degree advancement counts.
    if (degreeAdvancements != NULL) delete[] degreeAdvancements; degreeAdvancements = NULL;
    degreeAdvancements = new uint64_t[reactionModel->numberReactions];
    trackDegreeAdvancements = false;

    // Allocate space for the species counts.
    if (speciesCounts != NULL) delete[] speciesCounts; speciesCounts = NULL;
    if (speciesCountsPrevious != NULL) delete[] speciesCountsPrevious; speciesCountsPrevious = NULL;
    speciesCounts = new int[reactionModel->numberSpecies];
    speciesCountsPrevious = new int[reactionModel->numberSpecies];
}

void CMESolver::reset()
{
    lm::me::MESolver::reset();

    // Make sure we have a reaction model.
    if (reactionModel == NULL) throw Exception("Tried to reset state of CMESolver with no reaction model.");

    // Reset the total step count.
    totalSteps = 0;

    // Reset the species counts.
    for (uint i=0; i<reactionModel->numberSpecies; i++)
    {
        speciesCounts[i] = 0;
        speciesCountsPrevious[i] = 0;
    }

    // Reset the degree advancements.
    for (uint i=0; i<reactionModel->numberReactions; i++)
        degreeAdvancements[i] = 0;

    // Reset the order parameters.
    for (int i=0; i<numberOrderParameters; i++)
    {
        orderParameterValues[i] = 0.0;
        orderParameterValuesPrevious[i] = 0.0;
    }

    // Reset the time series lists.
    speciesTimeSeries.clear();
    degreeAdvancementTimeSeries.clear();
    orderParameterTimeSeries.clear();

    // Reset the fpt tracking lists.
    numberFptSpecies = 0;
    numberFptOrderParameters = 0;
    if (fptSpeciesValues != NULL) delete[] fptSpeciesValues; fptSpeciesValues = NULL;
    if (fptOrderParameterValues != NULL) delete[] fptOrderParameterValues; fptOrderParameterValues = NULL;

    // Reset the limits reached.
    limitReached.Clear();

    // Reset the barrier tracking counts.
    for (size_t i=0; i<numberTrackingBarriers; i++)
    {
        trackingBarrierPriorCrossings[i] = 0;
        trackingBarrierTimesCounts[i].clear();
        trackingBarrierTimesTimes[i].clear();
        trackingBarrierTimesSteps[i].clear();
    }

    // Reset the histograms.
    for (size_t i=0; i<tilingHistograms.size(); i++)
        tilingHistograms[i]->clear();

    // Reset the output.
    if (output != NULL) delete output;
    output = new lm::message::WorkUnitOutput();

    // Reset the status.
    status = lm::message::WorkUnitStatus::NONE;
    trajectoryId = std::numeric_limits<uint64_t>::max();
    previouslyStarted = false;

    // Reset the tiling histograms list.
//    numberTilingHists = 0;
//    if (tilingHists!=NULL) delete[] tilingHists; tilingHists = NULL;

    // Reset the time.
    time = 0.0;
    timeStep = 0.0;
}

void CMESolver::getState(lm::io::TrajectoryState* state, uint trajectoryNumber)
{
    if (trajectoryNumber >= getSimultaneousTrajectories()) throw lm::InvalidArgException("trajectoryNumber", "exceeded the maximum number of simultaneous trajectories",trajectoryNumber,getSimultaneousTrajectories());

    // Get the trajectory id.
    state->set_trajectory_id(trajectoryId);

    // Get the previously started flag.
    state->set_trajectory_started(true);

    // Get the total step count.
    state->mutable_cme_state()->set_total_steps(totalSteps);

    // Get the time.
    state->mutable_cme_state()->mutable_species_counts()->add_time(time);

    // Get the species counts.
    state->mutable_cme_state()->mutable_species_counts()->set_trajectory_id(trajectoryId);
    state->mutable_cme_state()->mutable_species_counts()->set_number_species((int)reactionModel->numberSpecies);
    state->mutable_cme_state()->mutable_species_counts()->set_number_entries(1);
    for (int i=0; i<(int)reactionModel->numberSpecies; i++)
    {
        state->mutable_cme_state()->mutable_species_counts()->add_species_count(speciesCounts[i]);
        state->mutable_cme_state()->mutable_species_counts()->add_species_count_previous(speciesCountsPrevious[i]);
    }

    // Get the degree advancements.
    if (trackDegreeAdvancements)
    {
        NDArraySerializer::serializeInto(state->mutable_cme_state()->mutable_degree_advancements(), degreeAdvancements, utuple(reactionModel->numberReactions));
    }

    // Get the first passage times.
    for (int i=0; i<numberFptSpecies; i++)
    {
        fptSpeciesValues[i].serializeInto(state->mutable_cme_state()->add_first_passage_times());
    }
    for (int i=0; i<numberFptSpecies; i++)
    {
        //TODOfptOrderParameterValues[i].serializeInto(state->mutable_cme_state()->add_order_parameter_first_passage_times());
    }

    // Get the order parameter values.
    if (numberOrderParameters > 0)
    {
        state->mutable_cme_state()->mutable_order_parameter_values()->set_trajectory_id(trajectoryId);
        state->mutable_cme_state()->mutable_order_parameter_values()->set_number_order_parameters(numberOrderParameters);
        state->mutable_cme_state()->mutable_order_parameter_values()->set_number_entries(1);
        for (int i=0; i<numberOrderParameters; i++)
        {
            state->mutable_cme_state()->mutable_order_parameter_values()->add_order_parameter_values(orderParameterValues[i]);
            state->mutable_cme_state()->mutable_order_parameter_values()->add_order_parameter_values_previous(orderParameterValuesPrevious[i]);
        }
        state->mutable_cme_state()->mutable_order_parameter_values()->add_time(time);
    }

    // Get the histogram values.
    for (size_t i=0; i<tilingHistograms.size(); i++)
    {
        tilingHistograms[i]->serializeInto(state->mutable_cme_state()->add_tiling_histogram());
    }

    // Get the previous barrier crossing counts, if any.
    for (size_t i=0; i<numberTrackingBarriers; i++)
    {
        state->mutable_cme_state()->add_barrier_crossings(trackingBarrierPriorCrossings[i]+trackingBarrierTimesTimes[i].size());
    }
}

void CMESolver::setState(const lm::io::TrajectoryState& state, uint trajectoryNumber)
{
    if (trajectoryNumber >= getSimultaneousTrajectories()) throw lm::InvalidArgException("trajectoryNumber", "exceeded the maximum number of simultaneous trajectories",trajectoryNumber,getSimultaneousTrajectories());

    // Validate the state.
    if (!state.has_cme_state()) throw Exception("State object does not contain the necessary data to initialize the solver.");
    if (state.cme_state().species_counts().number_species() != (int)reactionModel->numberSpecies) throw Exception("State object and reaction model have differing number of species",state.cme_state().species_counts().number_species(),reactionModel->numberSpecies);
    if (state.cme_state().species_counts().number_entries() != 1 || state.cme_state().species_counts().species_count_size() != (int)reactionModel->numberSpecies || state.cme_state().species_counts().time_size() != 1) throw Exception("State object has too many entries",state.cme_state().species_counts().number_entries());
    if (state.cme_state().tiling_histogram().size() != 0 && static_cast<size_t>(state.cme_state().tiling_histogram().size()) != tilingHistograms.size()) THROW_EXCEPTION(lm::RuntimeException, "CME state has an incorrect number of tiling histograms: %d/%d",state.cme_state().tiling_histogram().size(),tilingHistograms.size());

    // Load the trajectory id.
    trajectoryId = state.trajectory_id();

    // Load the previously started flag.
    previouslyStarted = state.trajectory_started();

    // Set the total step count.
    totalSteps = state.cme_state().total_steps();

    // Set the time.
    time = state.cme_state().species_counts().time(0);

    // Set the species counts. This has to happen first since some parts of the state (like the order parameter values) is calculated from the species counts
    for (int i=0; i<state.cme_state().species_counts().species_count_size(); i++)
    {
        speciesCounts[i] = state.cme_state().species_counts().species_count(i);

        // Initialize the previous count, either from the message or from the inital state.
        if (!previouslyStarted && state.cme_state().species_counts().species_count_previous().size() <= i)
            speciesCountsPrevious[i] =  speciesCounts[i];
        else
            speciesCountsPrevious[i] = state.cme_state().species_counts().species_count_previous(i);
    }

    // Set the degree advancements state.
    if (state.cme_state().has_degree_advancements())
    {
        NDArraySerializer::deserializeInto(degreeAdvancements, utuple(reactionModel->numberReactions), state.cme_state().degree_advancements());
        trackDegreeAdvancements = true;
        hasUpdateReactionCountsListeners = true;
    }

    // Set the first passage times.
    numberFptSpecies = state.cme_state().first_passage_times().size();
    if (numberFptSpecies > 0)
    {
        fptSpeciesValues = new lm::me::FPTDeque[numberFptSpecies];
        for (int i=0; i<numberFptSpecies; i++)
        {
            // Load the fpt table.
            fptSpeciesValues[i].deserializeFrom(state.cme_state().first_passage_times(i));

            // Add the initial value to the FPT table, if necessary.
            if (!previouslyStarted) fptSpeciesValues[i].insert(speciesCounts[fptSpeciesValues[i].species], time);
        }
        hasUpdateSpeciesCountsListeners = true;
    }
    numberFptOrderParameters = state.cme_state().order_parameter_first_passage_times().size();
    if (numberFptOrderParameters > 0)
    {
        fptOrderParameterValues = new lm::me::FPTDeque[numberFptOrderParameters];
        for (int i=0; i<numberFptOrderParameters; i++)
        {
            // Load the fpt table.
            //TODOfptOrderParameterValues[i].deserializeFrom(state.cme_state().order_parameter_first_passage_times(i));

            // Add the initial value to the FPT table, if necessary.
            //TODOif (!previouslyStarted) fptOrderParameterValues[i].insert(speciesCounts[fptSpeciesValues[i].species], time);
        }
        hasUpdateSpeciesCountsListeners = true;
    }



    // Set the order parameter first passage times.
//    numberFptTrackedOrderParameters = state.cme_state().order_parameter_first_passage_times_size();
//    if (numberFptTrackedOrderParameters > 0)
//    {
//        fptTrackedOrderParameters = new OParamFPTTracking[numberFptTrackedOrderParameters];
//        for (int i=0; i<numberFptTrackedOrderParameters; i++)
//        {
//            fptTrackedOrderParameters[i].deserializeFrom(state.cme_state().order_parameter_first_passage_times(i));
//        }
//        hasUpdateSpeciesCountsListeners = true;
//    }

    // Set the previous limit reached by the associated trajectory, if any.
//    if (state.has_limit_reached())
//    {
//        limitIDReached = state.limit_reached().id();
//        limitTypeReached = state.limit_reached().limit_type();
//    }

    // Set the order parameter values.
    if (state.cme_state().order_parameter_values().order_parameter_values_size() > 0)
    {
        // set the order parameter values from the passed in state, if any
        for (int i=0; i<state.cme_state().order_parameter_values().order_parameter_values_size(); i++)
        {
            orderParameterValues[i] = state.cme_state().order_parameter_values().order_parameter_values(i);

            // Initialize the previous value, either from the message or from the inital state.
            if (!previouslyStarted && state.cme_state().order_parameter_values().order_parameter_values_previous().size() <= i)
                orderParameterValuesPrevious[i] = orderParameterValues[i];
            else
                orderParameterValuesPrevious[i] = state.cme_state().order_parameter_values().order_parameter_values_previous(i);
        }
    }
    else
    {
        // Calculate the initial order parameters from the species counts.
        for (int i=0; i<numberOrderParameters; i++)
        {
            orderParameterValues[i] = orderParameterFunctions[i]->calculate(time, speciesCounts, reactionModel->numberSpecies);
            orderParameterValuesPrevious[i] = orderParameterValues[i];
        }
    }

    // Set the histogram values.
    for (size_t i=0; i<tilingHistograms.size(); i++)
    {
        // If the histogram was not in the state, intialize to zero.
        if (state.cme_state().tiling_histogram().size() == 0)
            tilingHistograms[i]->clear();

        // Otherwise, initialize from the state.
        else
            tilingHistograms[i]->deserializeFrom(state.cme_state().tiling_histogram(static_cast<int>(i)));
    }

    // Set the previous barrier crossing counts, if any.
    for (size_t i=0; i<numberTrackingBarriers && i<static_cast<size_t>(state.cme_state().barrier_crossings_size()); i++)
    {
        trackingBarrierPriorCrossings[i] = state.cme_state().barrier_crossings(static_cast<int>(i));
    }
}

lm::message::WorkUnitOutput* CMESolver::getOutput(uint trajectoryNumber)
{
    if (trajectoryNumber > 0) throw lm::InvalidArgException("trajectoryNumber", "exceeded the maximum number of simultaneous trajectories",trajectoryNumber,getSimultaneousTrajectories());

    // Get the output pointer.
    lm::message::WorkUnitOutput* ret = output;

    // Forget about the pointer, since the caller is now responsible for it.
    output = NULL;

    return ret;
}

lm::message::WorkUnitStatus::Status CMESolver::getStatus(uint trajectoryNumber)
{
    if (trajectoryNumber > 0) throw lm::InvalidArgException("trajectoryNumber", "exceeded the maximum number of simultaneous trajectories",trajectoryNumber,getSimultaneousTrajectories());
    return status;
}

const lm::types::TrajectoryLimit& CMESolver::getLimitReached(uint trajectoryNumber)
{
    if (trajectoryNumber > 0) throw lm::InvalidArgException("trajectoryNumber", "exceeded the maximum number of simultaneous trajectories",trajectoryNumber,getSimultaneousTrajectories());
    return limitReached;
}

void CMESolver::trajectoryStarted()
{
    //
    // Initialize any time series data.
    //

    // Get the interval for writing species counts.
    if (writeSpeciesTimeSeries)
    {        
        // Initialize the time series, filling in the first set of values if necessary.
        if (!previouslyStarted && writeInitialTrajectoryState)
            speciesTimeSeries.initialize(reactionModel->numberSpecies, speciesWriteInterval, time, speciesCounts, reactionModel->numberSpecies);
        else
            speciesTimeSeries.initialize(reactionModel->numberSpecies, speciesWriteInterval, time);
    }

    // Get the interval for writing degree advancements.
    if (writeDegreeAdvancementTimeSeries && trackDegreeAdvancements)
    {
        // Initialize the time series, filling in the first set of values if necessary.
        if (!previouslyStarted && writeInitialTrajectoryState)
            degreeAdvancementTimeSeries.initialize(reactionModel->numberReactions, degreeAdvancementWriteInterval, time, degreeAdvancements, reactionModel->numberReactions);
        else
            degreeAdvancementTimeSeries.initialize(reactionModel->numberReactions, degreeAdvancementWriteInterval, time);
    }

    // Get the interval for writing order parameters.
    if (writeOrderParameterTimeSeries)
    {
        // Initialize the time series, filling in the first set of values if necessary.
        if (!previouslyStarted && writeInitialTrajectoryState)
            orderParameterTimeSeries.initialize(static_cast<size_t>(numberOrderParameters), orderParameterWriteInterval, time, orderParameterValues, static_cast<size_t>(numberOrderParameters));
        else
            orderParameterTimeSeries.initialize(static_cast<size_t>(numberOrderParameters), orderParameterWriteInterval, time);
    }
}

bool CMESolver::performTimeIncrement(double deltaT)
{
    // Increment the time.
    timeStep = deltaT;
    time += timeStep;

    // See if we are outside of the time limit, and set the time to the max time if so.
    bool outsideTimeLimit = isTrajectoryOutsideTimeLimits();

    // Call the time updated method.
    timeUpdated();

    return outsideTimeLimit;
}

bool CMESolver::isTrajectoryOutsideTimeLimits()
{
    for (int i=0; i<timeLimits.limits().size(); i++)
    {
        const lm::types::TrajectoryLimit& l = timeLimits.limits(i);
        bool wasLimitReached = false;

        switch (l.type())
        {
        case lm::types::TrajectoryLimit::TIME:
            switch (l.stopping_condition())
            {
            case lm::types::TrajectoryLimit::MIN_INCLUSIVE:
                check_limit_MIN_INCLUSIVE(time, l.stopping_value_double(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MIN_EXCLUSIVE:
                check_limit_MIN_EXCLUSIVE(time, l.stopping_value_double(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MAX_INCLUSIVE:
                check_limit_MAX_INCLUSIVE(time, l.stopping_value_double(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MAX_EXCLUSIVE:
                check_limit_MAX_EXCLUSIVE(time, l.stopping_value_double(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::DECREASING_INCLUSIVE: throw Exception("unimplemented");
            case lm::types::TrajectoryLimit::DECREASING_EXCLUSIVE: throw Exception("unimplemented");
            case lm::types::TrajectoryLimit::INCREASING_INCLUSIVE: throw Exception("unimplemented");
            case lm::types::TrajectoryLimit::INCREASING_EXCLUSIVE: throw Exception("unimplemented");
            }
            break;
        default:
            break;
        }

        if (wasLimitReached)
        {
            timeStep -= time-l.stopping_value_double();
            time = l.stopping_value_double();
            status = lm::message::WorkUnitStatus::LIMIT_REACHED;
            limitReached.CopyFrom(l);
            return true;
        }
    }
    return false;
}

void CMESolver::setTrajectoryToMaxTime()
{
    // If we have a max time limit, say that we reached it.
    bool reachedMaxTime = false;
    for (int i=0; i<timeLimits.limits().size(); i++)
    {
        // See if we found a max time limit.
        if (timeLimits.limits(i).type() == lm::types::TrajectoryLimit::TIME && (timeLimits.limits(i).stopping_condition() == lm::types::TrajectoryLimit::MAX_EXCLUSIVE || timeLimits.limits(i).stopping_condition() == lm::types::TrajectoryLimit::MAX_INCLUSIVE))
        {
            // Set the time to the max limit.
            timeStep = timeLimits.limits(i).stopping_value_double()-time;
            time = timeLimits.limits(i).stopping_value_double();

            // We changed the time, so call the update time method.
            timeUpdated();

            // Mark the we found a time limit and stop the loop.
            status = lm::message::WorkUnitStatus::LIMIT_REACHED;
            limitReached.CopyFrom(timeLimits.limits(i));
            reachedMaxTime = true;
            break;
        }
    }

    // Otherwise, this is an error.
    if (!reachedMaxTime)
    {
        status = lm::message::WorkUnitStatus::ERROR;
    }
}

void CMESolver::timeUpdated()
{
    //
    // Record any time series data.
    //

    // If we are writing species time steps, write out any species time steps before this event occurred.
    if (writeSpeciesTimeSeries)
    {
        // Write time steps until the next write time is past the current time.
        speciesTimeSeries.append(speciesCounts, reactionModel->numberSpecies, time);
    }

    // If we are writing degree advancement time steps, write out any degree advancement time steps before this event occurred.
    if (writeDegreeAdvancementTimeSeries && trackDegreeAdvancements)
    {
        // Write degree advancement time steps until the next write time is past the current time.
        degreeAdvancementTimeSeries.append(degreeAdvancements, reactionModel->numberReactions, time);
    }

    // If we are writing order parameter time steps, write out any order parameter time steps before this event occurred.
    if (writeOrderParameterTimeSeries)
    {
        // Write order parameter time steps until the next write time is past the current time.
        orderParameterTimeSeries.append(orderParameterValues, static_cast<size_t>(numberOrderParameters), time);
    }
}

bool CMESolver::performReactionEvent(uint r)
{
    // Update the total step count.
    totalSteps++;

//    Print::printf(Print::DEBUG, "Before reaction %d: prev=%d,%d,%d,%d,%d,%d,%d cur=%d,%d,%d,%d,%d,%d,%d", r, speciesCountsPrevious[0],speciesCountsPrevious[1],speciesCountsPrevious[2],speciesCountsPrevious[3],speciesCountsPrevious[4],speciesCountsPrevious[5],speciesCountsPrevious[6],speciesCounts[0],speciesCounts[1],speciesCounts[2],speciesCounts[3],speciesCounts[4],speciesCounts[5],speciesCounts[6]);

    // Copy the previous species counts.
    memcpy(speciesCountsPrevious, speciesCounts, reactionModel->numberSpecies*sizeof(*speciesCounts));

    // Update the species counts according to the dependency tables.
    for (uint i=0; i<reactionModel->numberDependentSpecies[r]; i++)
    {
        // Get the species index.
        uint s = reactionModel->dependentSpecies[r][i];

        // Update the count.
        speciesCounts[s] += reactionModel->dependentSpeciesChange[r][i];

    }
//    Print::printf(Print::DEBUG, "After reaction %d: prev=%d,%d,%d,%d,%d,%d,%d cur=%d,%d,%d,%d,%d,%d,%d", r, speciesCountsPrevious[0],speciesCountsPrevious[1],speciesCountsPrevious[2],speciesCountsPrevious[3],speciesCountsPrevious[4],speciesCountsPrevious[5],speciesCountsPrevious[6],speciesCounts[0],speciesCounts[1],speciesCounts[2],speciesCounts[3],speciesCounts[4],speciesCounts[5],speciesCounts[6]);

    // Update any order parameters.
    for (int i=0; i<numberOrderParameters; i++)
    {
        orderParameterValuesPrevious[i] = orderParameterValues[i];
        orderParameterValues[i] = orderParameterFunctions[i]->calculate(time, speciesCounts, reactionModel->numberSpecies);
    }

    // If we have any reflecting barriers, see if they were hit.
    if (reflectingBarriers.barrier().size() > 0) applyReflectingBarriers();

    if (hasUpdateReactionCountsListeners) reactionCountsUpdated(r);
    if (hasUpdateSpeciesCountsListeners) speciesCountsUpdated();

    // If we are outside of the limits, stop the trajectory.
    bool outsideStateLimit = isTrajectoryOutsideStateLimits();

    return outsideStateLimit;
}

void CMESolver::applyReflectingBarriers()
{
    // Go through each reflecting barrier, and see if it was crossed.
    bool shouldReflectReaction = false;
    for (int i=0; i<reflectingBarriers.barrier().size(); i++)
    {
        const lm::types::TrajectoryBarrier& b = reflectingBarriers.barrier(i);
        bool crossingOccurred = false;
        switch (b.type())
        {
        case lm::types::TrajectoryBarrier::SPECIES:
            switch (b.behavior())
            {
            case lm::types::TrajectoryBarrier::REFLECTING:
                check_limit_CROSSED_INCLUSIVE(speciesCountsPrevious[b.type_arg()], speciesCounts[b.type_arg()], b.behavior_value_int(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::REFLECTING_DECREASING:
                check_limit_CROSSED_DECREASING_INCLUSIVE(speciesCountsPrevious[b.type_arg()], speciesCounts[b.type_arg()], b.behavior_value_int(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::REFLECTING_INCREASING:
                check_limit_CROSSED_INCREASING_INCLUSIVE(speciesCountsPrevious[b.type_arg()], speciesCounts[b.type_arg()], b.behavior_value_int(), crossingOccurred);
                break;
            default:
                THROW_EXCEPTION(lm::RuntimeException, "Invalid behavior in reflecting barrier: %d", b.behavior());
            }
            break;
        case lm::types::TrajectoryBarrier::ORDER_PARAMETER:
            switch (b.behavior())
            {
            case lm::types::TrajectoryBarrier::REFLECTING:
                check_limit_CROSSED_INCLUSIVE(orderParameterValuesPrevious[b.type_arg()], orderParameterValues[b.type_arg()], b.behavior_value_double(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::REFLECTING_DECREASING:
                check_limit_CROSSED_DECREASING_INCLUSIVE(orderParameterValuesPrevious[b.type_arg()], orderParameterValues[b.type_arg()], b.behavior_value_double(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::REFLECTING_INCREASING:
                check_limit_CROSSED_INCREASING_INCLUSIVE(orderParameterValuesPrevious[b.type_arg()], orderParameterValues[b.type_arg()], b.behavior_value_double(), crossingOccurred);
                break;
            default:
                THROW_EXCEPTION(lm::RuntimeException, "Invalid behavior in reflecting barrier: %d", b.behavior());
            }
        }

        // If this barrier was crossed, mark that we need to reflect the reaction.
        if (crossingOccurred) shouldReflectReaction = true;
    }

    // See if we need to reflect the reaction.
    if (shouldReflectReaction)
    {
        // Reverse the reaction.
        memcpy(speciesCounts, speciesCountsPrevious, sizeof(*speciesCounts)*reactionModel->numberSpecies);
        memcpy(orderParameterValues, orderParameterValuesPrevious, sizeof(*orderParameterValues)*static_cast<size_t>(numberOrderParameters));
    }
}

bool CMESolver::isTrajectoryOutsideStateLimits()
{
    for (int i=0; i<stateLimits.limits().size(); i++)
    {
        const lm::types::TrajectoryLimit& l = stateLimits.limits(i);
        bool wasLimitReached = false;
        size_t barrierIndex;
        int barrierCrossings;

        switch (l.type())
        {
        case lm::types::TrajectoryLimit::SPECIES:
            switch (l.stopping_condition())
            {
            case lm::types::TrajectoryLimit::MIN_INCLUSIVE:
                check_limit_MIN_INCLUSIVE(speciesCounts[l.type_arg()], l.stopping_value_int(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MIN_EXCLUSIVE:
                check_limit_MIN_EXCLUSIVE(speciesCounts[l.type_arg()], l.stopping_value_int(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MAX_INCLUSIVE:
                check_limit_MAX_INCLUSIVE(speciesCounts[l.type_arg()], l.stopping_value_int(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MAX_EXCLUSIVE:
                check_limit_MAX_EXCLUSIVE(speciesCounts[l.type_arg()], l.stopping_value_int(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::DECREASING_INCLUSIVE:
                check_limit_CROSSED_DECREASING_INCLUSIVE(speciesCountsPrevious[l.type_arg()], speciesCounts[l.type_arg()], l.stopping_value_int(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::DECREASING_EXCLUSIVE:
                check_limit_CROSSED_DECREASING_EXCLUSIVE(speciesCountsPrevious[l.type_arg()], speciesCounts[l.type_arg()], l.stopping_value_int(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::INCREASING_INCLUSIVE:
                check_limit_CROSSED_INCREASING_INCLUSIVE(speciesCountsPrevious[l.type_arg()], speciesCounts[l.type_arg()], l.stopping_value_int(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::INCREASING_EXCLUSIVE:
                check_limit_CROSSED_INCREASING_EXCLUSIVE(speciesCountsPrevious[l.type_arg()], speciesCounts[l.type_arg()], l.stopping_value_int(), wasLimitReached);
                break;
            }
            break;

        case lm::types::TrajectoryLimit::ORDER_PARAMETER:
            switch (l.stopping_condition())
            {
            case lm::types::TrajectoryLimit::MIN_INCLUSIVE:
                check_limit_MIN_INCLUSIVE(orderParameterValues[l.type_arg()], l.stopping_value_double(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MIN_EXCLUSIVE:
                check_limit_MIN_EXCLUSIVE(orderParameterValues[l.type_arg()], l.stopping_value_double(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MAX_INCLUSIVE:
                check_limit_MAX_INCLUSIVE(orderParameterValues[l.type_arg()], l.stopping_value_double(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MAX_EXCLUSIVE:
                check_limit_MAX_EXCLUSIVE(orderParameterValues[l.type_arg()], l.stopping_value_double(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::DECREASING_INCLUSIVE:
                check_limit_CROSSED_DECREASING_INCLUSIVE(orderParameterValuesPrevious[l.type_arg()], orderParameterValues[l.type_arg()], l.stopping_value_double(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::DECREASING_EXCLUSIVE:
                check_limit_CROSSED_DECREASING_EXCLUSIVE(orderParameterValuesPrevious[l.type_arg()], orderParameterValues[l.type_arg()], l.stopping_value_double(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::INCREASING_INCLUSIVE:
                check_limit_CROSSED_INCREASING_INCLUSIVE(orderParameterValuesPrevious[l.type_arg()], orderParameterValues[l.type_arg()], l.stopping_value_double(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::INCREASING_EXCLUSIVE:
                check_limit_CROSSED_INCREASING_EXCLUSIVE(orderParameterValuesPrevious[l.type_arg()], orderParameterValues[l.type_arg()], l.stopping_value_double(), wasLimitReached);
                break;
            }
            break;

        case lm::types::TrajectoryLimit::DEGREE_ADVANCEMENT:
            switch (l.stopping_condition())
            {
            case lm::types::TrajectoryLimit::MIN_INCLUSIVE:
                check_limit_MIN_INCLUSIVE(degreeAdvancements[l.type_arg()], uint64_t(l.stopping_value_int()), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MIN_EXCLUSIVE:
                check_limit_MIN_EXCLUSIVE(degreeAdvancements[l.type_arg()], uint64_t(l.stopping_value_int()), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MAX_INCLUSIVE:
                check_limit_MAX_INCLUSIVE(degreeAdvancements[l.type_arg()], uint64_t(l.stopping_value_int()), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MAX_EXCLUSIVE:
                check_limit_MAX_EXCLUSIVE(degreeAdvancements[l.type_arg()], uint64_t(l.stopping_value_int()), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::DECREASING_INCLUSIVE: THROW_EXCEPTION(lm::RuntimeException, "unimplemented");
            case lm::types::TrajectoryLimit::DECREASING_EXCLUSIVE: THROW_EXCEPTION(lm::RuntimeException, "unimplemented");
            case lm::types::TrajectoryLimit::INCREASING_INCLUSIVE: THROW_EXCEPTION(lm::RuntimeException, "unimplemented");
            case lm::types::TrajectoryLimit::INCREASING_EXCLUSIVE: THROW_EXCEPTION(lm::RuntimeException, "unimplemented");
            }
            break;

        case lm::types::TrajectoryLimit::BARRIER_CROSSING:
            barrierIndex = l.type_arg()-static_cast<uint32_t>(reflectingBarriers.barrier().size());
            if (barrierIndex >= numberTrackingBarriers) THROW_EXCEPTION(lm::RuntimeException, "Invalid tracking barrier index: %d,%d", barrierIndex, numberTrackingBarriers);
            barrierCrossings = static_cast<int32_t>(trackingBarrierTimesTimes[barrierIndex].size() + trackingBarrierPriorCrossings[barrierIndex]);
            switch (l.stopping_condition())
            {
            case lm::types::TrajectoryLimit::MIN_INCLUSIVE: THROW_EXCEPTION(lm::RuntimeException, "unimplemented");
            case lm::types::TrajectoryLimit::MIN_EXCLUSIVE: THROW_EXCEPTION(lm::RuntimeException, "unimplemented");
            case lm::types::TrajectoryLimit::MAX_INCLUSIVE:
                check_limit_MAX_INCLUSIVE(barrierCrossings, l.stopping_value_int(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::MAX_EXCLUSIVE:
                check_limit_MAX_EXCLUSIVE(barrierCrossings, l.stopping_value_int(), wasLimitReached);
                break;
            case lm::types::TrajectoryLimit::DECREASING_INCLUSIVE: THROW_EXCEPTION(lm::RuntimeException, "unimplemented");
            case lm::types::TrajectoryLimit::DECREASING_EXCLUSIVE: THROW_EXCEPTION(lm::RuntimeException, "unimplemented");
            case lm::types::TrajectoryLimit::INCREASING_INCLUSIVE: THROW_EXCEPTION(lm::RuntimeException, "unimplemented");
            case lm::types::TrajectoryLimit::INCREASING_EXCLUSIVE: THROW_EXCEPTION(lm::RuntimeException, "unimplemented");
            }
            break;

        default:
            break;
        }

        if (wasLimitReached)
        {
            Print::printf(Print::VERBOSE_DEBUG, "Reached limit %d (%d,%d) at %0.6e: %d -> %d", i, l.type(), l.stopping_value_int(), time, speciesCountsPrevious[0], speciesCounts[0]);
            status = lm::message::WorkUnitStatus::LIMIT_REACHED;
            limitReached.CopyFrom(l);
            return true;
        }
    }
    return false;
}

void CMESolver::reactionCountsUpdated(uint r)
{
    // Update the degree advancement, if enabled
    if (trackDegreeAdvancements)
    {
        degreeAdvancements[r]++;
    }
}

void CMESolver::speciesCountsUpdated()
{
    // Update the first passage time tables.
    for (int i=0; i<numberFptSpecies; i++)
    {
        int value = speciesCounts[fptSpeciesValues[i].species];
        if (value < fptSpeciesValues[i].minValue || value > fptSpeciesValues[i].maxValue)
            fptSpeciesValues[i].insert(value, time);
    }

    // Update the order parameter fpt tables.
    for (int i=0; i<numberFptOrderParameters; i++)
    {
        //TODO
//        int value = speciesCounts[fptSpeciesValues[i].species];
//        if (value < fptSpeciesValues[i].minValue || value > fptSpeciesValues[i].maxValue)
//            fptSpeciesValues[i].insert(value, time);
    }

    // Update the barrier crossing tables.
    for (size_t i=0; i<numberTrackingBarriers; i++)
    {
        bool crossingOccurred = false;
        const lm::types::TrajectoryBarrier& b = trackingBarriers.barrier(static_cast<int>(i));
        switch (b.type())
        {
        case lm::types::TrajectoryBarrier::SPECIES:
            switch (b.behavior())
            {
            case lm::types::TrajectoryBarrier::TRACKING:
                check_limit_CROSSED_INCLUSIVE(speciesCountsPrevious[b.type_arg()], speciesCounts[b.type_arg()], b.behavior_value_int(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::TRACKING_DECREASING_INCLUSIVE:
                check_limit_CROSSED_DECREASING_INCLUSIVE(speciesCountsPrevious[b.type_arg()], speciesCounts[b.type_arg()], b.behavior_value_int(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::TRACKING_DECREASING_EXCLUSIVE:
                check_limit_CROSSED_DECREASING_EXCLUSIVE(speciesCountsPrevious[b.type_arg()], speciesCounts[b.type_arg()], b.behavior_value_int(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::TRACKING_INCREASING_INCLUSIVE:
                check_limit_CROSSED_INCREASING_INCLUSIVE(speciesCountsPrevious[b.type_arg()], speciesCounts[b.type_arg()], b.behavior_value_int(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::TRACKING_INCREASING_EXCLUSIVE:
                check_limit_CROSSED_INCREASING_EXCLUSIVE(speciesCountsPrevious[b.type_arg()], speciesCounts[b.type_arg()], b.behavior_value_int(), crossingOccurred);
                break;
            default:
                THROW_EXCEPTION(lm::RuntimeException, "Invalid behavior in tracking barrier: %d", b.behavior());
            }
            break;
        case lm::types::TrajectoryBarrier::ORDER_PARAMETER:
            switch (b.behavior())
            {
            case lm::types::TrajectoryBarrier::TRACKING:
                check_limit_CROSSED_INCLUSIVE(orderParameterValuesPrevious[b.type_arg()], orderParameterValues[b.type_arg()], b.behavior_value_double(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::TRACKING_DECREASING_INCLUSIVE:
                check_limit_CROSSED_DECREASING_INCLUSIVE(orderParameterValuesPrevious[b.type_arg()], orderParameterValues[b.type_arg()], b.behavior_value_double(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::TRACKING_DECREASING_EXCLUSIVE:
                check_limit_CROSSED_DECREASING_EXCLUSIVE(orderParameterValuesPrevious[b.type_arg()], orderParameterValues[b.type_arg()], b.behavior_value_double(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::TRACKING_INCREASING_INCLUSIVE:
                check_limit_CROSSED_INCREASING_INCLUSIVE(orderParameterValuesPrevious[b.type_arg()], orderParameterValues[b.type_arg()], b.behavior_value_double(), crossingOccurred);
                break;
            case lm::types::TrajectoryBarrier::TRACKING_INCREASING_EXCLUSIVE:
                check_limit_CROSSED_INCREASING_EXCLUSIVE(orderParameterValuesPrevious[b.type_arg()], orderParameterValues[b.type_arg()], b.behavior_value_double(), crossingOccurred);
                break;
            default:
                THROW_EXCEPTION(lm::RuntimeException, "Invalid behavior in tracking barrier: %d", b.behavior());
            }
        }

        // If there was a crossing, add it to the list.
        if (crossingOccurred)
        {
            if (b.type() == lm::types::TrajectoryBarrier::ORDER_PARAMETER)
                Print::printf(Print::VERBOSE_DEBUG, "Crossed tracking barrier %d (%d) at %0.6e %ld: %0.3f -> %0.3f", i, b.behavior(), time, totalSteps, orderParameterValuesPrevious[0], orderParameterValues[0]);
            else
                Print::printf(Print::VERBOSE_DEBUG, "Crossed tracking barrier %d (%d) at %0.6e %ld: %d -> %d", i, b.behavior(), time, totalSteps, speciesCountsPrevious[0], speciesCounts[0]);
            // Record the time and species counts for the crossing event.
            for (uint j=0; j<reactionModel->numberSpecies; j++) trackingBarrierTimesCounts[i].push_back(speciesCounts[j]);
            trackingBarrierTimesTimes[i].push_back(time);
            trackingBarrierTimesSteps[i].push_back(totalSteps);
        }
    }

//        // Update any tilingHists.
//        if (tilings != NULL)
//        {
//            for (int i=0;i<numberTilingHists;i++)
//            {
//                tilingHists[i].tileVals[(*tilings)[tilingHists[i].tilingID]->getTileIndex((*oparams)[(*tilings)[tilingHists[i].tilingID]->getOrderParameterID()]->get())] += timeStep;
//            }
//        }
}

void CMESolver::trajectoryFinished()
{
    //
    // Record any final time series data.
    //

    // If we hit any limit, write out the final time if requested.
    if (status == lm::message::WorkUnitStatus::LIMIT_REACHED && writeFinalTrajectoryState)
    {
        // write out the last species count, or ensure that it already has been
        if (writeSpeciesTimeSeries)
        {
            speciesTimeSeries.appendFinal(speciesCounts, reactionModel->numberSpecies, time);
        }

        // write out the last degree advancement count, or ensure that it already has been
        if (writeDegreeAdvancementTimeSeries && trackDegreeAdvancements)
        {
            degreeAdvancementTimeSeries.appendFinal(degreeAdvancements, reactionModel->numberReactions, time);
        }

        // write out the last order parameter value, or ensure that it already has been
        if (writeOrderParameterTimeSeries)
        {
            orderParameterTimeSeries.appendFinal(orderParameterValues, static_cast<size_t>(numberOrderParameters), time);
        }
    }

    //
    // Add any saved records to the output.
    //

    // If we have any species time series data, add them to the output message.
    if (speciesTimeSeries.size() > 0)
    {
        // Mark that the message does contain some data.
        output->set_has_output(true);

        // Serialize the data into the message.
        lm::io::SpeciesTimeSeries* speciesTimeSeriesDataSet = output->mutable_species_time_series();
        speciesTimeSeriesDataSet->set_trajectory_id(trajectoryId);
        speciesTimeSeries.serializeInto(speciesTimeSeriesDataSet->mutable_counts(), speciesTimeSeriesDataSet->mutable_times());
    }

    // If we have any degree advancement time series data, add them to the output message.
    if (degreeAdvancementTimeSeries.size() > 0)
    {
        // Mark that the message does contain some data.
        output->set_has_output(true);

        // Serialize the data into the message.
        lm::io::DegreeAdvancementTimeSeries* degreeAdvancementTimeSeriesDataSet = output->mutable_degree_advancement_time_series();
        degreeAdvancementTimeSeriesDataSet->set_trajectory_id(trajectoryId);
        degreeAdvancementTimeSeries.serializeInto(degreeAdvancementTimeSeriesDataSet->mutable_counts(), degreeAdvancementTimeSeriesDataSet->mutable_times());
    }

    // If we have any order parameter time series data, add them to the output message.
    if (orderParameterTimeSeries.size())
    {
        // Mark that the message does contain some data.
        output->set_has_output(true);

        // Serialize the data into the message.
        lm::io::OrderParameterTimeSeries* orderParameterTimeSeriesDataSet = output->mutable_order_parameter_time_series();
        orderParameterTimeSeriesDataSet->set_trajectory_id(trajectoryId);
        orderParameterTimeSeries.serializeInto(orderParameterTimeSeriesDataSet->mutable_values(), orderParameterTimeSeriesDataSet->mutable_times());
    }

    // If we have any barrier tracking data, add it to the output message.
    for (size_t i=0; i<numberTrackingBarriers; i++)
    {
        // Make sure the arrays are of a consistent size.
        if (trackingBarrierTimesCounts[i].size() != trackingBarrierTimesTimes[i].size()*reactionModel->numberSpecies || trackingBarrierTimesTimes[i].size() != trackingBarrierTimesSteps[i].size()) THROW_EXCEPTION(lm::RuntimeException, "Tracking barrier %d species counts and time mismatch.", i);

        // See if we have any barrier tracking data to write.
        if (trackingBarrierTimesTimes[i].size() > 0)
        {
            // Mark that the message does contain some data.
            output->set_has_output(true);

            // Create the record.
            lm::io::BarrierCrossingTimes* bct = output->add_barrier_crossing_times();
            bct->set_trajectory_id(trajectoryId);
            bct->set_barrier_index(static_cast<uint32_t>(i));

            // Serialize the times.
            utuple shape(uint(trackingBarrierTimesTimes[i].size()));
            robertslab::pbuf::NDArraySerializer::serializeInto<double>(bct->mutable_times(), trackingBarrierTimesTimes[i].data(), shape);

            // Serialize the species counts.
            robertslab::pbuf::NDArraySerializer::serializeInto<int32_t>(bct->mutable_counts(), trackingBarrierTimesCounts[i].data(), utuple(uint(trackingBarrierTimesTimes[i].size()),reactionModel->numberSpecies));

            // Serialize the total steps.
            robertslab::pbuf::NDArraySerializer::serializeInto<uint64_t>(bct->mutable_total_steps(), trackingBarrierTimesSteps[i].data(), utuple(uint(trackingBarrierTimesTimes[i].size())));
        }
    }

    // If the simulation reached a limit and we are tracking first passage times, add them to the output message.
    if (status == lm::message::WorkUnitStatus::LIMIT_REACHED && numberFptSpecies > 0)
    {
        // Mark that the message does contain some data.
        output->set_has_output(true);

        for (int i=0; i<numberFptSpecies; i++)
        {
            fptSpeciesValues[i].serializeInto(output->add_first_passage_times());
        }
    }

    // If the simulation reached a limit and we are tracking order parameter first passage times, add them to the output message.
    if (status == lm::message::WorkUnitStatus::LIMIT_REACHED && numberFptOrderParameters > 0)
    {
        // Mark that the message does contain some data.
        output->set_has_output(true);

        for (int i=0; i<numberFptOrderParameters; i++)
        {
            //TODOfptOrderParameterValues[i].serializeTo(output->add_order_parameter_first_passage_times(), trajectoryId);
        }
    }
}

}
}
