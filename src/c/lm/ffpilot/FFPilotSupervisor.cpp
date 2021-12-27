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

////#include <algorithm>
////#include <iomanip>
////#include <istream>
////#include <iterator>
////#include <limits>
////#include <map>
#include <sstream>
#include <string>
////#include <valarray>
////#include <vector>

#include "hrtime.h"
#include "lm/ClassFactory.h"
#include "lm/Math.h"
#include "lm/Types.h"
#include "lm/Print.h"
#include "lm/cme/CMETrajectory.h"
#include "lm/ffpilot/TrajectoryData.h"
#include "lm/ffpilot/FFPilotMath.h"
#include "lm/ffpilot/FFPilotSupervisor.h"
#include "lm/ffpilot/FFPilotTrajectoryList.h"
#include "lm/ffpilot/FFPilotData.h"
#include "lm/ffpilot/PhaseData.h"
#include "lm/ffpilot/StageData.h"
#include "lm/io/BarrierCrossingTimes.pb.h"
#include "lm/io/ffpilot/FFPilotOutput.pb.h"
#include "lm/oparam/OrderParameterFunction.h"
////#include "lm/ffpilot/input/FFPilotInput.h"
////#include "lm/ffpilot/input/FFPilotStage.pb.h"
////#include "lm/ffpilot/input/FFPilotPhaseLimit.pb.h"
////#include "lm/ffpilot/FFPilotMath.h"
////#include "lm/input/Tilings.pb.h"
////#include "lm/io/OutputWriter.h"
////#include "lm/io/TrajectoryState.pb.h"
////#include "lm/main/Globals.h"
////#include "lm/simulation/SimulationSupervisor.h"
#include "lm/main/Globals.h"
#include "lm/message/Message.pb.h"
#include "lm/message/ProcessAggregatedOutput.pb.h"
////#include "lm/message/RunWorkUnit.pb.h"
////#include "lm/message/StartedOutputWriter.pb.h"
////#include "lm/message/StartedWorkUnit.pb.h"
////#include "lm/message/WorkUnitOutput.pb.h"
////#include "lm/protowrap/Repeated.h"
////#include "lm/resource/ResourceMap.h"
#include "lm/rng/RandomGenerator.h"
#include "lm/rng/XORShift.h"
#include "lm/types/Tilings.pb.h"
#include "lm/types/TrajectoryLimits.pb.h"
////#include "lm/VectorMath.h"
#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"
#include "robertslab/pbuf/NDArraySerializer.h"

//#include "lm/replicates/ReplicateTrajectoryList.h" //TODO

////using lm::protowrap::Repeated;
////using lm::resource::ResourceMap;
////using std::map;
////using std::setfill;
////using std::setw;
using std::string;
using std::stringstream;
////using std::valarray;
////using std::vector;
using robertslab::pbuf::NDArraySerializer;

namespace lm {
namespace ffpilot {

bool FFPilotSupervisor::registered=FFPilotSupervisor::registerClass();

bool FFPilotSupervisor::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::simulation::SimulationSupervisor","lm::ffpilot::FFPilotSupervisor",&FFPilotSupervisor::allocateObject);
    return true;
}

void* FFPilotSupervisor::allocateObject()
{
    return new FFPilotSupervisor();
}

////// if >0, we use a hand-rolled mpi receive polling scheme in order to reduce the supervisor cpu%
////int FFPilotSupervisor::getRecvSleepMilliseconds()
////{
////    return -1;
////}

FFPilotSupervisor::FFPilotSupervisor()
:random(NULL),orderParameterIndex(-1),orderParameterFunction(NULL),tilingIndex(-1),tilingEdges(NULL),
nextTrajectoryId(0),trajectoryList(NULL)
////:ffpilotPhaseOutputListsWrap(&ffpilotPhaseOutputListsMsg),previousFFPilotPhaseOutputWrapPtr(&_ffpilotPhaseOutputWrap_0),currentFFPilotPhaseOutputWrapPtr(&_ffpilotPhaseOutputWrap_1),ffpilotStageOutputsWrap(&ffpilotStageOutputsMsg),
//// ffpilotPhaseTerminated(false),simulationPhaseOutputSent(false),simulationStageOutputSent(false),input(NULL),trajectoryList(NULL),ffpilotProgress_lastPrintTime(getHrTime())
{
    random = new lm::rng::XORShift(0,0);

////    ffpilotPhaseOutputContainingMsg.mutable_process_work_unit_output()->set_work_unit_id(0);
////    ffpilotPhaseOutputContainingMsg.mutable_process_work_unit_output()->add_part_output();

////    ffpilotStageOutputRawContainingMsg.mutable_process_work_unit_output()->set_work_unit_id(0);
////    ffpilotStageOutputRawContainingMsg.mutable_process_work_unit_output()->add_part_output();

////    ffpilotStageOutputSummaryContainingMsg.mutable_process_work_unit_output()->set_work_unit_id(0);
////    ffpilotStageOutputSummaryContainingMsg.mutable_process_work_unit_output()->add_part_output();
}

FFPilotSupervisor::~FFPilotSupervisor()
{
    if (tilingEdges != NULL) delete tilingEdges; tilingEdges = NULL;
    if (trajectoryList != NULL) delete trajectoryList; trajectoryList = NULL;
    if (random != NULL) delete random; random = NULL;
}

std::string FFPilotSupervisor::getClassName()
{
    return "lm::ffpilot::FFPilotSupervisor";
}

void FFPilotSupervisor::setInput(lm::input::Input* input)
{
    // Call the base class method.
    StagedSimulationSupervisor::setInput(input);

    // Make sure we have the necessary ffpilot input.
    if (!input->has_reaction_model()) THROW_EXCEPTION(lm::InputException, "No reaction model was present in the input.");
    if (!input->has_order_parameters() || input->order_parameters().order_parameter().size() < 1) THROW_EXCEPTION(lm::InputException, "No order parameters were present in the input.");
    if (!input->has_tilings() || input->tilings().tiling().size() < 1) THROW_EXCEPTION(lm::InputException, "No tilings were present in the input.");
    if (!input->has_ffpilot_options()) THROW_EXCEPTION(lm::InputException, "No FFPilot option were present in the input.");

    // Get the order parameter to use.
    if (input->ffpilot_options().order_parameter_index() >= input->order_parameters().order_parameter_size()) THROW_EXCEPTION(lm::InputException, "The order parameter specified for ffpilot was not present in the input.");
    orderParameterIndex = input->ffpilot_options().order_parameter_index();

    // Create the order parameter function.
    lm::oparam::OrderParameterFunctionFactory fs;
    orderParameterFunction = fs.createOrderParameterFunction(input->order_parameters().order_parameter(orderParameterIndex));
    Print::printf(Print::DEBUG, "FFPilot created order parameter function object: %d", orderParameterFunction->getType());

    // Get the tilings to use.
    if (input->ffpilot_options().tiling_index() >= input->tilings().tiling_size()) THROW_EXCEPTION(lm::InputException, "The tiling specified for ffpilot was not present in the input.");
    tilingIndex = input->ffpilot_options().tiling_index();
    tilingEdges = robertslab::pbuf::NDArraySerializer::deserializeAllocate<double>(input->tilings().tiling(tilingIndex).edges());
    ndarray<int32_t>opIndices = robertslab::pbuf::NDArraySerializer::deserialize<int32_t>(input->tilings().tiling(tilingIndex).order_parameter_indices());
    if (opIndices.shape.len != 1 || opIndices.shape[0] != tilingEdges->shape.len) THROW_EXCEPTION(lm::InputException, "The tiling specified had inconsistent sizes.");
    if (opIndices.shape[0] != 1) THROW_EXCEPTION(lm::InputException, "Multidimensional tilings are not currently supported by ffpilot.");
    if (orderParameterIndex != opIndices[0]) THROW_EXCEPTION(lm::InputException, "The order parameter in the specified tiling did not match the ffpilot order parameter.");

    // Check the basins to use when starting phase 0 trajectories.
    ndarray<int32_t> phase0Basins = robertslab::pbuf::NDArraySerializer::deserialize<int32_t>(input->ffpilot_options().phase_zero_basins());
    if (phase0Basins.shape.len != 2 || phase0Basins.shape[0] != 2 || phase0Basins.shape[1] != input->reaction_model().number_species()) THROW_EXCEPTION(lm::InputException, "The basins array was the wrong shape: %d.",phase0Basins.shape.len);

    // See which mode we are running in.
    if (input->ffpilot_options().pilot_skip())
    {
        // We are skipping the pilot phase.
        data.initialize(2);
        data.getStageData(0).initialize(tilingEdges->shape[0], false, true);
        data.getStageData(1).initialize(tilingEdges->shape[0], false, false);

    }
    else if (input->ffpilot_options().prod_skip())
    {
        // We are skipping the production phase.
        data.initialize(2);
        data.getStageData(0).initialize(tilingEdges->shape[0], true, true);
        data.getStageData(1).initialize(tilingEdges->shape[0], true, false);

    }
    else
    {
        // Initialize the ffpilot data structure with pilot and production phases.
        data.initialize(4);
        data.getStageData(0).initialize(tilingEdges->shape[0], true, true, 2);
        data.getStageData(1).initialize(tilingEdges->shape[0], true, false, 3);
        data.getStageData(2).initialize(tilingEdges->shape[0], false, true);
        data.getStageData(3).initialize(tilingEdges->shape[0], false, false);
    }

    Print::printf(Print::DEBUG, "FFPilot supervisor configuration complete: %d stages with %d phases each", data.getNumberStages(), tilingEdges->shape[0]);
}

int FFPilotSupervisor::getNumberStages()
{
    return static_cast<int>(data.getNumberStages());
}

int FFPilotSupervisor::getNumberPhases(int stageIndex)
{
    if (data.getNumberStages() == 0) THROW_EXCEPTION(lm::RuntimeException, "FFPilotSupervisor::getNumberPhases called before setup was completed");
    return static_cast<int>(data.getStageData(static_cast<int>(stageIndex)).getNumberPhases());
}

void FFPilotSupervisor::startSimulation()
{
    Print::printf(Print::INFO, "FFPilot supervisor started.");

    // Call the base class method.
    StagedSimulationSupervisor::startSimulation();

////    if (input->hasFFPilotStages())
////    {
////        initSimulationStageListCustom();
////    }
////    else
////    {
////        initSimulationStageList();
////    }


////    // print out the stage plans for debugging
////    //for (FFPilotStageRepeated::const_iterator it = ffpilotStages.begin(); it != ffpilotStages.end(); it++) std::cout << it->DebugString();

}

////void FFPilotSupervisor::initSimulationStageList()
////{
////    // sanity check the input (make we sure we have at least one tiling, etc)
////    input->initSanityCheck();

////    // build the stage list
////    for (vector<uint64_t>::const_iterator replicateIt = replicates.begin();replicateIt!=replicates.end();++replicateIt)
////    {
////        for (lm::tiling::Tilings::const_iterator tilingIt=input->getTilings().begin();tilingIt!=input->getTilings().end();++tilingIt)
////        {
////            for (int basinIndex=0;basinIndex<tilingIt->second->basins().size();basinIndex++)
////            {
////                // initialize a production stage (and possibly also its pilot stage)
////                FFPilotStageMsg productionStage;
////                buildProductionStage(&productionStage, *replicateIt, *tilingIt->second, basinIndex);

////                // place the stage in the to-be-executed list (the pilot stage, if any, will be placed before the production stage)
////                FFPilotStageMsg* productionStageStored = ffpilotStages.Add();
////                productionStageStored->CopyFrom(productionStage);
////            }
////        }
////    }
////    currentFFPilotStageIter = ffpilotStages.begin();
////    _currentStageIndex = 0;
////}

////void FFPilotSupervisor::initSimulationStageListCustom()
////{
////    for (vector<uint64_t>::const_iterator replicateIt = replicates.begin();replicateIt!=replicates.end();++replicateIt)
////    {
////        // iterate over the stages in input->getFFPilotStages()
////        for (FFPilotStageRepeated::const_iterator stageIt=input->getFFPilotStages().begin(); stageIt!=input->getFFPilotStages().end(); stageIt++)
////        {
////            // copy the stage from the input
////            FFPilotStageMsg ffpilotStage(*stageIt);

////            ffpilotStage.set_replicate_id(*replicateIt);
////            // TODO: implement pilot stages in initSimulationStageListCustom
////            //bool pilotRunRequired = false;

////            // if the stage doesn't already have a tiling set, use the one from the simulation input
////            if (not ffpilotStage.has_tiling())
////            {
////                addTiling(&ffpilotStage, input->getTilings().at(ffpilotStage.tiling_id()), ffpilotStage.basin_id());
////            }

////            // iterate over the phases in ffpilotStage.ffpilot_phases()
////            for (FFPilotPhasesWrap::iterator phaseIt=ffpilotStage.mutable_ffpilot_phases()->begin(); phaseIt!=ffpilotStage.mutable_ffpilot_phases()->end(); phaseIt++)
////            {
////                if (not phaseIt->has_output_options())
////                {
////                    // if the phase doesn't already have output options set, set them in the standard way based on the simulation input
////                    addOptions(&*phaseIt, *stageIt);
////                }

////                if (not phaseIt->has_ffpilot_phase_limit())
////                {
////                    // If the phase doesn't already have a limit set, add it
////                    addFFPilotPhaseLimit(&*phaseIt, FFPhaseLimEnums::TRAJECTORY_COUNT, input->productionStageCountMinimum());
////                }
////                else if (not phaseIt->ffpilot_phase_limit().has_events_per_trajectory())
////                {
////                    // If the phase has a limit but the number of runs/runners is not set, automatically figure it out
////                    buildFFPilotPhaseLimitTrajectoriesToRun(phaseIt->mutable_ffpilot_phase_limit(), *phaseIt, slots.getSimultaneousWorkUnits());
////                }
////            }

////            // place the stage in the to-be-executed list (the pilot stage, if any, will be placed before the production stage)
////            FFPilotStageMsg* stageStored = ffpilotStages.Add();
////            stageStored->CopyFrom(ffpilotStage);
////        }
////    }
////    currentFFPilotStageIter = ffpilotStages.begin();
////    _currentStageIndex = 0;
////}

////FFPilotSupervisor::FFPilotStageMsg* FFPilotSupervisor::buildProductionStage(FFPilotStageMsg* productionStage, uint64_t replicateID, const lm::tiling::Tiling& tiling, int64_t basinIndex)
////{
////    productionStage->set_replicate_id(replicateID);
////    productionStage->set_name("Production");

////    productionStage->set_is_pilot_stage(false);
////    productionStage->set_needs_pilot_stage(true);

////    addTiling(productionStage, tiling, basinIndex);

////    if (input->hasErrorGoal())
////    {
////        // initialize the pilot stage
////        FFPilotStageMsg pilotStage;
////        buildPilotStage(&pilotStage, productionStage);

////        // add the pilot stage to the execution order
////        FFPilotStageMsg* pilotStageStored = ffpilotStages.Add();
////        pilotStageStored->CopyFrom(pilotStage);
////    }

////    addFFPilotPhases(productionStage, FFPhaseEnums::LAZY, FFPhaseEnums::UNIFORM_RANDOM);

////    if (not productionStage->needs_pilot_stage())
////    {
////        addFFPilotPhaseLimitsFromInput(productionStage);
////    }

////    return productionStage;
////}

////void FFPilotSupervisor::addTiling(FFPilotStageMsg* stage, const lm::tiling::Tiling& tiling, int64_t basinIndex)
////{
////    // TODO: encapsulate this mess in a FFPilotStage wrapper
////    // shallow copy the tiling wrapper. We're done with the passed in tiling wrapper
////    lm::tiling::Tiling tilingWrapCopy = tiling;

////    // copy the actual tiling message over to the production stage message
////    stage->mutable_tiling()->CopyFrom(tilingWrapCopy.getTilingMsg());

////    // reseat the tiling wrapper copy around the tiling message copy
////    tilingWrapCopy.setTilingMsg(stage->mutable_tiling());

////    // use the tiling wrapper copy to set the appropriate basin_id in the tiling. This will also reverse the tiling, if needed
////    tilingWrapCopy.setBasin(basinIndex);

////    stage->set_basin_id(basinIndex);
////    stage->set_tiling_id(tilingWrapCopy.id());
////}

////FFPilotSupervisor::FFPilotStageMsg* FFPilotSupervisor::buildPilotStage(FFPilotStageMsg* pilotStage, FFPilotStageMsg* productionStage)
////{
////    pilotStage->set_name("Pilot");

////    pilotStage->set_is_pilot_stage(true);
////    pilotStage->set_needs_pilot_stage(false);

////    pilotStage->set_replicate_id(productionStage->replicate_id());
////    pilotStage->mutable_tiling()->CopyFrom(productionStage->tiling());
////    pilotStage->set_basin_id(productionStage->basin_id());

////    addFFPilotPhases(pilotStage, FFPhaseEnums::LAZY, FFPhaseEnums::UNIFORM_RANDOM);

////    addFFPilotPhaseLimitsForPilotStage(pilotStage, FFPhaseLimEnums::FORWARD_FLUXES, static_cast<uint64_t>(
////        round(input->getFFPilotOptionsMsg().pilot_stage_count()*input->phaseZeroSamplingMultiplier())),
////        input->getFFPilotOptionsMsg().pilot_stage_count()
////    );

////    return pilotStage;
////}

////void FFPilotSupervisor::addFFPilotPhases(FFPilotStageMsg* stage, FFPhaseEnums::TrajectoryGeneration trajGeneration, FFPhaseEnums::TrajectoryDuplication trajDuplication)
////{
////    for (int i=0;i<stage->tiling().edges_size();i++)
////    {
////        FFPilotPhaseMsg* ffpilotPhase = stage->add_ffpilot_phases();

////        ffpilotPhase->set_phase_id(i);
////        ffpilotPhase->set_basin_id(stage->basin_id());
////        ffpilotPhase->set_tiling_id(stage->tiling().id());
////        ffpilotPhase->set_tile_id(i);

////        // set ffpilotPhase values that depend on whether phaseID==0 or phaseID > 0
////        if (i==0)
////        {
////            ffpilotPhase->set_trajectory_duplication(FFPhaseEnums::NONE);
////            ffpilotPhase->set_trajectory_generation(FFPhaseEnums::LAZY);
////        }
////        else
////        {
////            ffpilotPhase->set_trajectory_duplication(trajDuplication);
////            ffpilotPhase->set_trajectory_generation(trajGeneration);
////        }

////        // (re)initialize the relevant output options
////        addOptions(ffpilotPhase, *stage);
////    }
////}

////void FFPilotSupervisor::addOptions(FFPilotPhaseMsg* phase, const FFPilotStageMsg& stage)
////{
////    // (re)initialize the general simulation options
////    input->reinitOptions(phase->phase_id());
////    phase->mutable_options()->CopyFrom(input->getOptionsMsg());

////    // (re)initialize the relevant output options
////    input->reinitOutputOptions(phaseInfo(true, phase, &stage), stage.is_pilot_stage());
////    phase->mutable_output_options()->CopyFrom(input->getOutputOptions());
////}

void FFPilotSupervisor::startStage(int stage)
{
    if (data.getStageData(stage).pilot)
    {
        startPilotStage(stage);
    }
    else
    {
        startProductionStage(stage);
    }

    // Print a log message.
    if (ffpilotPrintStageMessages) Print::printf(Print::INFO, "FFPilot stage %d started (stage_type=%s,direction=%s).", stage, data.getStageData(stage).pilot?"Pilot":"Production", data.getStageData(stage).forward?"forward":"backward");


////    // set the current tiling to wrap the current stage's tiling msg
////    setCurrentTiling(mutableCurrentStage()->mutable_tiling());

////    // if the current stage needs pilot stage output, copy over the output from the previous stage
////    if (currentStage().needs_pilot_stage())
////    {
////        mutableCurrentStage()->mutable_pilot_stage_output()->CopyFrom(currentStageOutput().wrappedMsg());
////    }

////    // set the ffpilotPhaseLimits for this stage from its pilot stage output, if it has it
////    if (currentStage().has_pilot_stage_output())
////    {
////        addFFPilotPhaseLimitsFromStageOutput(mutableCurrentStage(), currentStageOutput());
////    }

////    // add a new stage output
////    addFFPilotStageOutput();

////    // set the first phase of the new stage as the currentFFPilotPhase
////    currentFFPilotPhaseIter = mutableCurrentStage()->mutable_ffpilot_phases()->begin();

////    // set the stage output flag
////    simulationStageOutputSent = false;

////    // start the new phase
////    startSimulationPhase();
}

void FFPilotSupervisor::startPilotStage(int stage)
{
    // Initialize each phase in the stage.
    for (size_t i=0; i<data.getStageData(stage).getNumberPhases(); i++)
    {
        data.getStageData(stage).getPhaseData(i).initialize(i==0, data.getStageData(stage).forward, static_cast<size_t>(input->ffpilot_options().pilot_stage_crossings()));
    }
}

void FFPilotSupervisor::startProductionStage(int stage)
{
    // If we skipped the pilot stage, initialize the counts from the input.
    if (input->ffpilot_options().pilot_skip())
    {
        // Make sure we have an array of trajectory counts in the input.
        if (!input->ffpilot_options().has_prod_trajectory_counts()) THROW_EXCEPTION(lm::RuntimeException, "must specifiy the trajectory counts in the input file when skipping the pilot stage");

        // Get the trajectory counts.
        ndarray<uint64_t> counts = robertslab::pbuf::NDArraySerializer::deserialize<uint64_t>(input->ffpilot_options().prod_trajectory_counts());

        // Make sure the array is the right size.
        if (counts.shape.len != 2) THROW_EXCEPTION(lm::RuntimeException, "invalid number of dimensions for trajectory counts array: %d (%d)", counts.shape.len, 2);
        if (counts.shape[0] != 2 || counts.shape[1] != tilingEdges->shape[0]) THROW_EXCEPTION(lm::RuntimeException, "invalid shape for trajectory counts: %d,%d (%d,%d)",counts[0],counts[1],2,tilingEdges->shape[0]);

        // Initialize each phase in the production stage.
        for (size_t i=0; i<data.getStageData(stage).getNumberPhases(); i++)
        {
            data.getStageData(stage).getPhaseData(i).initialize(i==0, data.getStageData(stage).forward, counts[utuple(data.getStageData(stage).forward?0:1,static_cast<uint>(i))]);
        }
    }
    else
    {
        // Otherwise, production phases are initialized after the corresponding pilot stage has finished.
    }
}


////void FFPilotSupervisor::addFFPilotStageOutput()
////{
////    // add a new phase output
////    lm::ffpilot::io::FFPilotStageOutput* newStageOutputMsg = ffpilotStageOutputsWrap.Add();

////    // set the new stage output to be the current phase output
////    currentFFPilotStageOutputWrap.setWrappedMsg(newStageOutputMsg);

////    // add a new FFPilotPhaseOutputList to go with this stage
////    currentFFPilotPhaseOutputsWrap.setWrappedField(ffpilotPhaseOutputListsWrap.Add()->mutable_ffpilot_phase_outputs());
////}

////template <typename Value>
////void* FFPilotSupervisor::addFFPilotPhaseLimit(FFPilotPhaseMsg* ffpilotPhase, FFPhaseLimEnums::StopCondition stopCondition, Value value)
////{
////    FFPilotPhaseLimitMsg* ffpilotPhaseLimit = ffpilotPhase->mutable_ffpilot_phase_limit();
    
////    ffpilotPhaseLimit->set_stop_condition(stopCondition);

////    switch (stopCondition)
////    {
////    case FFPhaseLimEnums::FORWARD_FLUXES: ffpilotPhaseLimit->set_uvalue(value); break;
////    case FFPhaseLimEnums::TRAJECTORY_COUNT: ffpilotPhaseLimit->set_uvalue(value); break;
////    case FFPhaseLimEnums::TIME: ffpilotPhaseLimit->set_dvalue(value); break;
////    }

////    buildFFPilotPhaseLimitTrajectoriesToRun(ffpilotPhaseLimit, *ffpilotPhase, slots.getSimultaneousWorkUnits());

////    return ffpilotPhaseLimit;
////}

////void FFPilotSupervisor::buildFFPilotPhaseLimitTrajectoriesToRun(FFPilotPhaseLimitMsg* ffpilotPhaseLimit, const FFPilotPhaseMsg& ffpilotPhase, uint simultaneousWorkUnits)
////{
////    // builds the events_per_trajectory field for this limit if it hasn't already been set
////    buildFFPilotPhaseLimitEventsPerTrajectory(ffpilotPhaseLimit, ffpilotPhase, simultaneousWorkUnits);

////    uint64_t simulataneousActiveTrajectories = simultaneousWorkUnits*ffpilotPhase.options().parts_per_work_unit();

////    if (ffpilotPhase.trajectory_generation()==FFPhaseEnums::EAGER)
////    {
////        // EAGER is only implemented for certain ffpilotPhaseLimit.stop_condition() values
////        if (ffpilotPhaseLimit->stop_condition()==FFPhaseLimEnums::TRAJECTORY_COUNT or (ffpilotPhaseLimit->stop_condition()==FFPhaseLimEnums::FORWARD_FLUXES and ffpilotPhase.phase_id()==0))
////        {
////            // given that our trajectory limits are set up to observe x events per trajectory, run ceil(y/x) trajectories to ensure that we observe at least y events total
////            ffpilotPhaseLimit->set_trajectories_per_phase(ceilDiv(ffpilotPhaseLimit->uvalue(), ffpilotPhaseLimit->events_per_trajectory()));
////        }
////        else THROW_EXCEPTION(UnimplementedException, "In Forward Flux phase %d, ffpilotPhase.trajectory_generation()==EAGER is only implemented for certain ffpilotPhaseLimit.stop_condition() values (ie those that let us calculate the necessary trajectory count up front). Attempting to use unimplemented ffpilotPhaseLimit.stop_condition(): %s", ffpilotPhase.phase_id(), FFPhaseLimEnums::StopCondition_Name(ffpilotPhaseLimit->stop_condition()).c_str());
////    }
////    else if (ffpilotPhase.trajectory_generation()==FFPhaseEnums::LAZY)
////    {
////        return ffpilotPhaseLimit->set_trajectories_per_phase(simulataneousActiveTrajectories);
////    }
////    else THROW_EXCEPTION(UnimplementedException, "unimplemented");
////}

////void FFPilotSupervisor::buildFFPilotPhaseLimitEventsPerTrajectory(FFPilotPhaseLimitMsg* ffpilotPhaseLimit, const FFPilotPhaseMsg& ffpilotPhase, uint simultaneousWorkUnits)
////{
////    uint64_t simulataneousActiveTrajectories = simultaneousWorkUnits*ffpilotPhase.options().parts_per_work_unit();

////    // events_per_trajectory can be set before running this function
////    if (not ffpilotPhaseLimit->has_events_per_trajectory())
////    {
////        if (ffpilotPhase.phase_id()==0)
////        {
////            if (ffpilotPhaseLimit->stop_condition()==FFPhaseLimEnums::FORWARD_FLUXES)
////            {
////                if (ffpilotPhase.trajectory_generation()==FFPhaseEnums::EAGER)
////                {
////                    ffpilotPhaseLimit->set_events_per_trajectory(ceilDiv(ffpilotPhaseLimit->uvalue(), simulataneousActiveTrajectories));
//////                    ffpilotPhaseLimit->set_events_per_trajectory(-1);
////                }
////                else if (ffpilotPhase.trajectory_generation()==FFPhaseEnums::LAZY)
////                {
////                    ffpilotPhaseLimit->set_events_per_trajectory(ffpilotPhaseLimit->uvalue());
////                }
////            }
////            else THROW_EXCEPTION(UnimplementedException, "ffpilotPhaseLimit->stop_condition()==TIME and ffpilotPhaseLimit->stop_condition()==TRAJECTORY_COUNT currently unimplemented for ffpilot phase 0");
////        }
////        else
////        {
////            ffpilotPhaseLimit->set_events_per_trajectory(1);
////        }
////    }
////}

////template <typename Value>
////void FFPilotSupervisor::addFFPilotPhaseLimitsForPilotStage(FFPilotStageMsg* stage, FFPhaseLimEnums::StopCondition stopCondition, Value phaseZeroValue, Value value)
////{
////    // get iterator over the phases
////    FFPilotPhasesWrap::iterator phase_it=stage->mutable_ffpilot_phases()->begin();

////    // add limits to phase 0
////    addFFPilotPhaseLimit(&*phase_it, stopCondition, phaseZeroValue);
////    phase_it++;

////    // add limits to the other phases
////    for (; phase_it!=stage->mutable_ffpilot_phases()->end(); phase_it++)
////    {
////        addFFPilotPhaseLimit(&*phase_it, stopCondition, value);
////    }
////}

////void FFPilotSupervisor::addFFPilotPhaseLimitsFromInput(FFPilotStageMsg* productionStage)
////{
////    // TODO: implement manually specified ffpilotPhaseLimits
////    //productionStage->mutable_ffpilot_phase_limits()->CopyFrom(input->getFFPilotPhaseLimits(productionStage->tiling().id(), productionStage->basin_id()));

////    // temporary placeholder
////    addFFPilotPhaseLimitsForPilotStage(productionStage, FFPhaseLimEnums::FORWARD_FLUXES, input->getFFPilotOptionsMsg().pilot_stage_count(), input->getFFPilotOptionsMsg().pilot_stage_count());
////}

////void FFPilotSupervisor::addFFPilotPhaseLimitsFromStageOutput(FFPilotStageMsg* productionStage, const lm::protowrap::FFPilotStageOutputWrap& stageOutput)
////{
////    // calculate the optimum trajectory counts for the production stage from the pilot stage output
////    vector<uint64_t> trajectoryCounts(optimizeTrajectoryCounts(input->errorGoal(), input->errorGoalConfidence(), stageOutput, input->productionStageCountMinimum(), input->phaseZeroSamplingMultiplier(), input->minimizeCost()));

////    // print some info about the pilot stage (costs, probabilities, etc) to the log (ie stdout)
////    Print::printf(Print::INFO, stageLogPilot(stageOutput, input->errorGoal(), input->errorGoalConfidence(), trajectoryCounts).c_str());

////    // zip over the trajectory counts and the phase specification msgs
////    vector<uint64_t>::const_iterator trajcount_it=trajectoryCounts.begin();
////    FFPilotPhasesWrap::iterator phase_it=productionStage->mutable_ffpilot_phases()->begin();

////    // special treatment for phase zero
////    addFFPilotPhaseLimit(&*phase_it, FFPhaseLimEnums::FORWARD_FLUXES, *trajcount_it);
////    trajcount_it++, phase_it++;

////    // all phases n>0
////    for (; trajcount_it!=trajectoryCounts.end() and phase_it!=productionStage->mutable_ffpilot_phases()->end(); trajcount_it++, phase_it++)
////    {
////        addFFPilotPhaseLimit(&*phase_it, FFPhaseLimEnums::TRAJECTORY_COUNT, *trajcount_it);
////    }
////}

////void FFPilotSupervisor::repeatFFPilotPhaseLimits(FFPilotStageMsg* stage, const FFPilotPhaseLimitMsg& limitToRepeat)
////{
////    // add copies of limitToRepeat for every ffpilotPhase that's missing a corresponding ffpilotPhaseLimit
////    for (FFPilotPhasesWrap::iterator phase_it = stage->mutable_ffpilot_phases()->begin(); phase_it != stage->mutable_ffpilot_phases()->end(); phase_it++)
////    {
////        if (not phase_it->has_ffpilot_phase_limit())
////        {
////            phase_it->mutable_ffpilot_phase_limit()->CopyFrom(limitToRepeat);
////        }
////    }
////}

void FFPilotSupervisor::startPhase(int stageIndex, int phaseIndex)
{
    StageData& stage = data.getStageData(stageIndex);
    PhaseData& phase = stage.getPhaseData(phaseIndex);

    // Create the trajectories for the phase.
    buildPhase(stageIndex, phaseIndex);

//    // add a new phase output
//    addFFPilotPhaseOutput();

//    // set the first trajectory id of this phase in the phase output. If there is an old trajectoryList, get the trajectory count from that. Otherwise, we're at the very start of the simulation so count is 0
//    currentFFPilotPhaseOutputWrapPtr->set_first_trajectory_id(trajectoryList != NULL ? trajectoryList->count() : 0);

//    // reset the phase output and termination flag
//    simulationPhaseOutputSent = false;
//    ffpilotPhaseTerminated = false;

    // Print a log message.
    if (phaseIndex == 0)
    {
        if (ffpilotPrintStageMessages) Print::printf(Print::INFO, "FFPilot stage %d phase %d started (goal_edge=%0.3e, phase_limit=%s>=%d).", stageIndex, phaseIndex, phase.goalEdge, stage.pilot?"FORWARD_FLUXES":"TRAJECTORIES", phase.phaseLimit);
    }
    else
    {
        if (ffpilotPrintStageMessages) Print::printf(Print::INFO, "FFPilot stage %d phase %d started (starting_edge=%0.3e, goal_edge=%0.3e, phase_limit=%s>=%d).", stageIndex, phaseIndex, stage.getPhaseData(phaseIndex-1).goalEdge, phase.goalEdge, stage.pilot?"FORWARD_FLUXES":"TRAJECTORIES", phase.phaseLimit);
    }
}

void FFPilotSupervisor::buildPhase(int stage, int phase)
{
    // Create a new trajectory list.
    trajectoryList = new FFPilotTrajectoryList();

    if (phase == 0)
        buildPhase0(stage);
    else
        buildPhaseN(stage, phase);
//    trajectoryList = new lm::replicates::ReplicateTrajectoryList(*input, 1, 1);

    // if there is an old trajectoryList, get the trajectory count from that. Otherwise, we're at the very start of the simulation so count is 0
//    uint64_t currentTrajectoryCount = trajectoryList != NULL ? trajectoryList->count() : 0;

    // set the trajectory limits/tracking for this phase
//    input->reinitTrajectoryLimits(currentPhase(), currentTiling());

//    if (currentPhase().start_points_size() > 0)
//    {
//        setTrajectoryList(new FFPilotTrajectoryList(currentTrajectoryCount, currentFFPilotPhaseID(), currentPhase(), slots.getSimultaneousWorkUnits(), *input));  //, &ffpilotPhaseOutputMsgCustom));
//    }
//    else if(currentFFPilotPhaseID()==0)
//    {
//        setTrajectoryList(new FFPilotTrajectoryList(currentTrajectoryCount, currentFFPilotPhaseID(), currentPhase(), slots.getSimultaneousWorkUnits(), *input, currentTiling().currentBasin()));
//    }
//    else
//    {
//        setTrajectoryList(new FFPilotTrajectoryList(currentTrajectoryCount, currentFFPilotPhaseID(), currentPhase(), slots.getSimultaneousWorkUnits(), *input, previousPhaseOutput()));
//    }
}

void FFPilotSupervisor::buildPhase0(int stageIndex)
{
    if (input->ffpilot_options().fallback_method() == lm::input::ffpilot::FFPilotOptions::FALLBACK_TILING_EDGE)
        buildPhase0FallbackTiling(stageIndex);
    else if (input->ffpilot_options().fallback_method() == lm::input::ffpilot::FFPilotOptions::FALLBACK_BASIN)
        buildPhase0FallbackBasin(stageIndex);
    else
        THROW_EXCEPTION(lm::RuntimeException, "unkown ffpilot fallback method: %d", input->ffpilot_options().fallback_method());
}

void FFPilotSupervisor::buildPhase0FallbackTiling(int stageIndex)
{
    StageData& stage = data.getStageData(stageIndex);
    PhaseData& phase = stage.getPhaseData(0);

    // Get the initial species counts from the basin.
    ndarray<int32_t> basins = robertslab::pbuf::NDArraySerializer::deserialize<int32_t>(input->ffpilot_options().phase_zero_basins());
    ndarray<int32_t> initialSpeciesCounts(basins.shape[1]);
    uint row = stage.forward?0:1;
    for (uint i=0; i<initialSpeciesCounts.shape[0]; i++)
        initialSpeciesCounts[i] = basins[utuple(row,i)];

    // Create a trajectory for each simultaneous work unit we can process.
    for (uint i=0; i<slotList->getSimultaneousWorkUnits(); i++)
    {
        // Create a trajectory.
        lm::cme::CMETrajectory* trajectory;
        trajectory = new lm::cme::CMETrajectory(nextTrajectoryId++);

        // Initialize the CME state.
        trajectory->initializeSpeciesCountsState(initialSpeciesCounts, 0.0);

        // Initialize any order parameters state.
        if (input->has_order_parameters()) trajectory->initializeOrderParametersState(input->order_parameters());

        // Add the trajectory to the list.
        trajectoryList->addTrajectory(trajectory->getID(), trajectory);

        // Initialize the phase zero data for this trajectory.
        phase.phase0NumberCrossings[trajectory->getID()]  = 0;
        phase.phase0LastCrossingDirection[trajectory->getID()] = 1;
        phase.phase0LastCrossingTime[trajectory->getID()]  = 0.0;
        phase.phase0LastCrossingStep[trajectory->getID()]  = 0;
        phase.phase0LastCrossingState[trajectory->getID()] = initialSpeciesCounts;
    }
    Print::printf(Print::DEBUG, "Built trajectories for phase 0: %d", trajectoryList->size());

    // Set the reflecting barrier at the transition point.
    lm::types::TrajectoryBarrier* b1 = trajectoryBarriers.add_barrier();
    b1->set_type(lm::types::TrajectoryBarrier::ORDER_PARAMETER);
    b1->set_type_arg(static_cast<uint32_t>(orderParameterIndex));
    b1->set_behavior(lm::types::TrajectoryBarrier::REFLECTING);
    if (input->ffpilot_options().has_phase_zero_transition_barrier())
        b1->set_behavior_value_double(input->ffpilot_options().phase_zero_transition_barrier());
    else
        b1->set_behavior_value_double(tilingEdges->get(tilingEdges->shape[0]/2));
    Print::printf(Print::DEBUG, "Set the reflecting barrier for phase 0: %0.3e", b1->behavior_value_double());

    // Set the tracking barriers at the phase zero edge.
    lm::types::TrajectoryBarrier* b2 = trajectoryBarriers.add_barrier();
    b2->set_type(lm::types::TrajectoryBarrier::ORDER_PARAMETER);
    b2->set_type_arg(static_cast<uint32_t>(orderParameterIndex));
    lm::types::TrajectoryBarrier* b3 = trajectoryBarriers.add_barrier();
    b3->set_type(lm::types::TrajectoryBarrier::ORDER_PARAMETER);
    b3->set_type_arg(static_cast<uint32_t>(orderParameterIndex));
    if (stage.forward)
    {
        b2->set_behavior(lm::types::TrajectoryBarrier::TRACKING_INCREASING_INCLUSIVE);
        b2->set_behavior_value_double(tilingEdges->get(0));
        b3->set_behavior(lm::types::TrajectoryBarrier::TRACKING_DECREASING_EXCLUSIVE);
        b3->set_behavior_value_double(tilingEdges->get(0));
    }
    else
    {
        b2->set_behavior(lm::types::TrajectoryBarrier::TRACKING_DECREASING_INCLUSIVE);
        b2->set_behavior_value_double(tilingEdges->get(tilingEdges->shape[0]-1));
        b3->set_behavior(lm::types::TrajectoryBarrier::TRACKING_INCREASING_EXCLUSIVE);
        b3->set_behavior_value_double(tilingEdges->get(tilingEdges->shape[0]-1));
    }
    phase.goalEdge = b2->behavior_value_double();
    Print::printf(Print::DEBUG, "Set the tracking barrier for phase 0: %0.3e (type %d)", b2->behavior_value_double(), static_cast<int>(b2->behavior()));

    // Set a limits for this trajectory list as the number of tracking barrier crossings.
    lm::types::TrajectoryLimit* l = trajectoryLimits.add_limits();
    l->set_id(trajectoryLimits.limits_size()-1);
    l->set_type(lm::types::TrajectoryLimit::BARRIER_CROSSING);
    l->set_type_arg(2); // Use the "backward" crossing as the limit so that we have an equal number of forward and backward crossings.
    l->set_stopping_condition(lm::types::TrajectoryLimit::MAX_INCLUSIVE);
    l->set_stopping_value_int(static_cast<int32_t>((phase.phaseLimit/trajectoryList->size())) + 1 + input->ffpilot_options().phase_zero_equilibration_crossings());
    Print::printf(Print::DEBUG, "Set the trajectory limit for phase 0: %d x %d", trajectoryList->size(), l->stopping_value_int());
}

void FFPilotSupervisor::buildPhase0FallbackBasin(int stageIndex)
{
    StageData& stage = data.getStageData(stageIndex);
    PhaseData& phase = stage.getPhaseData(0);

    // Get the initial species counts from the basin.
    ndarray<int32_t> basins = robertslab::pbuf::NDArraySerializer::deserialize<int32_t>(input->ffpilot_options().phase_zero_basins());
    ndarray<int32_t> initialSpeciesCounts(basins.shape[1]);
    uint row = stage.forward?0:1;
    for (uint i=0; i<initialSpeciesCounts.shape[0]; i++)
        initialSpeciesCounts[i] = basins[utuple(row,i)];

    // Create a trajectory for each crossing we need.
    for (size_t i=0; i<phase.phaseLimit; i++)
    {
        // Create a trajectory.
        lm::cme::CMETrajectory* trajectory;
        trajectory = new lm::cme::CMETrajectory(nextTrajectoryId++);

        // Initialize the CME state.
        trajectory->initializeSpeciesCountsState(initialSpeciesCounts, 0.0);

        // Initialize any order parameters state.
        if (input->has_order_parameters()) trajectory->initializeOrderParametersState(input->order_parameters());

        // Add the trajectory to the trajectory list.
        trajectoryList->addTrajectory(trajectory->getID(), trajectory);

        // Add the trajectory to the ffpilot data.
        TrajectoryData t(trajectory->getID(), orderParameterFunction->calculate(initialSpeciesCounts), initialSpeciesCounts);
        phase.trajectoryData[t.id] = t;
    }
    Print::printf(Print::DEBUG, "Built trajectories for phase 0: %d", trajectoryList->size());

    // Set the tracking barriers at the phase zero edge.
    lm::types::TrajectoryBarrier* b0 = trajectoryBarriers.add_barrier();
    b0->set_type(lm::types::TrajectoryBarrier::ORDER_PARAMETER);
    b0->set_type_arg(static_cast<uint32_t>(orderParameterIndex));
    if (stage.forward)
    {
        b0->set_behavior(lm::types::TrajectoryBarrier::TRACKING_INCREASING_INCLUSIVE);
        b0->set_behavior_value_double(tilingEdges->get(0));
    }
    else
    {
        b0->set_behavior(lm::types::TrajectoryBarrier::TRACKING_DECREASING_INCLUSIVE);
        b0->set_behavior_value_double(tilingEdges->get(tilingEdges->shape[0]-1));
    }
    phase.goalEdge = b0->behavior_value_double();
    Print::printf(Print::DEBUG, "Set the tracking barrier for phase 0: %0.3e (type %d)", b0->behavior_value_double(), static_cast<int>(b0->behavior()));

    // Set a limits for this trajectory list as the number of tracking barrier crossings.
    lm::types::TrajectoryLimit* l = trajectoryLimits.add_limits();
    l->set_id(trajectoryLimits.limits_size()-1);
    l->set_type(lm::types::TrajectoryLimit::BARRIER_CROSSING);
    l->set_type_arg(0);
    l->set_stopping_condition(lm::types::TrajectoryLimit::MAX_INCLUSIVE);
    l->set_stopping_value_int(1);
    Print::printf(Print::DEBUG, "Set the trajectory limit for phase 0: %d", trajectoryList->size(), l->stopping_value_int());
}


void FFPilotSupervisor::buildPhaseN(int stageIndex, int phaseIndex)
{
    StageData& stage = data.getStageData(stageIndex);
    PhaseData& phase = stage.getPhaseData(phaseIndex);

    // Create a trajectory for each crossing we need.
    for (size_t i=0; i<phase.phaseLimit; i++)
    {
        createPhaseNTrajectory(stageIndex, phaseIndex);
    }
    Print::printf(Print::DEBUG, "Built trajectories for phase %d: %d", phaseIndex, trajectoryList->size());

    // Set the forward limit barrier at this interface edge.
    lm::types::TrajectoryLimit* l0 = trajectoryLimits.add_limits();
    l0->set_id(trajectoryLimits.limits_size()-1);
    l0->set_type(lm::types::TrajectoryLimit::ORDER_PARAMETER);
    l0->set_type_arg(static_cast<uint32_t>(orderParameterIndex));
    if (stage.forward)
    {
        l0->set_stopping_condition(lm::types::TrajectoryLimit::INCREASING_INCLUSIVE);
        l0->set_stopping_value_double(tilingEdges->get(static_cast<uint>(phaseIndex)));
    }
    else
    {
        l0->set_stopping_condition(lm::types::TrajectoryLimit::DECREASING_INCLUSIVE);
        l0->set_stopping_value_double(tilingEdges->get(tilingEdges->shape[0]-1-static_cast<uint>(phaseIndex)));
    }
    phase.goalEdge = l0->stopping_value_double();
    Print::printf(Print::DEBUG, "Set the forward limit for phase %d: %0.3e (type %d)", phaseIndex, l0->stopping_value_double(), static_cast<int>(l0->stopping_condition()));

    // Create the backward limit.
    lm::types::TrajectoryLimit* l1 = trajectoryLimits.add_limits();

    if (input->ffpilot_options().fallback_method() == lm::input::ffpilot::FFPilotOptions::FALLBACK_TILING_EDGE)
    {
        // Set the backward limit barrier at the phase zero edge.
        l1->set_id(trajectoryLimits.limits_size()-1);
        l1->set_type(lm::types::TrajectoryLimit::ORDER_PARAMETER);
        l1->set_type_arg(static_cast<uint32_t>(orderParameterIndex));
        if (stage.forward)
        {
            l1->set_stopping_condition(lm::types::TrajectoryLimit::DECREASING_EXCLUSIVE);
            l1->set_stopping_value_double(tilingEdges->get(0));
        }
        else
        {
            l1->set_stopping_condition(lm::types::TrajectoryLimit::INCREASING_EXCLUSIVE);
            l1->set_stopping_value_double(tilingEdges->get(tilingEdges->shape[0]-1));
        }
    }
    else if (input->ffpilot_options().fallback_method() == lm::input::ffpilot::FFPilotOptions::FALLBACK_BASIN)
    {
        // Set the backward limit barrier at the basin.
        ndarray<int32_t> basins = robertslab::pbuf::NDArraySerializer::deserialize<int32_t>(input->ffpilot_options().phase_zero_basins());
        l1->set_id(trajectoryLimits.limits_size()-1);
        l1->set_type(lm::types::TrajectoryLimit::ORDER_PARAMETER);
        l1->set_type_arg(static_cast<uint32_t>(orderParameterIndex));
        if (stage.forward)
        {
            ndarray<int32_t> initialSpeciesCounts(basins.shape[1]);
            for (uint i=0; i<initialSpeciesCounts.shape[0]; i++)
                initialSpeciesCounts[i] = basins[utuple(0,i)];
            double basinLimitValue = orderParameterFunction->calculate(0.0, initialSpeciesCounts.values, initialSpeciesCounts.shape[0]);
            l1->set_stopping_condition(lm::types::TrajectoryLimit::DECREASING_INCLUSIVE);
            l1->set_stopping_value_double(basinLimitValue);
        }
        else
        {
            ndarray<int32_t> initialSpeciesCounts(basins.shape[1]);
            for (uint i=0; i<initialSpeciesCounts.shape[0]; i++)
                initialSpeciesCounts[i] = basins[utuple(1,i)];
            double basinLimitValue = orderParameterFunction->calculate(0.0, initialSpeciesCounts.values, initialSpeciesCounts.shape[0]);
            l1->set_stopping_condition(lm::types::TrajectoryLimit::INCREASING_INCLUSIVE);
            l1->set_stopping_value_double(basinLimitValue);
        }
    }
    else
    {
        THROW_EXCEPTION(lm::RuntimeException, "unkown ffpilot fallback method: %d", input->ffpilot_options().fallback_method());
    }

    Print::printf(Print::DEBUG, "Set the backward limit for phase %d: %0.3e (type %d)", phaseIndex, l1->stopping_value_double(), static_cast<int>(l1->stopping_condition()));
}

void FFPilotSupervisor::createPhaseNTrajectory(int stageIndex, int phaseIndex)
{
    StageData& stage = data.getStageData(stageIndex);
    PhaseData& phase = stage.getPhaseData(phaseIndex);
    PhaseData& previousPhase = stage.getPhaseData(phaseIndex-1);

    // Create a trajectory.
    lm::cme::CMETrajectory* trajectory;
    trajectory = new lm::cme::CMETrajectory(nextTrajectoryId++);

    // Sample a random crossing from the previous phase.
    size_t previousCrossingIndex = random->getRandomInRange(0, static_cast<uint32_t>(previousPhase.successfulCrossings.size()));
    uint64_t previousTrajectoryId = previousPhase.successfulCrossings[previousCrossingIndex];

    // Get the ending location of the crossing.
    TrajectoryData previousCrossing = previousPhase.trajectoryData[previousTrajectoryId];
    ndarray<int32_t> previousState = previousCrossing.endingState;

    // Initialize the CME state.
    trajectory->initializeSpeciesCountsState(previousState, 0.0);

    // Initialize any order parameters state.
    if (input->has_order_parameters()) trajectory->initializeOrderParametersState(input->order_parameters());

    // Add the trajectory to the trajectory list.
    trajectoryList->addTrajectory(trajectory->getID(), trajectory);

    // Add the trajectory to the ffpilot data.
    TrajectoryData t(previousTrajectoryId, trajectory->getID(), orderParameterFunction->calculate(previousState), previousState);
    phase.trajectoryData[t.id] = t;
}

bool FFPilotSupervisor::assignWork(int stage, int phase)
{
    // While we have free slots and waiting trajectories, assign them.
    while (slotList->hasFreeSlots() && trajectoryList->anyWaiting())
    {
        // Create the run work unit message.
        lm::message::Message msg;
        lm::message::RunWorkUnit* rwuMsg = msg.mutable_run_work_unit();

        // Build the run work units message.
        buildRunWorkUnitHeader(rwuMsg);

        // Build the work unit parts.
        buildRunWorkUnitParts(rwuMsg);

        // Run the work unit.
        slotList->runWorkUnit(communicator, &msg);
    }

    return !trajectoryList->allFinished();
}

void FFPilotSupervisor::buildRunWorkUnitHeader(lm::message::RunWorkUnit* msg)
{
    // Call the base class method.
    StagedSimulationSupervisor::buildRunWorkUnitHeader(msg);

    // Set the solver to be master equation.
    msg->set_solver_type(lm::types::SolverType::ME);

    // Set the limits.
    msg->mutable_trajectory_limits()->CopyFrom(trajectoryLimits);

    // Set the barriers.
    msg->mutable_trajectory_barriers()->CopyFrom(trajectoryBarriers);

    // Set the output options.
    msg->mutable_output_options()->CopyFrom(outputOptions);

    // Override the output address so that the ffpilot supervisor receives the output messages.
    msg->mutable_output_address()->CopyFrom(communicator->getSourceAddress());

    // Set the maximum number of steps for the work unit.
    msg->set_max_steps(input->simulation_options().steps_per_work_unit_part());

    // Add the model.
    if (input->has_reaction_model()) msg->mutable_reaction_model()->CopyFrom(input->reaction_model());
    if (input->has_diffusion_model()) msg->mutable_diffusion_model()->CopyFrom(input->diffusion_model());
    if (input->has_order_parameters()) msg->mutable_order_parameters()->CopyFrom(input->order_parameters());
    //if (input->has_tilings()) msg->mutable_()->CopyFrom(input->order_parameters());
}

void FFPilotSupervisor::buildRunWorkUnitParts(lm::message::RunWorkUnit* msg)
{
    // Figure out how many parts if we equally divide the remaining work units by the remaining slots.
    size_t waitingTrajectories = trajectoryList->numberWaiting();
    size_t freeSlots = slotList->numberFreeSlots();
    size_t trajectoriesPerSlot = waitingTrajectories/freeSlots;

    // Make sure we are less than the maximum number of parts per work unit from the configuration.
//    trajectoriesPerSlot = std::min(input->ffpilot_options().max_trajectories_per_work_unit(), trajectoriesPerSlot);
	trajectoriesPerSlot = std::min(static_cast<size_t>(input->ffpilot_options().max_trajectories_per_work_unit()), trajectoriesPerSlot);

    // Make sure we are greater than the minimum number of parts needed for the next slot.
    trajectoriesPerSlot = std::max(static_cast<size_t>(slotList->getFreeSlot().getSimultaneousWorkUnits()), trajectoriesPerSlot);

    // Build the work units.
    trajectoryList->buildWorkUnitParts(msg->work_unit_id(), msg, trajectoriesPerSlot);

    //Print::printf(Print::DEBUG, "FFPilot built work unit with %d parts", msg->part_size());
}

bool FFPilotSupervisor::processEvent(lm::message::Message& msg)
{
    StageData& stage = data.getStageData(currentStage);
    PhaseData& phase = stage.getPhaseData(currentPhase);

    // See if this is a message we want to process.
    if (msg.has_process_work_unit_output())
    {
        int64_t workUnitId = msg.process_work_unit_output().work_unit_id();
        for (int i=0; i<msg.process_work_unit_output().part_output().size(); i++)
        {
            if (input->ffpilot_options().fallback_method() == lm::input::ffpilot::FFPilotOptions::FALLBACK_TILING_EDGE)
            {
                // See if we need to process the phase zero barrier crossing.
                if (currentPhase == 0 && msg.process_work_unit_output().part_output(i).barrier_crossing_times().size() == 2)
                {
                    processPhaseZeroBarrierCrossingTimesFallbackTiling(stage, phase, workUnitId, msg.process_work_unit_output().part_output(i).barrier_crossing_times(0), msg.process_work_unit_output().part_output(i).barrier_crossing_times(1));
                }
            }
            else if (input->ffpilot_options().fallback_method() == lm::input::ffpilot::FFPilotOptions::FALLBACK_BASIN)
            {
                // See if we need to process the phase zero barrier crossing.
                if (currentPhase == 0 && msg.process_work_unit_output().part_output(i).barrier_crossing_times().size() == 1)
                {
                    processPhaseZeroBarrierCrossingTimesFallbackBasin(stage, phase, workUnitId, msg.process_work_unit_output().part_output(i).barrier_crossing_times(0));
                }
            }
            else
            {
                THROW_EXCEPTION(lm::RuntimeException, "unkown ffpilot fallback method: %d", input->ffpilot_options().fallback_method());
            }
        }
        return true;
    }

    // Othwerwise, see if our parent can process the message.
    return StagedSimulationSupervisor::processEvent(msg);
}

void FFPilotSupervisor::processPhaseZeroBarrierCrossingTimesFallbackTiling(StageData& stage, PhaseData& phase, int64_t workUnitId, const lm::io::BarrierCrossingTimes& data0, const lm::io::BarrierCrossingTimes& data1)
{
    if (data0.trajectory_id() != data1.trajectory_id()) THROW_EXCEPTION(lm::RuntimeException, "inconsistent trajectory id between barrier corssing times: %d, %d", static_cast<int>(data0.trajectory_id()), static_cast<int>(data1.trajectory_id()));
    if (data0.barrier_index() != 0 || data1.barrier_index() != 1) THROW_EXCEPTION(lm::RuntimeException, "unexpected barrier crossing indices: %d, %d", static_cast<int>(data0.barrier_index()), static_cast<int>(data1.barrier_index()));

    // Get the trajectory id.
    uint64_t trajectoryId = data0.trajectory_id();

    // Get the counts, times, and steps.
    ndarray<int32_t> counts[2];
    counts[0] = robertslab::pbuf::NDArraySerializer::deserialize<int32_t>(data0.counts());
    counts[1] = robertslab::pbuf::NDArraySerializer::deserialize<int32_t>(data1.counts());
    ndarray<double> times[2];
    times[0] = robertslab::pbuf::NDArraySerializer::deserialize<double>(data0.times());
    times[1] = robertslab::pbuf::NDArraySerializer::deserialize<double>(data1.times());
    ndarray<uint64_t> steps[2];
    steps[0] = robertslab::pbuf::NDArraySerializer::deserialize<uint64_t>(data0.total_steps());
    steps[1] = robertslab::pbuf::NDArraySerializer::deserialize<uint64_t>(data1.total_steps());

    // Walk through the crossings and add the crossing times to the phase.
    size_t i[2] = {0,0};
    size_t previousDirection = phase.phase0LastCrossingDirection[trajectoryId];
    double previousTime = phase.phase0LastCrossingTime[trajectoryId];
    uint64_t previousStep = phase.phase0LastCrossingStep[trajectoryId];
    ndarray<int32_t> previousState = phase.phase0LastCrossingState[trajectoryId];

    while (true)
    {
        // Get the next direction.
        size_t currentDirection = previousDirection?0:1;

        // See if we are past the length of the crossing list.
        if (i[currentDirection] >= static_cast<size_t>(times[currentDirection].shape[0]))
        {
            // If the the other crossing list isn't empty, something has gone wrong.
            if (i[previousDirection] != times[previousDirection].shape[0]) THROW_EXCEPTION(lm::RuntimeException, "inconsistency in the crossing time list: %lu,%lu", i[previousDirection], times[previousDirection].shape[0]);

            // Save the last direction and time for the phase and stop the loop.
            phase.phase0LastCrossingDirection[trajectoryId] = previousDirection;
            phase.phase0LastCrossingTime[trajectoryId] = previousTime;
            phase.phase0LastCrossingState[trajectoryId] = previousState;
            phase.phase0LastCrossingStep[trajectoryId] = previousStep;
            break;
        }

        // Get the next time and state.
        double currentTime = times[currentDirection].get(static_cast<uint>(i[currentDirection]));
        uint64_t currentStep = steps[currentDirection].get(static_cast<uint>(i[currentDirection]));
        ndarray<int32_t> currentState(utuple(counts[currentDirection].shape[1]));
        for (uint j=0; j<currentState.shape[0]; j++)
            currentState[utuple(j)] = counts[currentDirection].get(utuple(static_cast<uint>(i[currentDirection]),j));

        // Make sure the current time is after the last time.
        if (currentTime < previousTime)
        {
            times[0].shape.print();
            times[1].shape.print();
            Print::printf(Print::INFO, "Invalid time times 0: %0.30e,%0.30e,%0.30e", times[0][i[0]-2],times[0][i[0]-1],times[0][i[0]]);
            Print::printf(Print::INFO, "Invalid time times 1: %0.30e,%0.30e,%0.30e", times[1][i[1]-2],times[1][i[1]-1],times[1][i[1]]);
            Print::printf(Print::INFO, "Invalid time steps 0: %ld,%ld,%ld", steps[0][i[0]-2],steps[0][i[0]-1],steps[0][i[0]]);
            Print::printf(Print::INFO, "Invalid time steps 1: %ld,%ld,%ld", steps[1][i[1]-2],steps[1][i[1]-1],steps[1][i[1]]);
            Print::printf(Print::INFO, "Invalid time indices: %d,%d,%d", currentDirection,i[0],i[1]);
            Print::printf(Print::INFO, "Invalid time last: %ld,%d,%0.8e", phase.phase0LastCrossingStep[trajectoryId],phase.phase0LastCrossingDirection[trajectoryId],phase.phase0LastCrossingTime[trajectoryId]);
            THROW_EXCEPTION(lm::RuntimeException, "invalid crossing time received: %0.30e,%0.30e %ld,%ld", previousTime, currentTime, previousStep, currentStep);
        }

        // See if this is a forward flux event.
        if (currentDirection == 0)
        {
            // See if we need to skip this event.
            if (phase.phase0NumberCrossings[trajectoryId] < input->ffpilot_options().phase_zero_equilibration_crossings())
            {
            }
            else if (phase.trajectoryData.size() < phase.phaseLimit)
            {
                // Add the event to the list.
                uint64_t trajectoryPartId = (trajectoryId+1)*100000000000000UL + static_cast<uint64_t>(phase.phase0NumberCrossings[trajectoryId]);
                TrajectoryData t(trajectoryPartId, orderParameterFunction->calculate(previousState), previousState);
                t.endingValue = orderParameterFunction->calculate(currentState);
                t.endingState = currentState;
                t.time = currentTime-previousTime;
                t.steps = currentStep-previousStep;
                phase.trajectoryData[trajectoryPartId] = t;
                phase.successfulCrossings.push_back(trajectoryPartId);
                Print::printf(Print::VERBOSE_DEBUG, "Added a phase zero crossing: %ld,%ld,%0.6e,%ld %d->%d %0.2f->%0.2f", phase.trajectoryData.size(),t.id, t.time, t.steps, t.startingState[0], t.endingState[0], t.startingValue, t.endingValue);
            }
            else
            {
                Print::printf(Print::DEBUG, "Warning: skipped an extraneous phase zero crossing");
            }

            // Increment the crossing counter.
            phase.phase0NumberCrossings[trajectoryId]++;
        }

        // Update the counter.
        i[currentDirection]++;

        // Update the previous positions.
        previousDirection = currentDirection;
        previousTime = currentTime;
        previousStep = currentStep;
        previousState = currentState;
    }
}

void FFPilotSupervisor::processPhaseZeroBarrierCrossingTimesFallbackBasin(StageData& stage, PhaseData& phase, int64_t workUnitId, const lm::io::BarrierCrossingTimes& data)
{
    // Get the trajectory id.
    uint64_t trajectoryId = data.trajectory_id();

    // Get the counts and the times.
    ndarray<int32_t> counts = robertslab::pbuf::NDArraySerializer::deserialize<int32_t>(data.counts());
    ndarray<double> times = robertslab::pbuf::NDArraySerializer::deserialize<double>(data.times());
    ndarray<uint64_t> steps = robertslab::pbuf::NDArraySerializer::deserialize<uint64_t>(data.total_steps());

    // Check the data.
    if (counts.shape[0] != 1 || times.shape[0] != 1 || steps.shape[0] != 1) THROW_EXCEPTION(lm::RuntimeException, "invalid number of rows in phase 0 crossing data: %d,%d,%d",counts.shape[0], times.shape[0], steps.shape[0]);

    // Update the trajectory data with the crossing time and steps.
    phase.trajectoryData[trajectoryId].time = times[0];
    phase.trajectoryData[trajectoryId].steps = steps[0];

    // Update the trajectory data with the crossing state.
    ndarray<int32_t> crossingState(utuple(counts.shape[1]));
    for (uint j=0; j<crossingState.shape[0]; j++)
        crossingState[j] = counts[utuple(0,j)];
    phase.trajectoryData[trajectoryId].endingValue = orderParameterFunction->calculate(crossingState);
    phase.trajectoryData[trajectoryId].endingState = crossingState;
    phase.successfulCrossings.push_back(trajectoryId);

    Print::printf(Print::VERBOSE_DEBUG, "Added a phase zero crossing: %ld,%ld,%0.6e,%ld %d->%d %0.2f->%0.2f", phase.trajectoryData.size(), phase.trajectoryData[trajectoryId].id, phase.trajectoryData[trajectoryId].time, phase.trajectoryData[trajectoryId].steps, phase.trajectoryData[trajectoryId].startingState[0], phase.trajectoryData[trajectoryId].endingState[0],phase.trajectoryData[trajectoryId].startingValue, phase.trajectoryData[trajectoryId].endingValue);
}

void FFPilotSupervisor::receivedFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg)
{
    // Update the trajectory list.
    trajectoryList->processWorkUnitFinished(msg.work_unit_id(), msg);

    // See if this is a message we want to process.
    if (currentPhase > 0)
    {
        // Go through the parts that finished.
        for (int i=0; i<msg.part_status().size(); i++)
        {
            // See if an order parameter limit was reached.
            if (msg.part_status(i).status() == lm::message::WorkUnitStatus::LIMIT_REACHED && msg.part_status(i).limit_reached().type() == lm::types::TrajectoryLimit::ORDER_PARAMETER)
            {
                processPhaseNTrajectoryFinished(currentStage, currentPhase, msg.part_status(i).final_state(), msg.part_status(i).limit_reached().id() == 0);
            }
        }
    }

    // Call the base class method.
    StagedSimulationSupervisor::receivedFinishedWorkUnit(msg);
}

void FFPilotSupervisor::processPhaseNTrajectoryFinished(int stageIndex, int phaseIndex, const lm::io::TrajectoryState& state, bool success)
{
    StageData& stage = data.getStageData(stageIndex);
    PhaseData& phase = stage.getPhaseData(phaseIndex);

    // Get the trajectory id.
    uint64_t trajectoryId = state.trajectory_id();

    TrajectoryData& t = phase.trajectoryData[trajectoryId];

    // Update the trajectory data with the crossing time and steps.
    t.time = state.cme_state().species_counts().time(0);
    t.steps = state.cme_state().total_steps();

    // Update the trajectory data with the final state.
    ndarray<int32_t> finalState(utuple(static_cast<uint>(state.cme_state().species_counts().number_species())));
    for (uint i=0; i<finalState.shape[0]; i++)
        finalState[i] = state.cme_state().species_counts().species_count(i);
    phase.trajectoryData[trajectoryId].endingValue = orderParameterFunction->calculate(finalState);
    t.endingState = finalState;

    // See if the trajectory crossed successfully.
    if (success)
    {
        // Add it to the successful list.
        phase.successfulCrossings.push_back(trajectoryId);
    }
    else
    {
        // Add it to the failed list.
        phase.failedCrossings.push_back(trajectoryId);

        // If this is the pilot stage we are counting only successful trajectories, so start another one.
        if (stage.pilot)
        {
            createPhaseNTrajectory(stageIndex, phaseIndex);
        }
    }

    Print::printf(Print::VERBOSE_DEBUG, "Finished phase %d trajectory: %ld(<-%ld) %0.6e,%ld %d->%d %0.2f->%0.2f (%s)", currentPhase, t.id, t.parentId, t.time, t.steps, t.startingState[0], t.endingState[0], t.startingValue, t.endingValue, success?"success":"failure");
}


////void FFPilotSupervisor::addFFPilotPhaseOutput()
////{
////    // hand the previous FFPilotPhaseOutput message off to the storage list (if this isn't the first or second phase of a simulation stage)
////    if (previousFFPilotPhaseOutputWrapPtr->wrappedMsg()!=NULL)
////    {
////        currentFFPilotPhaseOutputsWrap.AddAllocated(previousFFPilotPhaseOutputWrapPtr->wrappedMsg());
////        previousFFPilotPhaseOutputWrapPtr->setWrappedMsgNull();
////    }

////    // swap the subjects of the current and previous phase output wrapper pointers
////    lm::protowrap::FFPilotPhaseOutputWrap* tmpFFPilotPhaseOutputWrapPtr = previousFFPilotPhaseOutputWrapPtr;
////    previousFFPilotPhaseOutputWrapPtr = currentFFPilotPhaseOutputWrapPtr;
////    currentFFPilotPhaseOutputWrapPtr = tmpFFPilotPhaseOutputWrapPtr;

////    // add a new phase output and set some informational fields
////    lm::ffpilot::io::FFPilotPhaseOutput* newFFPilotPhaseOutputPtr = currentFFPilotPhaseOutputsWrap.Add();
////    newFFPilotPhaseOutputPtr->set_phase_id(currentFFPilotPhaseID());
////    newFFPilotPhaseOutputPtr->set_basin_id(currentStage().basin_id());
////    newFFPilotPhaseOutputPtr->set_tiling_id(currentStage().tiling_id());

////    // set the new phase output to be the current phase output
////    currentFFPilotPhaseOutputWrapPtr->setWrappedMsg(currentFFPilotPhaseOutputsWrap.ReleaseLast());
////}


////bool FFPilotSupervisor::_terminateSimulationPhase()
////{
////    if (not ffpilotPhaseTerminated)
////    {
////        switch(mutableCurrentPhaseLimit()->stop_condition())
////        {
////        case FFPhaseLimEnums::FORWARD_FLUXES:
////            ffpilotPhaseTerminated = (currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->successful_trajectories_launched_count()>=mutableCurrentPhaseLimit()->uvalue());
////            break;
////        case FFPhaseLimEnums::TRAJECTORY_COUNT:
////            // all phases during a forward flux pilot simulation need to record at least one forward crossing or else it can't continue
////            ffpilotPhaseTerminated = ((currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->successful_trajectories_launched_count() > 0) and \
////                                         (currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->successful_trajectories_launched_count() + currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->failed_trajectories_launched_count()>=mutableCurrentPhaseLimit()->uvalue()));
////            break;
////        case FFPhaseLimEnums::TIME:
////            // all phases during a forward flux pilot simulation need to record at least one forward crossing or else it can't continue
////            ffpilotPhaseTerminated = ((currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->successful_trajectories_launched_count() > 0) and \
////                                         (currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->successful_trajectories_launched_total_time() + currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->failed_trajectories_launched_total_time()>=mutableCurrentPhaseLimit()->dvalue()));
////            break;
////        default: THROW_EXCEPTION(UnimplementedException, "unimplemented");
////        }
////    }
////    printFFPilotLimitProgress();

////    return ffpilotPhaseTerminated;
////}

////void FFPilotSupervisor::printFFPilotLimitProgress()
////{
////    if (not ffpilotPhaseTerminated)
////    {
////        // Print some performance statistics, if it has been a while.
////        hrtime currentTime = getHrTime();
////        if (convertHrToSeconds(currentTime-ffpilotProgress_lastPrintTime) > 610.0)
////        {
////            Print::printf(Print::INFO, "Forward Flux Phase Limit Progress");
////            Print::printf(Print::INFO, "  Phase_ID Limit_Type       Progress    Limit");
////            Print::printf(Print::INFO, "----------------------------------------------------");

////            switch(mutableCurrentPhaseLimit()->stop_condition())
////            {
////            case FFPhaseLimEnums::FORWARD_FLUXES:
////                Print::printf(Print::INFO, "%10lld %-17s %12d %12d",
////                    currentFFPilotPhaseID(),
////                    FFPhaseLimEnums::StopCondition_Name(mutableCurrentPhaseLimit()->stop_condition()).c_str(),
////                    currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->successful_trajectories_launched_count(),
////                    mutableCurrentPhaseLimit()->uvalue());
////                break;
////            case FFPhaseLimEnums::TRAJECTORY_COUNT:
////                Print::printf(Print::INFO, "%10lld %-17s %12d %12d",
////                    currentFFPilotPhaseID(),
////                    FFPhaseLimEnums::StopCondition_Name(mutableCurrentPhaseLimit()->stop_condition()).c_str(),
////                    currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->successful_trajectories_launched_count() + currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->failed_trajectories_launched_count(),
////                    mutableCurrentPhaseLimit()->uvalue());
////                break;
////            case FFPhaseLimEnums::TIME:
////                Print::printf(Print::INFO, "%10lld %-17s %12.2e %12.2e",
////                    currentFFPilotPhaseID(),
////                    FFPhaseLimEnums::StopCondition_Name(mutableCurrentPhaseLimit()->stop_condition()).c_str(),
////                    currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->successful_trajectories_launched_total_time() + currentFFPilotPhaseOutputWrapPtr->wrappedMsg()->failed_trajectories_launched_total_time(),
////                    mutableCurrentPhaseLimit()->dvalue());
////                break;
////            default: THROW_EXCEPTION(UnimplementedException, "unimplemented");
////            }

////            ffpilotProgress_lastPrintTime = getHrTime();
////        }
////    }
////}

void FFPilotSupervisor::finishPhase(int stageIndex, int phaseIndex)
{
    StageData& stage = data.getStageData(stageIndex);
    PhaseData& phase = stage.getPhaseData(phaseIndex);

    if (ffpilotPrintStageMessages) Print::printf(Print::INFO, "FFPilot stage %d phase %d finished (trajectories=%d, success=%d).", stageIndex, phaseIndex, phase.trajectoryData.size(), phase.successfulCrossings.size());

    // Make sure that all the crossings finished.
    if (stage.pilot && phase.successfulCrossings.size() != phase.phaseLimit) THROW_EXCEPTION(lm::RuntimeException, "not all crossings finished for the phase: %d,%d,%d,%lu",stageIndex,phaseIndex,phase.successfulCrossings.size(),phase.phaseLimit);
    if (!stage.pilot && phase.successfulCrossings.size()+phase.failedCrossings.size() != phase.phaseLimit) THROW_EXCEPTION(lm::RuntimeException, "not all crossings finished for the phase: %d,%d,%d,%d,%lu",stageIndex,phaseIndex,phase.successfulCrossings.size(),phase.failedCrossings.size(),phase.phaseLimit);

//    for (size_t i=0; i<phase.successfulCrossings.size(); i++)
//    {
//        TrajectoryData t = phase.trajectoryData[phase.successfulCrossings[i]];
//        Print::printf(Print::DEBUG, "FFPilot Crossing %d: %06e,%ld %d->%d",i,t.time,t.steps,t.startingState[0],t.endingState[0]);
//    }

    // Clear any objects used by this phase.
    trajectoryLimits.Clear();
    trajectoryBarriers.Clear();
    outputOptions.Clear();
    if (trajectoryList == NULL) THROW_EXCEPTION(lm::RuntimeException, "Trajectory list was null.")
    delete trajectoryList;
    trajectoryList = NULL;
}

////void FFPilotSupervisor::sendSimulationPhaseOutput()
////{
////    if (not simulationPhaseOutputSent)
////    {
////        // first set the final trajectory id of this phase in the phase output (subtracting one from count since trajectoryList postcrements to get a trajectory_id)
////        currentFFPilotPhaseOutputWrapPtr->set_final_trajectory_id(trajectoryList->count() - 1);

////        // send the phase output to the output writer
////        if ((not currentStage().is_pilot_stage()) or input->getFFPilotOptionsMsg().pilot_stage_output())
////        {
////            if (input->getFFPilotOptionsMsg().phase_output())
////            {
////                // (re)initialize the relevant output options
////                input->reinitOutputOptions(phaseInfo(true), currentStage().is_pilot_stage());

////                // create a handle to the relevant work unit output part
////                lm::message::WorkUnitOutput* wuoPart(ffpilotPhaseOutputContainingMsg.mutable_process_work_unit_output()->mutable_part_output(0));

////                // set the data and output options the work unit output part
////                wuoPart->set_condense_output(input->getOutputOptions().condense_output());
////                wuoPart->set_record_name_prefix(input->getOutputOptions().record_name_prefix());

////                // temporarily hand off the allocated phase output and send it
////                wuoPart->mutable_work_unit_output_generic()->mutable_ffpilot_phase_outputs()->AddAllocated(currentFFPilotPhaseOutputWrapPtr->wrappedMsg());
////                communicator->sendMessage(outputWriterAddress, &ffpilotPhaseOutputContainingMsg);
////                wuoPart->mutable_work_unit_output_generic()->mutable_ffpilot_phase_outputs()->ReleaseLast();
////            }
////        }
////        // even if we're not actually sending phase output, make sure the sent flag is set after sendSimulationPhaseOutput has run. Needed for setting the phase's final trajectory id
////        simulationPhaseOutputSent = true;
////    }
////}

////bool FFPilotSupervisor::incrementSimulationPhase()
////{
////    if (isCurrentPhaseLast())
////    {
////        return false;
////    }
////    else
////    {
////        // increment the currentFFPilotPhase iterator
////        currentFFPilotPhaseIter++;

////        // call the base class method
////        lm::simulation::SimulationSupervisor::incrementSimulationPhase();

////        return true;
////    }
////}

void FFPilotSupervisor::finishStage(int stageIndex)
{
    if (data.getStageData(stageIndex).pilot)
        finishPilotStage(stageIndex);
    else
        finishProductionStage(stageIndex);
}

void FFPilotSupervisor::finishPilotStage(int stageIndex)
{
    StageData& pilotStage = data.getStageData(stageIndex);
    if (ffpilotPrintStageMessages) Print::printf(Print::INFO, "FFPilot stage %d finished (stage_type=Pilot).", stageIndex);

    // Calculate all of the statistics for the phase.
    pilotStage.calculateStageStatistics();

    // Print a summary of the pilots stage.
    if (ffpilotPrintStageMessages)
    {
        Print::printf(Print::INFO, "FFPilot stage %d phase costs are: %s", stageIndex, pilotStage.phaseCosts.toString("%0.2e").c_str());
        Print::printf(Print::INFO, "FFPilot stage %d phase weight sample variances are: %s", stageIndex, pilotStage.phaseWeightVariances.toString("%0.2e").c_str());
        Print::printf(Print::INFO, "FFPilot stage %d conservative estimates of the phase weights are: %s", stageIndex, pilotStage.phaseWeights.toString("%0.2g").c_str());
    }

    // Calculate the optimized trajectory counts.
    ndarray<uint64_t> optimizedTrajectoryCounts = calculateOptimizedTrajectoryCounts(pilotStage);

    // If we have a linked production stage, use the counts to initialize the phase limits.
    if (pilotStage.linkedProductionStage != std::numeric_limits<size_t>::max())
    {
        StageData& prodStage = data.getStageData(pilotStage.linkedProductionStage);
        for (size_t i=0; i<prodStage.getNumberPhases(); i++)
        {
            prodStage.getPhaseData(i).initialize(i==0, prodStage.forward, optimizedTrajectoryCounts[static_cast<uint>(i)]);
        }
    }

    // Print the optimzied counts.
    if (ffpilotPrintStageMessages) Print::printf(Print::INFO, "FFPilot stage %d optimized trajectory counts for error goal %0.2f at confidence %0.2f are: %s", stageIndex, input->ffpilot_options().error_goal(), input->ffpilot_options().error_goal_confidence(), optimizedTrajectoryCounts.toString().c_str());

    // Save the data to the output.
    savePilotStageOutput(stageIndex, pilotStage, optimizedTrajectoryCounts);
}

ndarray<uint64_t> FFPilotSupervisor::calculateOptimizedTrajectoryCounts(StageData& pilotStageData)
{
    // Get some optimization options.
    double errorGoal = input->ffpilot_options().error_goal();
    double errorGoalConfidence = input->ffpilot_options().error_goal_confidence();
    uint64_t minimumCount = input->ffpilot_options().optimize_minimum_trajectories();
    bool minimizeCost = input->ffpilot_options().optimize_minimize_cost();

    // Get the estimated phase cost, weights, and variances.
    ndarray<double>& phaseCosts = pilotStageData.phaseCosts;
    ndarray<double>& phaseWeights = pilotStageData.phaseWeights;
    ndarray<double>& phaseWeightVariances = pilotStageData.phaseWeightVariances;

    // Perform the optimization.
    ndarray<uint64_t> trajectoryCounts;
    if (minimizeCost)
        trajectoryCounts = calculateCountsMinimizeCost(errorGoal, errorGoalConfidence, phaseWeights, phaseWeightVariances, phaseCosts);
    else
        trajectoryCounts = calculateCountsMinimizeCounts(errorGoal, errorGoalConfidence, phaseWeights, phaseWeightVariances);


    // Make sure the counts are lower bounded.
    for (uint i=0; i<trajectoryCounts.shape[0]; i++)
        trajectoryCounts[i] = max(trajectoryCounts[i], minimumCount);

    // See if we need to do any oversampling.
    if (input->ffpilot_options().has_optimize_sampling_multipliers())
    {
        // Get the multipliers.
        ndarray<double> multipliers = robertslab::pbuf::NDArraySerializer::deserialize<double>(input->ffpilot_options().optimize_sampling_multipliers());

        // Make sure the array is the right size.
        if (multipliers.shape.len != 2) THROW_EXCEPTION(lm::RuntimeException, "invalid number of dimensions for optimize_sampling_multipliers: %d (%d)", multipliers.shape.len, 2);
        if (multipliers.shape[0] != 2 || multipliers.shape[1] != trajectoryCounts.shape[0]) THROW_EXCEPTION(lm::RuntimeException, "invalid shape for optimize_sampling_multipliers: %d,%d (%d,%d)",multipliers[0],multipliers[1],2,trajectoryCounts.shape[0]);

        // Apply the multiplier to each phase.
        for (uint i=0; i<trajectoryCounts.shape[0]; i++)
            trajectoryCounts[i] = static_cast<uint64_t>(multipliers[utuple(pilotStageData.forward?0:1,i)]*static_cast<double>(trajectoryCounts[i]));
    }

    return trajectoryCounts;
}


ndarray<uint64_t> FFPilotSupervisor::calculateCountsMinimizeCost(double errorGoal, double errorGoalConfidence, const ndarray<double>& weights, const ndarray<double>& variances, const ndarray<double>& costs)
{
    // Create an array for storing the optimzied counts.
    ndarray<uint64_t> trajectoryCounts(weights.shape[0]);

    // Calculate the per phase factors.
    ndarray<double> phaseFactors = calculatePhaseFactors(weights, variances, costs);
    Print::printf(Print::VERBOSE_DEBUG, "FFPilot phase factors are %s = %0.3e",phaseFactors.toString().c_str(), phaseFactors.sum());

    // Go through each phase and calculate the optimum trajectory count.
    for (uint i=0; i<trajectoryCounts.shape[0]; i++)
    {
        if (i == 0)
        {
            double x = pow(normalZ(errorGoalConfidence)/errorGoal, 2) * sqrt(variances[0]/(pow(weights[0], 2)*costs[0])) * phaseFactors.sum();
            trajectoryCounts[0] = static_cast<uint64_t>(ceil(x));
            Print::printf(Print::VERBOSE_DEBUG, "FFPilot phase %d terms are %0.3e * %0.3e * %0.3e = %0.3e (%ld)",i,pow(normalZ(errorGoalConfidence)/errorGoal, 2),sqrt(variances[0]/(pow(weights[0], 2)*costs[0])),phaseFactors.sum(),x,trajectoryCounts[i]);
        }
        else
        {
            double x = pow(normalZ(errorGoalConfidence)/errorGoal, 2) * sqrt((1.0-weights[i])/(weights[i]*costs[i])) * phaseFactors.sum();
            trajectoryCounts[i] = static_cast<uint64_t>(ceil(x));
            Print::printf(Print::VERBOSE_DEBUG, "FFPilot phase %d terms are %0.3e * %0.3e * %0.3e = %0.3e (%ld)",i,pow(normalZ(errorGoalConfidence)/errorGoal, 2),sqrt(variances[0]/(pow(weights[0], 2)*costs[0])),phaseFactors.sum(),x,trajectoryCounts[i]);
        }
    }

    return trajectoryCounts;
}

ndarray<uint64_t> FFPilotSupervisor::calculateCountsMinimizeCounts(double errorGoal, double errorGoalConfidence, const ndarray<double>& weights, const ndarray<double>& variances)
{
    // Create an array for storing the optimzied counts.
    ndarray<uint64_t> trajectoryCounts(weights.shape[0]);

//    valarray<double> constantFactors(getConstantFactors(weights, variances));

//    constantFactors *= pow(normalZ(errorGoalConfidence)/errorGoal, 2)*(constantFactors.sum());

//    vector<uint64_t> trajectoryCounts;
//    for (int i=0;i<constantFactors.size();i++)
//    {
//        trajectoryCounts.push_back(static_cast<uint64_t>(ceil(constantFactors[i])));
//    }

    return trajectoryCounts;
}

ndarray<double> FFPilotSupervisor::calculatePhaseFactors(const ndarray<double>& weights, const ndarray<double>& variances, const ndarray<double>& costs)
{
    ndarray<double> phaseFactors(weights.shape[0]);
    for (uint i=0; i<phaseFactors.shape[0]; i++)
    {
        if (i == 0)
        {
            // Calculate sqrt(variance*cost/(weight^2)), the phase zero constant factor.
            phaseFactors[0] = sqrt(variances[0]*costs[0]/pow(weights[0], 2));
        }
        else
        {
            // calculate sqrt((1-p)*cost/p), the constant factor in the FFPilot optimizing equation. Doesn't produce the correct value for phase zero
            phaseFactors[i] = sqrt((1.0-weights[i])*costs[i]/weights[i]);
        }
    }
    return phaseFactors;
}

void FFPilotSupervisor::savePilotStageOutput(int stageIndex, StageData& stageData, ndarray<uint64_t>& optimizedCounts)
{
    lm::message::Message msg;
    msg.mutable_process_aggregated_output()->set_record_name_prefix(input->mutable_output_options()->record_name_prefix());
    lm::io::ffpilot::FFPilotStageOutput* output = msg.mutable_process_aggregated_output()->mutable_ffpilot_output()->add_stage_output();

    // Build the output message.
    output->set_id(stageIndex);
    output->set_direction(stageData.forward?lm::io::ffpilot::FFPilotStageOutput::FORWARD:lm::io::ffpilot::FFPilotStageOutput::BACKWARD);
    output->set_type(lm::io::ffpilot::FFPilotStageOutput::PILOT);

    // Add the edges.
    ndarray<double> edges(tilingEdges->shape[0]);
    for (uint i=0; i<edges.shape[0]; i++)
        edges[i] = tilingEdges->get((stageData.forward)?(i):(edges.shape[0]-1-i));
    NDArraySerializer::serializeInto(output->mutable_edges(), edges);

    // Add the counts.
    ndarray<uint64_t> counts(tilingEdges->shape[0]);
    for (uint i=0; i<counts.shape[0]; i++)
        counts[i] = stageData.getPhaseData(static_cast<size_t>(i)).trajectoryData.size();
    NDArraySerializer::serializeInto(output->mutable_trajectory_counts(), counts);

    // Add the costs.
    NDArraySerializer::serializeInto(output->mutable_costs(), stageData.phaseCosts);

    // Add the weights.
    NDArraySerializer::serializeInto(output->mutable_weights(), stageData.phaseWeights);

    // Add the weight variances.
    NDArraySerializer::serializeInto(output->mutable_weight_variances(), stageData.phaseWeightVariances);

    // Add the optimzied trajectory counts.
    NDArraySerializer::serializeInto(output->mutable_optimized_trajectory_counts(), optimizedCounts);

    // Send the output message.
    communicator->sendMessage(outputWriterAddress, &msg);
}

void FFPilotSupervisor::finishProductionStage(int stageIndex)
{
    StageData& prodStage = data.getStageData(stageIndex);
    if (ffpilotPrintStageMessages) Print::printf(Print::INFO, "FFPilot stage %d finished (stage_type=Production).", stageIndex);

    // Calculate all of the statistics for the phase.
    prodStage.calculateStageStatistics();

    // Print a summary of the pilots stage.
    if (ffpilotPrintStageMessages)
    {
        Print::printf(Print::INFO, "FFPilot stage %d phase costs are: %s", stageIndex, prodStage.phaseCosts.toString("%0.2e").c_str());
        Print::printf(Print::INFO, "FFPilot stage %d phase weights are: %s", stageIndex, prodStage.phaseWeights.toString("%0.2g").c_str());
    }

    // Calculate the first passage times.
    ndarray<double> fpts = calculateFirstPassageTimes(prodStage);

    // Print the FPTs.
    if (ffpilotPrintStageMessages)
    {
        Print::printf(Print::INFO, "FFPilot stage %d first passage times are: %s", stageIndex, fpts.toString("%0.2e").c_str());
        Print::printf(Print::INFO, "FFPilot stage %d overall first passage time to the last tile edge is: %0.4e", stageIndex, fpts[fpts.shape[0]-1]);
    }

    // Save the data to the output.
    saveProductionStageOutput(stageIndex, prodStage, fpts);
}

ndarray<double> FFPilotSupervisor::calculateFirstPassageTimes(StageData& stageData)
{
    // Get the phase weights.
    ndarray<double>& weights = stageData.phaseWeights;

    // Calculate the cummulative product of the phase weights for phase >= 1.
    ndarray<double> cumulativeProducts(weights.shape[0]);
    cumulativeProducts[0] = 1.0;
    for (uint i=1; i<weights.shape[0]; i++)
    {
        cumulativeProducts[i] = 1.0;
        for (uint j=1; j<=i; j++)
        {
            cumulativeProducts[i] = cumulativeProducts[i] * weights[j];
        }
    }

    // Create an array for storing the FPTs.
    ndarray<double> fpts(weights.shape[0]);

    // Calculate the FTPs.
    for (uint i=0; i<fpts.shape[0]; i++)
        fpts[i] = weights[0] * (1/cumulativeProducts[i]);

    return fpts;
}

void FFPilotSupervisor::saveProductionStageOutput(int stageIndex, StageData& stageData, ndarray<double>& fpts)
{
    lm::message::Message msg;
    msg.mutable_process_aggregated_output()->set_record_name_prefix(input->mutable_output_options()->record_name_prefix());
    lm::io::ffpilot::FFPilotStageOutput* output = msg.mutable_process_aggregated_output()->mutable_ffpilot_output()->add_stage_output();

    // Build the output message.
    output->set_id(stageIndex);
    output->set_direction(stageData.forward?lm::io::ffpilot::FFPilotStageOutput::FORWARD:lm::io::ffpilot::FFPilotStageOutput::BACKWARD);
    output->set_type(lm::io::ffpilot::FFPilotStageOutput::PRODUCTION);

    // Add the edges.
    ndarray<double> edges(tilingEdges->shape[0]);
    for (uint i=0; i<edges.shape[0]; i++)
        edges[i] = tilingEdges->get((stageData.forward)?(i):(edges.shape[0]-1-i));
    NDArraySerializer::serializeInto(output->mutable_edges(), edges);

    // Add the counts.
    ndarray<uint64_t> counts(tilingEdges->shape[0]);
    for (uint i=0; i<counts.shape[0]; i++)
        counts[i] = stageData.getPhaseData(static_cast<size_t>(i)).trajectoryData.size();
    NDArraySerializer::serializeInto(output->mutable_trajectory_counts(), counts);

    // Add the costs.
    NDArraySerializer::serializeInto(output->mutable_costs(), stageData.phaseCosts);

    // Add the weights.
    NDArraySerializer::serializeInto(output->mutable_weights(), stageData.phaseWeights);

    // Add the weight variances.
    NDArraySerializer::serializeInto(output->mutable_weight_variances(), stageData.phaseWeightVariances);

    // Add the fpts.
    NDArraySerializer::serializeInto(output->mutable_first_passage_times(), fpts);

    // Send the output message.
    communicator->sendMessage(outputWriterAddress, &msg);
}

////    if (not simulationStageOutputSent)
////    {
////        // hand off the final ffpilotPhaseOutputs to the repeated field wrapped by currentFFPilotPhaseOutputsWrap
////        if (previousFFPilotPhaseOutputWrapPtr->wrappedMsg()!=NULL)
////        {
////            currentFFPilotPhaseOutputsWrap.AddAllocated(previousFFPilotPhaseOutputWrapPtr->wrappedMsg());
////            previousFFPilotPhaseOutputWrapPtr->setWrappedMsgNull();
////        }
////        if (currentFFPilotPhaseOutputWrapPtr->wrappedMsg()!=NULL)
////        {
////            currentFFPilotPhaseOutputsWrap.AddAllocated(currentFFPilotPhaseOutputWrapPtr->wrappedMsg());
////            currentFFPilotPhaseOutputWrapPtr->setWrappedMsgNull();
////        }

////        // build the stage output from the stage input and the phase outputs
////        currentFFPilotStageOutputWrap.build(currentStage(), currentFFPilotPhaseOutputsWrap);

////        // add some stage output to the log file
////        if (not currentStage().is_pilot_stage())
////        {
////            // pilot stage log output is handled in .addFFPilotPhaseLimitsFromStageOutput(...)
////            Print::printf(Print::INFO, stageLogProduction(currentFFPilotStageOutputWrap).c_str());
////        }

////        if ((not currentStage().is_pilot_stage()) or input->getFFPilotOptionsMsg().pilot_stage_output())
////        {
////            // (re)initialize the relevant output options
////            input->reinitOutputOptions(currentStageInfo(true), currentStage().is_pilot_stage());

////            if (input->getFFPilotOptionsMsg().stage_output_raw())
////            {
////                // create a handle to the relevant work unit output part
////                lm::message::WorkUnitOutput* wuoPart(ffpilotStageOutputRawContainingMsg.mutable_process_work_unit_output()->mutable_part_output(0));

////                // set the data and output options the work unit output part
////                wuoPart->set_condense_output(input->getOutputOptions().condense_output());
////                wuoPart->set_record_name_prefix(input->getOutputOptions().record_name_prefix());

////                // temporarily hand off the allocated stage output and send it
////                wuoPart->mutable_work_unit_output_generic()->mutable_ffpilot_stage_output_raws()->AddAllocated(currentFFPilotStageOutputWrap.mutable_ffpilot_stage_output_raw()->mutableWrappedMsg());
////                communicator->sendMessage(outputWriterAddress, &ffpilotStageOutputRawContainingMsg);
////                wuoPart->mutable_work_unit_output_generic()->mutable_ffpilot_stage_output_raws()->ReleaseLast();
////                simulationStageOutputSent = true;
////            }
////            if (input->getFFPilotOptionsMsg().stage_output_summary())
////            {
////                // create a handle to the relevant work unit output part
////                lm::message::WorkUnitOutput* wuoPart(ffpilotStageOutputSummaryContainingMsg.mutable_process_work_unit_output()->mutable_part_output(0));

////                // set the data and output options the work unit output part
////                wuoPart->set_condense_output(input->getOutputOptions().condense_output());
////                wuoPart->set_record_name_prefix(input->getOutputOptions().record_name_prefix());

////                // temporarily hand off the allocated stage output and send it
////                wuoPart->mutable_work_unit_output_generic()->mutable_ffpilot_stage_output_summaries()->AddAllocated(currentFFPilotStageOutputWrap.mutable_ffpilot_stage_output_summary()->mutableWrappedMsg());
////                communicator->sendMessage(outputWriterAddress, &ffpilotStageOutputSummaryContainingMsg);
////                wuoPart->mutable_work_unit_output_generic()->mutable_ffpilot_stage_output_summaries()->ReleaseLast();
////                simulationStageOutputSent = true;
////            }
////        }
////    }
////}

////bool FFPilotSupervisor::incrementSimulationStage()
////{
////    if (isCurrentStageLast())
////    {
////        return false;
////    }
////    else
////    {
////        // increment the current stage iterator
////        currentFFPilotStageIter++;
////        _currentStageIndex++;

////        // increment the base phase id counter
////        lm::simulation::SimulationSupervisor::incrementSimulationPhase();

////        return true;
////    }
////}

void FFPilotSupervisor::finishSimulation()
{
    // Call the base class method.
    StagedSimulationSupervisor::finishSimulation();

    Print::printf(Print::INFO, "FFPilot finished in %0.3f seconds.", convertHrToSeconds(getHrTime()-simulationStartTime));
}

////void FFPilotSupervisor::receivedFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg)
////{
////    PROF_BEGIN(PROF_FFLUX_RECEIVED_FINISHED_WORK_UNIT);

////    // deal with the individual parts of the work unit at the ffpilot supervisor level
////    if (currentFFPilotPhaseID()==0)
////    {
////        for (int i=0;i<msg.part_status_size();i++)
////        {
////            receivedFinishedWorkUnitPartPhaseZero(msg.part_status(i));
////        }
////    }
////    else
////    {
////        for (int i=0;i<msg.part_status_size();i++)
////        {
////            receivedFinishedWorkUnitPart(msg.part_status(i));
////        }
////    }

//////    // If the phase "plan" calls for it, generate replacement trajectories
//////    if (currentFFPilotPhaseIter->trajectory_generation()==FFPhaseEnums::LAZY)
//////    {
//////        trajectoryList->initTrajectories(msg.part_status_size());
//////    }

////    // call the base class function
////    lm::simulation::SimulationSupervisor::receivedFinishedWorkUnit(msg);

////    PROF_END(PROF_FFLUX_RECEIVED_FINISHED_WORK_UNIT);
////}

////void FFPilotSupervisor::receivedFinishedWorkUnitPartPhaseZero(const lm::message::WorkUnitStatus& wusMsg)
////{
////    if (not trajectoryList->isTrajectoryAborted(wusMsg.final_state().trajectory_id()))
////    {
////        PROF_BEGIN(PROF_FFLUX_RECEIVED_FINISHED_WORK_UNIT_PHASE_ZERO);

////        // get the relevant Trajectory instance
////        lm::trajectory::Trajectory* trajectory = trajectoryList->getTrajectoryForFinishedWorkUnit(wusMsg.final_state().trajectory_id());

////        // keep track of how much time each phase 0 trajectory spent in the region of a basin other than its initial basin
////        lm::ffpilot::FFPilotPhaseZeroTrajectory* phaseZeroTrajectory = static_cast<lm::ffpilot::FFPilotPhaseZeroTrajectory*>(trajectory);
////        phaseZeroTrajectory->processState(wusMsg.final_state());

////        // update the phase output
////        currentFFPilotPhaseOutputWrapPtr->addEndPointPhaseZero(wusMsg.final_state(), trajectory, input->getFFPilotOptionsMsg().phase_zero_burn_in_count());

////        PROF_END(PROF_FFLUX_RECEIVED_FINISHED_WORK_UNIT_PHASE_ZERO);
////    }
////}

////void FFPilotSupervisor::receivedFinishedWorkUnitPart(const lm::message::WorkUnitStatus& wusMsg)
////{
////    if (not trajectoryList->isTrajectoryAborted(wusMsg.final_state().trajectory_id()) and wusMsg.status()==lm::message::WorkUnitStatus::LIMIT_REACHED)
////    {
////        PROF_BEGIN(PROF_FFLUX_RECEIVED_FINISHED_WORK_UNIT_PHASE_N);

////        // update the phase output
////        currentFFPilotPhaseOutputWrapPtr->addEndPoint(wusMsg.final_state(), *trajectoryList->getTrajectoryForFinishedWorkUnit(wusMsg.final_state().trajectory_id()));

////        PROF_END(PROF_FFLUX_RECEIVED_FINISHED_WORK_UNIT_PHASE_N);
////    }
////}



////std::string FFPilotSupervisor::stageLogPilot(const lm::protowrap::FFPilotStageOutputWrap& stageOutput, double errorGoal, double errorGoalConfidence, vector<uint64_t>& trajectoryCounts) const
////{
////    stringstream stageLog;
////    stageLog.unsetf(std::ios::floatfield);                  // allow for dynamic choice between float and sci format
////    //stageLog.setf(std::ios::fixed, std::ios::floatfield); // force float format
////    stageLog.precision(3);

////    stageLog << "Pilot stage output:\n";
////    stageLog << "The phase costs are:\n" << stageOutput.ffpilot_stage_output_summary().costs() << "\n";
////    stageLog << "The phase weight sample variances are:\n" << stageOutput.ffpilot_stage_output_raw().variances() << "\n";
////    stageLog << "Conservative estimates of the phase weights are:\n" << estimateBernoulliProbabilities(stageOutput) << "\n";
////    stageLog << "Attempting to achieve error goal " << errorGoal << " (confidence level " << errorGoalConfidence << ") with the following optimized trajectory counts:\n" << trajectoryCounts << "\n";

////    return stageLog.str();
////}

////std::string FFPilotSupervisor::stageLogProduction(const lm::protowrap::FFPilotStageOutputWrap& stageOutput) const
////{
////    stringstream stageLog;
////    stageLog.unsetf(std::ios::floatfield);
////    stageLog.precision(3);

////    stageLog << "Production stage output:\n";
////    stageLog << "The phase costs are:\n" << stageOutput.ffpilot_stage_output_summary().costs() << "\n";
////    stageLog << "The phase weights are:\n" << stageOutput.ffpilot_stage_output_summary().weights() << "\n";
////    stageLog << "The first passage times to each tile edge are:\n" << stageOutput.ffpilot_stage_output_summary().first_passage_times() << "\n";
////    stageLog << "The overall first passage time from the starting basin to the last tile edge is:\n" << stageOutput.ffpilot_stage_output_summary().first_passage_times(stageOutput.ffpilot_stage_output_summary().first_passage_times_size() - 1) << "\n";

////    return stageLog.str();
////}

////// setters
////void FFPilotSupervisor::setInput(lm::input::Input* newInput)
////{
////    lm::simulation::SimulationSupervisor::setInput(newInput);
////    input = static_cast<lm::ffpilot::input::FFPilotInput*>(lm::simulation::SimulationSupervisor::input);
////}

void FFPilotSupervisor::printPerformanceStatistics()
{
    SimulationSupervisor::printPerformanceStatistics();
    Print::printf(Print::INFO, "FFPilot status (%s): stage %d of %d, phase %d of %d", lm::message::Communicator::printableAddress(communicator->getSourceAddress()).c_str(), currentStage, getNumberStages(), currentPhase, getNumberPhases(currentStage));
    Print::printf(Print::INFO, "%5s %5s %9s %10s %10s %10s %10s", "Stage", "Phase", "Goal", "Successful", "Failed", "Total", "Limit");
    Print::printf(Print::INFO, "-------------------------------------------");
    for (size_t i=0; i<data.getNumberStages(); i++)
    {
        StageData& stageData = data.getStageData(i);
        for (size_t j=0; j<stageData.getNumberPhases(); j++)
        {
            PhaseData& p = stageData.getPhaseData(j);
            Print::printf(Print::INFO, "%5ld %5ld %9.2e %10ld %10ld %10ld %10ld", i, j, p.goalEdge, p.successfulCrossings.size(), p.failedCrossings.size(), p.trajectoryData.size(), p.phaseLimit);
        }
    }
}

}
}



/*
 * sample optimziation data
 *

phaseCosts[0] = 4.62;
phaseCosts[1] = 0.151;
phaseCosts[2] = 1.61;
phaseCosts[3] = 2.41;
phaseCosts[4] = 2.72;
phaseCosts[5] = 2.09;
phaseCosts[6] = 1.65;
phaseCosts[7] = 1.22;
phaseCosts[8] = 0.967;
phaseCosts[9] = 1.01;
phaseCosts[10] = 0.951;
phaseCosts[11] = 1.25;
phaseCosts[12] = 1.5;
phaseWeightVariances[0] = 152;
phaseWeightVariances[1] = 0.018;
phaseWeightVariances[2] = 0.12;
phaseWeightVariances[3] = 0.228;
phaseWeightVariances[4] = 0.245;
phaseWeightVariances[5] = 0.155;
phaseWeightVariances[6] = 0.0665;
phaseWeightVariances[7] = 0.0211;
phaseWeightVariances[8] = 0.00473;
phaseWeightVariances[9] = 0.00296;
phaseWeightVariances[10] = 0.000693;
phaseWeightVariances[11] = 0.000495;
phaseWeightVariances[12] = 0.000396;
phaseWeights[0] = 4.62;
phaseWeights[1] = 0.0179;
phaseWeights[2] = 0.137;
phaseWeights[3] = 0.345;
phaseWeights[4] = 0.56;
phaseWeights[5] = 0.8;
phaseWeights[6] = 0.922;
phaseWeights[7] = 0.974;
phaseWeights[8] = 0.993;
phaseWeights[9] = 0.995;
phaseWeights[10] = 0.998;
phaseWeights[11] = 0.998;
phaseWeights[12] = 0.999;

//optimized trajectory counts for error goal 0.05 at confidence 0.95 are: [32361,496845,51558,23135,14009,9015,5902,3856,2226,1839,1197,1044,1000]

*/
