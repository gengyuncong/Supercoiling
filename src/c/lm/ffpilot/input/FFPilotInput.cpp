///*
// * Copyright 2012-2019 Johns Hopkins University
// *
// * Licensed under the Apache License, Version 2.0 (the "License");
// * you may not use this file except in compliance with the License.
// * You may obtain a copy of the License at
// *
// *     http://www.apache.org/licenses/LICENSE-2.0
// *
// * Unless required by applicable law or agreed to in writing, software
// * distributed under the License is distributed on an "AS IS" BASIS,
// * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// * See the License for the specific language governing permissions and
// * limitations under the License.
// *
// * Developed by: Roberts Group
// *               Johns Hopkins University
// *               http://biophysics.jhu.edu/roberts/
// *
// * Author(s): Elijah Roberts, Max Klein
// */

//#include <cmath>
//#include <limits>
//#include <map>
//#include <string>
//#include <vector>

//#include "lm/ClassFactory.h"
//#include "lm/EnumHelper.h"
//#include "lm/ffpilot/FFPilotMath.h"
//#include "lm/ffpilot/input/FFPilotInput.h"
//#include "lm/ffpilot/input/FFPilotOptions.pb.h"
//#include "lm/io/hdf5/SimulationFile.h"
//#include "lm/input/Options.pb.h"
//#include "lm/input/OutputOptions.pb.h"
//#include "lm/input/TrajectoryLimits.pb.h"
//#include "lm/main/Globals.h"
//#include "lm/Print.h"
//#include "lm/limit/TrajectoryLimits.h"
//#include "lm/Types.h"

//using std::string;
//using std::vector;

//using lm::input::Options;
//using lm::input::OutputOptions;
//using lm::ffpilot::input::FFPilotOptions;
//using lm::limit::LimitElement;

//namespace lm {
//namespace ffpilot {
//namespace input {

//bool FFPilotInput::registered=FFPilotInput::registerClass();

//bool FFPilotInput::registerClass()
//{
//    lm::ClassFactory::getInstance().registerClass("lm::input::Input","lm::ffpilot::input::FFPilotInput",(ClassAllocator)&FFPilotInput::allocateObject);
//    return true;
//}

//void* FFPilotInput::allocateObject()
//{
//    return new FFPilotInput();
//}

//FFPilotInput::FFPilotInput(): Input()
//{
//}

//void FFPilotInput::init()
//{
//    // Call the base class method.
//    lm::input::Input::init();

//    // initialize FFPilot specific settings
//    initFFPilotOptions();

//    // although reinitOutputOptions() will be run at least once more before any related values are actually used, run it once here so the sanity check works correctly
//    reinitOutputOptions("", false);

////    // warn the user about any unrecognized/unparsed simulation parameters
////    initSanityCheck();

////    // set up the initial simulation "plan", in the form of a sequence of stages
////    initFFPilotStages();
//}

//// Get the Forward Flux specific options.
//void FFPilotInput::initFFPilotOptions()
//{
////    FFPilotOptions* ffpopt = inputMsg.mutable_ffpilot_options();
    
////    parseAndSet("errorGoal",           &FFPilotOptions::set_error_goal,            ffpopt);
////    parseAndSet("errorGoalConfidence", &FFPilotOptions::set_error_goal_confidence, ffpopt);

////    parseAndSet("productionStageCountMinimum", &FFPilotOptions::set_production_stage_count_minimum, ffpopt);
////    parseAndSet("phaseZeroBurnInCount",        &FFPilotOptions::set_phase_zero_burn_in_count,       ffpopt);

////    parseAndSet("ffpilotPilotOutput",        &FFPilotOptions::set_pilot_stage_output,   ffpopt);
////    parseAndSet("ffpilotPhaseOutput",        &FFPilotOptions::set_phase_output,         ffpopt);
////    parseAndSet("ffpilotStageOutputRaw",     &FFPilotOptions::set_stage_output_raw,     ffpopt);
////    parseAndSet("ffpilotStageOutputSummary", &FFPilotOptions::set_stage_output_summary, ffpopt);

////    parseAndSet("ffpilotMinimizeCost", &FFPilotOptions::set_minimize_cost, ffpopt);

////    parseAndSet("simultaneousTrajectoriesPhaseZero", &FFPilotOptions::set_simultaneous_trajectories_phase_zero, ffpopt);
////    parseAndSet("simultaneousTrajectories",          &FFPilotOptions::set_simultaneous_trajectories,            ffpopt);

////    // TODO: figure out how to properly control landscape error and remove this
////    parseAndSet("phaseZeroSamplingMultiplier", &FFPilotOptions::set_phase_zero_sampling_multiplier, ffpopt);

////    // determine the default value of pilot_stage_count based on confidence and twice the specified error goal
////    uint64_t pscd = npilotConservative(2*getFFPilotOptionsMsg().error_goal(), getFFPilotOptionsMsg().error_goal_confidence());
////    parseAndSet("pilotStageCount", &FFPilotOptions::set_pilot_stage_count, ffpopt, &pscd);
//}

//void FFPilotInput::initSanityCheck()
//{
//    // run the base class sanity check (covers parameter parsing)
//    Input::initSanityCheck();
    
//    // ensure that we got a tiling
//    if (getTilings().size() <= 0)
//    {
//        THROW_EXCEPTION(InputException, "FFPilot simulation requested (via cmd line arguments), but no tilings were provided in the input.");
//    }

//    // ensure that the tiling has at least one basin
//    int totalBasinCount = 0;
//    for (lm::tiling::Tilings::const_iterator tilingIt=getTilings().begin();tilingIt!=getTilings().end();++tilingIt)
//    {
//        totalBasinCount += tilingIt->second->basins().size();
//    }
//    if (totalBasinCount <= 0)
//    {
//        THROW_EXCEPTION(InputException, "FFPilot simulation requested (via cmd line arguments), but there were no basins in any of the tilings provided in the input.");
//    }
//}

////void FFPilotInput::initFFPilotStages()
////{
////    if (hasFFPilotStages())
////    {
////        initFFPilotStagesCustom();
////    }
////    else
////    {
////        initSimulationStageList();
////    }
////}
////
////void FFPilotInput::initFFPilotStagesCustom()
////{
////
////}
////
////void FFPilotInput::initFFPilotStagesFromScratch()
////{
////
////}

//void FFPilotInput::reinitOptions(int64_t phaseID)
//{
////    Options* opt = inputMsg.mutable_options();
////    opt->Clear();

////    if (phaseID==0)
////    {
////        // phase zero always uses 1 part per work unit
////        opt->set_parts_per_work_unit(1);
////        parseAndSet("stepsPerWorkUnitPart", &Options::set_steps_per_work_unit_part, opt);
////    }
////    else
////    {
////        uint64_t defaultPartsPerWorkUnit = partsPerWorkUnit > 0 ? partsPerWorkUnit : 100;
////        parseAndSet("partsPerWorkUnit", &Options::set_parts_per_work_unit, opt, &defaultPartsPerWorkUnit);

////        parseAndSet("stepsPerWorkUnitPart", &Options::set_steps_per_work_unit_part, opt);
////    }
//}

//void FFPilotInput::reinitOutputOptions(const std::string& recordNamePrefix, bool isPilotStage)
//{
////    OutputOptions* oopt = inputMsg.mutable_output_options();
////    oopt->Clear();

////    // set output options for the pilot stage only if pilot stage output is explicitly requested
////    if ((not isPilotStage) or getFFPilotOptionsMsg().pilot_stage_output())
////    {
////        initOutputOptions(recordNamePrefix);

////        bool condenseOutputDefault = true;
////        parseAndSet("condenseOutput", &OutputOptions::set_condense_output, oopt, &condenseOutputDefault);

////        // Flags that control whether output is recorded for the initial and/or the final state of every trajectory.
////        bool writeInitialFinalStateDefault = false;
////        parseAndSet("writeInitialTrajectoryState", &OutputOptions::set_write_initial_trajectory_state, oopt, &writeInitialFinalStateDefault);
////        parseAndSet("writeFinalTrajectoryState",   &OutputOptions::set_write_final_trajectory_state,   oopt, &writeInitialFinalStateDefault);

////        // If output of initial or final state has been requested but none of the output intervals have been set, assume the user wants just the species counts output from just the initial and final states
////        if (getOutputOptions().write_initial_trajectory_state() or getOutputOptions().write_final_trajectory_state())
////        {
////            if (not (getOutputOptions().has_degree_advancement_write_interval() or getOutputOptions().has_lattice_write_interval() or getOutputOptions().has_order_parameter_write_interval() or getOutputOptions().has_species_write_interval()))
////            {
////                // set species_write_interval to the largest value possible. This ensures that only the first and/or last moments of a trajectory will be written out
////                oopt->set_species_write_interval(std::numeric_limits<double>::max());
////            }
////        }
////    }
////    else
////    {
////        // no output requested, make sure that it's all turned off
////        oopt->set_write_initial_trajectory_state(false);
////        oopt->set_write_final_trajectory_state(false);
////    }
//}

//void FFPilotInput::reinitTrajectoryLimits(const FFPilotPhaseMsg& ffpilotPhase, const lm::tiling::Tiling& tiling)
//{
//    // get ref to phase limit
//    const FFPilotPhaseLimitMsg& ffpilotPhaseLimit = ffpilotPhase.ffpilot_phase_limit();

//    if (ffpilotPhase.phase_id()==0)
//    {
//        reinitTrajectoryLimitsPhaseZero(ffpilotPhase, tiling);
//    }
//    else
//    {
//        trajectoryLimits.Clear();
//        limitTrackingListWrap.Clear();

//        // - if currentFFPilotPhaseID() > 0, we can use addTileExitLimitsMsg() in a straightforward way to set the needed limits. Two limits are set:
//        //     - if limit id==0 is triggered, this indicates that the trajectory fluxed backwards
//        //     - if limit id==1 is triggered, this indicates that the trajectory fluxed forwards
//        trajectoryLimits.addTileExitLimitsMsg(tiling, 0, ffpilotPhase.phase_id());
//        limitTrackingListWrap.addTrackingMsg(trajectoryLimits.findMsg(0), true, true, ffpilotPhaseLimit.events_per_trajectory());
//        limitTrackingListWrap.addTrackingMsg(trajectoryLimits.findMsg(1), true, true, ffpilotPhaseLimit.events_per_trajectory());
//    }
//}

//void FFPilotInput::reinitTrajectoryLimitsPhaseZero(const FFPilotPhaseMsg& ffpilotPhase, const lm::tiling::Tiling& tiling)
//{
//    // get ref to phase limit
//    const FFPilotPhaseLimitMsg& ffpilotPhaseLimit = ffpilotPhase.ffpilot_phase_limit();

//    trajectoryLimits.Clear();
//    limitTrackingListWrap.Clear();

//    // - first we set a limit with id==0
//    //     - this limit is the important one. a triggering of this limit corresponds to one of the flux events that we're trying to sample during phase 0
//    trajectoryLimits.addTileExitLimitsMsg(tiling, -1, 0, false, true);
//    if (ffpilotPhaseLimit.events_per_trajectory() == std::numeric_limits<uint64_t>::max())
//    {
//        limitTrackingListWrap.addTrackingMsgNonterminating(trajectoryLimits.findMsg(0), true, true);
//    }
//    else
//    {
//        limitTrackingListWrap.addTrackingMsg(trajectoryLimits.findMsg(0), true, true, ffpilotPhaseLimit.events_per_trajectory());
//    }

//    // - next, we set two more limits with id==1 and id==2
//    //     - these limits are used to help track which basin was last visited by a trajectory
//    //     - limit_id==1: tracks flux back into the starting basin
//    //     - limit_id==2: tracks flux into the basin opposite from the starting basin
//    trajectoryLimits.addTileExitLimitsMsg(tiling, 0, tiling.edges().lastIndex());
//    limitTrackingListWrap.addTrackingMsgNonterminating(trajectoryLimits.findMsg(1), true, true);
//    limitTrackingListWrap.addTrackingMsgNonterminating(trajectoryLimits.findMsg(2), true, true);
//}

//}
//}
//}
