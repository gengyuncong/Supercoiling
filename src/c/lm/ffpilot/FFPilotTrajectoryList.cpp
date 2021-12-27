/*
 * Copyright 2012-2018 Johns Hopkins University
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

//#include <algorithm>
//#include <cmath>
//#include <csignal>
//#include <functional>
//#include <list>
//#include <map>
//#include <numeric>
//#include <string>
//#include <vector>

//#include "lm/EnumHelper.h"
//#include "lm/Exceptions.h"
#include "lm/ffpilot/FFPilotTrajectoryList.h"
//#include "lm/ffpilot/input/FFPilotPhase.pb.h"
//#include "lm/io/CMEState.pb.h"
//#include "lm/io/FirstPassageTimes.pb.h"
//#include "lm/io/hdf5/SimulationFile.h"
//#include "lm/input/ReactionModel.pb.h"
//#include "lm/io/SpeciesCounts.pb.h"
//#include "lm/io/SpeciesTimeSeries.pb.h"
//#include "lm/io/TrajectoryState.pb.h"
//#include "lm/main/Globals.h"
//#include "lm/message/WorkUnitStatus.pb.h"
//#include "lm/Print.h"
//#include "lm/ffpilot/io/FFPilotPhaseOutputWrap.h"
//#include "lm/protowrap/Repeated.h"
//#include "lm/trajectory/Trajectory.h"
//#include "lm/tiling/Tilings.h"
//#include "lptf/Profile.h"
//#include "lptf/ProfileCodes.h"

//using lm::protowrap::FFPilotPhaseOutputWrap;
//using lm::protowrap::Repeated;
//using std::map;
//using std::string;
//using std::vector;

namespace lm {
namespace ffpilot {

//// ffpilotPhase n==0 constructor
//FFPilotTrajectoryList::FFPilotTrajectoryList(
//    uint64_t count, uint64_t newSimulationPhaseID, const FFPilotPhaseMsg& ffpilotPhase,
//    uint simultaneousWorkUnits, const FFPilotInput& input, const lm::input::Basin& basin
//)
//:TrajectoryList(count, newSimulationPhaseID),input(input),ffpilotPhase(ffpilotPhase),
// previousPhaseOutputPtr(NULL),cyclicCounter(0),simultaneousWorkUnits(simultaneousWorkUnits)
//{
//    // consistency check
//    if (ffpilotPhase.phase_id()!=0)
//    {
//        throw ConsistencyException("Forward Flux phase 0 version of FFPilotTrajectoryList constructor called durring phase %d", ffpilotPhase.phase_id());
//    }

//    // figure out how many trajectories we need to start right now
//    uint64_t trajectoriesToStart = getTrajectoriesToStart();

//    // initialize trajectories based on simulation input files
//    for (uint64_t i=0;i<trajectoriesToStart;i++)
//    {
//        initFFPilotPhaseZeroTrajectory(input, basin.species_count().begin(), basin.species_count().end(), 0.0, simulationPhaseID(), DEFAULT_TRAJECTORY_ID);
//    }
//}

//// ffpilotPhase n>0 constructor
//FFPilotTrajectoryList::FFPilotTrajectoryList(
//    uint64_t count, uint64_t newSimulationPhaseID, const FFPilotPhaseMsg& ffpilotPhase,
//    uint simultaneousWorkUnits, const FFPilotInput& input,
//    const FFPilotPhaseOutputWrap& previousPhaseOutput
//)
//:TrajectoryList(count, newSimulationPhaseID),input(input),ffpilotPhase(ffpilotPhase),
// previousPhaseOutputPtr(&previousPhaseOutput),cyclicCounter(0),simultaneousWorkUnits(simultaneousWorkUnits)
//{
//    // consistency check
//    if (ffpilotPhase.phase_id()==0) throw ConsistencyException("Forward Flux phase n>0 version of FFPilotTrajectoryList constructor called durring phase 0. phase_id: %d", ffpilotPhase.phase_id());

//    // figure out how many trajectories we need to start right now
//    uint64_t trajectoriesToStart = getTrajectoriesToStart();

//    initTrajectories(trajectoriesToStart);
//}

//// ffpilotPhase custom constructor
//FFPilotTrajectoryList::FFPilotTrajectoryList(
//    uint64_t count, uint64_t newSimulationPhaseID, const FFPilotPhaseMsg& ffpilotPhase,
//    uint simultaneousWorkUnits, const FFPilotInput& input
//)
//:TrajectoryList(count, newSimulationPhaseID),input(input),ffpilotPhase(ffpilotPhase),
// previousPhaseOutputPtr(&previousPhaseOutputCustomWrap), cyclicCounter(0),simultaneousWorkUnits(simultaneousWorkUnits)
//{
//    previousPhaseOutputCustom.mutable_successful_trajectory_end_points()->CopyFrom(ffpilotPhase.start_points());
//    previousPhaseOutputCustomWrap.setWrappedMsg(&previousPhaseOutputCustom);

//    // consistency check
//    if (ffpilotPhase.phase_id()==0)
//    {
//        throw ConsistencyException("Forward Flux phase custom version of FFPilotTrajectoryList constructor called durring phase 0. phase_id: %d", ffpilotPhase.phase_id());
//    }

//    // figure out how many trajectories we need to start right now
//    uint64_t trajectoriesToStart = getTrajectoriesToStart();

//    initTrajectories(trajectoriesToStart);
//}

//uint64_t FFPilotTrajectoryList::getSimultaneousTrajectories()
//{
//    if (ffpilotPhase.phase_id()==0)
//    {
//        if (input.getFFPilotOptionsMsg().simultaneous_trajectories_phase_zero() != 0)
//        {
//            return input.getFFPilotOptionsMsg().simultaneous_trajectories_phase_zero();
//        }
//        else
//        {
//            return simultaneousWorkUnits * ffpilotPhase.options().parts_per_work_unit();
//        }
//    }
//    else
//    {
//        if (input.getFFPilotOptionsMsg().simultaneous_trajectories() != 0)
//        {
//            return input.getFFPilotOptionsMsg().simultaneous_trajectories();
//        }
//        else
//        {
//            return simultaneousWorkUnits * ffpilotPhase.options().parts_per_work_unit();
//        }
//    }
//}

//uint64_t FFPilotTrajectoryList::getTrajectoriesToStart()
//{
//    const FFPilotPhaseLimitMsg& ffpilotPhaseLimit = ffpilotPhase.ffpilot_phase_limit();
    
//    if (ffpilotPhase.trajectory_generation()==FFPhaseEnums::EAGER)
//    {
//        // EAGER is only implemented for certain ffpilotPhaseLimit.stop_condition() values
//        if (ffpilotPhaseLimit.stop_condition()==FFPhaseLimEnums::TRAJECTORY_COUNT or (ffpilotPhaseLimit.stop_condition()==FFPhaseLimEnums::FORWARD_FLUXES and ffpilotPhase.phase_id()==0))
//        {
//            if (ffpilotPhaseLimit.has_events_per_trajectory())
//            {
//                // given that our trajectory limits are set up to observe x events per trajectory, run ceil(y/x) trajectories to ensure that we observe at least y events total
//                return (uint64_t)ceil(ffpilotPhaseLimit.uvalue()/(double)ffpilotPhaseLimit.events_per_trajectory());
//            }
//            else
//            {
//                // in this case assume that we want to observe the maximum number of events per trajectory, so just run enough trajectories for one "round" (ie one trajectory per work unit runner)
//                return getSimultaneousTrajectories();
//            }
//        }
//        else THROW_EXCEPTION(UnimplementedException, "In Forward Flux phase %d, ffpilotPhase.trajectory_generation()==EAGER is only implemented for certain ffpilotPhaseLimit.stop_condition() values (ie those that let us calculate the necessary trajectory count up front). Attempting to use unimplemented ffpilotPhaseLimit.stop_condition(): %s", ffpilotPhase.phase_id(), FFPhaseLimEnums::StopCondition_Name(ffpilotPhaseLimit.stop_condition()).c_str());
//    }
//    else if (ffpilotPhase.trajectory_generation()==FFPhaseEnums::LAZY)
//    {
//        return getSimultaneousTrajectories();
//    }
//    else THROW_EXCEPTION(UnimplementedException, "unimplemented");
//}

//void FFPilotTrajectoryList::workUnitPartFinished(const message::WorkUnitStatus& wusMsg, lm::trajectory::Trajectory* traj)
//{
//    // Call the base class method.
//    lm::trajectory::TrajectoryList::workUnitPartFinished(wusMsg, traj);

//    // If the work unit stopped because it hit a terminating limit...
//    if (wusMsg.status()==lm::message::WorkUnitStatus::LIMIT_REACHED)
//    {
//        // FFPilotSupervisor will have already extracted the necessary information, so we need to dispose of the trajectory here
//        if (ffpilotPhase.trajectory_generation()==FFPhaseEnums::LAZY)
//        {
//            // If the phase "plan" calls for it, generate a replacement trajectory by recycling the old one
//            recycleFFPilotTrajectory(traj);
//        }
//        else
//        {
//            // otherwise just delete the trajectory
//            deleteTrajectory(traj->getID());
//        }
//    }
//    else if (ffpilotPhase.phase_id()==0)
//    {
//        // clear out any limit tracking time series data. prevents a major slowdown on long runs
//        limitTrackingListWrap.setWrappedMsg(traj->getStateMutable()->mutable_limit_tracking_list());
//        limitTrackingListWrap.clearTimeSeriesData();
//    }
//}

//void FFPilotTrajectoryList::initTrajectories(uint64_t trajectoriesToStart)
//{
//    switch(ffpilotPhase.trajectory_duplication())
//    {
//    case FFPhaseEnums::NONE:
//        // do nothing
//        break;
//    case FFPhaseEnums::CYCLIC:
//        // initialize trajectories by cycling through EndPoints from a previous phase
//        initTrajectoriesCyclic(trajectoriesToStart);
//        break;
//    case FFPhaseEnums::UNIFORM_RANDOM:
//        // initialize trajectories based on randomly selected EndPoints from a previous phase
//        initTrajectoriesUniformRandom(trajectoriesToStart);
//        break;
//    default: THROW_EXCEPTION(UnimplementedException, "Unimplemented");
//    }
//}

//void FFPilotTrajectoryList::recycleFFPilotTrajectory(lm::trajectory::Trajectory* traj)
//{
//    limitTrackingListWrap.setWrappedMsg(traj->getStateMutable()->mutable_limit_tracking_list());
//    limitTrackingListWrap.clearStateData();

//    traj->clearLimitReached();

//    switch(ffpilotPhase.trajectory_duplication())
//    {
//    case FFPhaseEnums::NONE:
//        // do nothing
//        break;
//    case FFPhaseEnums::CYCLIC:
//        // initialize trajectories by cycling through EndPoints from a previous phase
//        recycleTrajectoryCyclic(traj->getID());
//        break;
//    case FFPhaseEnums::UNIFORM_RANDOM:
//        // initialize trajectories based on randomly selected EndPoints from a previous phase
//        recycleTrajectoryUniformRandom(traj->getID());
//        break;
//    default: THROW_EXCEPTION(UnimplementedException, "Unimplemented");
//    }
//}

//void FFPilotTrajectoryList::initTrajectoriesCyclic(uint64_t trajectoriesToStart)
//{
//    for (uint64_t i=0;i<trajectoriesToStart;i++)
//    {
//        const lm::protowrap::EndPointVector::Pair& endPointPair(previousPhaseOutputPtr->getEndPointCyclic(cyclicCounter++));
//        initFFPilotTrajectory(input, endPointPair.first->species_coordinates().begin(), endPointPair.first->species_coordinates().end(), endPointPair.first->times(endPointPair.second), simulationPhaseID(), DEFAULT_TRAJECTORY_ID);
//    }
//}

//void FFPilotTrajectoryList::initTrajectoriesUniformRandom(uint64_t trajectoriesToStart)
//{
//    for (uint64_t i=0;i<trajectoriesToStart;i++)
//    {
//        PROF_BEGIN(PROF_TRAJECTORY_LIST_WORK_UNIT_FINISHED_ONE);
//        const lm::protowrap::EndPointVector::Pair& endPointPair(previousPhaseOutputPtr->getEndPointUniformRandom());
//        PROF_END(PROF_TRAJECTORY_LIST_WORK_UNIT_FINISHED_ONE);

//        PROF_BEGIN(PROF_TRAJECTORY_LIST_WORK_UNIT_FINISHED_TWO);
//        initFFPilotTrajectory(input, endPointPair.first->species_coordinates().begin(), endPointPair.first->species_coordinates().end(), endPointPair.first->times(endPointPair.second), simulationPhaseID(), DEFAULT_TRAJECTORY_ID);
//        PROF_END(PROF_TRAJECTORY_LIST_WORK_UNIT_FINISHED_TWO);
//    }
//}

//lm::trajectory::Trajectory* FFPilotTrajectoryList::recycleTrajectoryCyclic(uint64_t oldID)
//{
//    // choose an endpoint (by cycling through the list of endpoints) from the previous phase to use as a starting state
//    const lm::protowrap::EndPointVector::Pair& endPointPair(previousPhaseOutputPtr->getEndPointCyclic(cyclicCounter++));

//    // reuse as much of the existing trajectory as possible
//    lm::trajectory::Trajectory* recycTraj = recycleTrajectory(endPointPair.first->species_coordinates().begin(), endPointPair.first->species_coordinates().end(), endPointPair.first->times(endPointPair.second), oldID, DEFAULT_TRAJECTORY_ID);

//    // record the correct "initial" state in the recycled trajectory
//    static_cast<lm::ffpilot::FFPilotTrajectory*>(recycTraj)->setInitialStateToCurrentState();

//    // return the refurbished trajectory
//    return recycTraj;
//}
//lm::trajectory::Trajectory* FFPilotTrajectoryList::recycleTrajectoryUniformRandom(uint64_t oldID)
//{
//    // choose an endpoint (at random from the list of endpoints) from the previous phase to use as a starting state
//    const lm::protowrap::EndPointVector::Pair& endPointPair(previousPhaseOutputPtr->getEndPointUniformRandom());

//    // reuse as much of the existing trajectory as possible
//    lm::trajectory::Trajectory* recycTraj = recycleTrajectory(endPointPair.first->species_coordinates().begin(), endPointPair.first->species_coordinates().end(), endPointPair.first->times(endPointPair.second), oldID, DEFAULT_TRAJECTORY_ID);

//    // record the correct "initial" state in the recycled trajectory
//    static_cast<lm::ffpilot::FFPilotTrajectory*>(recycTraj)->setInitialStateToCurrentState();

//    // return the refurbished trajectory
//    return recycTraj;
//}

}
}
