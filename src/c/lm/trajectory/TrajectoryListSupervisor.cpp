/*
 * Copyright 2018 Johns Hopkins University
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

#include <algorithm>

#include "lm/simulation/SimulationSupervisor.h"
#include "lm/trajectory/TrajectoryListSupervisor.h"
#include "lm/trajectory/TrajectoryList.h"
#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"

namespace lm {
namespace trajectory {

TrajectoryListSupervisor::TrajectoryListSupervisor()
:trajectoryList(NULL)
{
}

TrajectoryListSupervisor::~TrajectoryListSupervisor()
{
    if (trajectoryList != NULL) delete trajectoryList; trajectoryList = NULL;
}

void TrajectoryListSupervisor::startSimulation()
{
    // Build the trajectory list.
    buildTrajectoryList();

    // Build the trajectory options.
    buildTrajectoryOptions();

    // Make sure we got a trajectory list.
    if (trajectoryList == NULL) THROW_EXCEPTION(lm::RuntimeException, "No trajectory list created for the trajectory simulation supervisor.");

    // Call the base class method.
    SimulationSupervisor::startSimulation();
}

bool TrajectoryListSupervisor::assignWork()
{
    // While we have free slots and waiting trajectories, assign them.
    while (slotList->hasFreeSlots() && trajectoryList->anyWaiting())
    {
        // Create the run work unit message.
        lm::message::Message msg;
        lm::message::RunWorkUnit* rwuMsg = msg.mutable_run_work_unit();

        // Build the run work units message.
        buildRunWorkUnitHeader(rwuMsg);

        // Get the free slot.
        const lm::slot::Slot slot = slotList->getFreeSlot();

        // Build the work unit parts.
        buildRunWorkUnitParts(rwuMsg, slot.getSimultaneousWorkUnits());

        //rwuMsg->PrintDebugString();  //DEBUGTODO

        // Run the work unit.
        slotList->runWorkUnit(communicator, &msg);
    }

    return !trajectoryList->allFinished();
}

void TrajectoryListSupervisor::buildRunWorkUnitHeader(lm::message::RunWorkUnit* msg)
{
    // Call the base class method.
    SimulationSupervisor::buildRunWorkUnitHeader(msg);

    // Set the solver to be master equation.
    msg->set_solver_type(lm::types::SolverType::ME);

    // Set the limits.
    if (trajectoryLimits.limits().size() > 0)
        msg->mutable_trajectory_limits()->CopyFrom(trajectoryLimits);

    // Set the barriers.
    if (trajectoryBarriers.barrier().size() > 0)
        msg->mutable_trajectory_barriers()->CopyFrom(trajectoryBarriers);

    // Set the output options.
    msg->mutable_output_options()->CopyFrom(input->output_options());

    // Set the maximum number of steps for the work unit.
    msg->set_max_steps(input->simulation_options().steps_per_work_unit_part());

    // Add the model.
    if (input->has_reaction_model()) msg->mutable_reaction_model()->CopyFrom(input->reaction_model());
    if (input->has_diffusion_model()) msg->mutable_diffusion_model()->CopyFrom(input->diffusion_model());
    if (input->has_order_parameters()) msg->mutable_order_parameters()->CopyFrom(input->order_parameters());
}

void TrajectoryListSupervisor::buildRunWorkUnitParts(lm::message::RunWorkUnit* msg, uint minWorkUnits)
{
    trajectoryList->buildWorkUnitParts(msg->work_unit_id(), msg, std::max(uint64_t(minWorkUnits), input->simulation_options().parts_per_work_unit()));
}


void TrajectoryListSupervisor::receivedFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg)
{
    PROF_BEGIN(PROF_TRAJECTORY_LIST_WORK_UNIT_FINISHED);

    // Update the trajectory list.
    trajectoryList->processWorkUnitFinished(msg.work_unit_id(), msg);

    PROF_END(PROF_TRAJECTORY_LIST_WORK_UNIT_FINISHED);

    // Call the base class method.
    SimulationSupervisor::receivedFinishedWorkUnit(msg);
}

void TrajectoryListSupervisor::finishSimulation()
{
    // Make sure all the trajectories were finished.
    if (!trajectoryList->allFinished()) THROW_EXCEPTION(lm::ConsistencyException, "at end of simulation there were trajectories not finished");

    // Call the base class method.
    SimulationSupervisor::finishSimulation();
}

void TrajectoryListSupervisor::printPerformanceStatistics()
{
    SimulationSupervisor::printPerformanceStatistics();
    trajectoryList->printTrajectoryStatistics();
}



}
}
