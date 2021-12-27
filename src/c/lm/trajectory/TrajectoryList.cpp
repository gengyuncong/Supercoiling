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

#include <list>
#include <map>
#include <set>

#include "lm/Exceptions.h"
#include "lm/Types.h"
#include "lm/message/RunWorkUnit.pb.h"
#include "lm/trajectory/Trajectory.h"
#include "lm/trajectory/TrajectoryList.h"

using std::list;
using std::map;
using std::set;

namespace lm {
namespace trajectory {

TrajectoryList::TrajectoryList()
{
}

TrajectoryList::~TrajectoryList()
{
    removeAllTrajectories();
}

void TrajectoryList::addTrajectory(uint64_t id, Trajectory* trajectory)
{
    if (id != trajectory->getID()) THROW_EXCEPTION(lm::RuntimeException, "Trajectory id's do not match: %d %d", id, trajectory->getID());
    if (trajectories.count(id) > 0) THROW_EXCEPTION(lm::RuntimeException, "Tried to add a trajectory that already exists: %d", id);

    // Add the trajectory to the main list.
    trajectories[id] = trajectory;

    // Add the trajectory to the correct status list.
    getTrajectorySet(trajectory->getStatus()).insert(id);
}

void TrajectoryList::removeTrajectory(uint64_t id)
{
    if (trajectories.count(id) > 0)
    {
        Trajectory* trajectory =  trajectories[id];

        // Remove the trajectory from the appropriate status list.
        getTrajectorySet(trajectory->getStatus()).erase(id);

        // Remove the trajectory from the main list.
        trajectories.erase(id);

        // Delete the trajectory.
        delete trajectory;
    }
}

void TrajectoryList::removeAllTrajectories()
{
    // Delete all of the trajectories.
    for (std::map<uint64_t,Trajectory*>::iterator it=trajectories.begin(); it!=trajectories.end(); it++)
        delete it->second;

    // Clear all of the lists.
    trajectories.clear();
    abortedTrajectories.clear();
    finishedTrajectories.clear();
    runningTrajectories.clear();
    waitingTrajectories.clear();
    workUnitTrajectories.clear();
}

void TrajectoryList::restartFinishedTrajectories()
{
    for (std::set<uint64_t>::iterator it=finishedTrajectories.begin(); it!=finishedTrajectories.end(); it++)
    {
        trajectories.at(*it)->setStatus(Trajectory::WAITING);
        trajectories.at(*it)->setLimitReached(lm::types::TrajectoryLimit());
        waitingTrajectories.insert(*it);
    }
    finishedTrajectories.clear();
}

uint64_t TrajectoryList::buildWorkUnitParts(uint64_t workUnitId, lm::message::RunWorkUnit* msg, uint64_t numberParts)
{
    set<uint64_t> trajectoriesAdded;
    for (uint64_t i=0; i<numberParts; i++)
    {
        // See if any trajectories are waiting.
        if (waitingTrajectories.size() > 0)
        {
            // Get the next trajectory.
            uint64_t id = pickNextTrajectoryToRun(waitingTrajectories);
            Trajectory* t = trajectories[id];

            // Validate that it really needs to be run.
            if (t->getStatus() != Trajectory::NOT_STARTED && t->getStatus() != Trajectory::WAITING) THROW_EXCEPTION(lm::ConsistencyException, "invalid trajectory picked from the waiting list: id %d, status: %d", id, t->getStatus());

            // Track that we added it to the list for this work unit.
            trajectoriesAdded.insert(id);

            // Change the trajectory status and move it to the running list.
            setTrajectoryStatus(id, Trajectory::RUNNING);

            // Add it to the message.
            lm::message::WorkUnit* wu = msg->add_part();
            wu->mutable_initial_state()->CopyFrom(t->getState());
        }
    }

    // Add these trajectories to the work units running map.
    workUnitTrajectories[workUnitId] = trajectoriesAdded;

    return uint64_t(trajectoriesAdded.size());
}

void TrajectoryList::processWorkUnitFinished(uint64_t workUnitId, const lm::message::FinishedWorkUnit& fwuMsg)
{
    //fwuMsg.PrintDebugString(); //DEBUGTODO

    if (workUnitTrajectories.count(workUnitId) == 0) THROW_EXCEPTION(ConsistencyException, "finished work unit not found in the list of running work units: %d", workUnitId);

    // Get the list of trajectories associated with this work unit.
    set<uint64_t> involvedTrajectories = workUnitTrajectories[workUnitId];

    // Remove the work unit from the map.
    workUnitTrajectories.erase(workUnitId);

    //Make sure the sizes between the list and the message are consistent.
    if (uint64_t(involvedTrajectories.size()) != uint64_t(fwuMsg.part_status().size())) THROW_EXCEPTION(ConsistencyException, "number of involved trajectories differed from work units finished message: %d, %d, %d", workUnitId, involvedTrajectories.size(), fwuMsg.part_status_size());

    // Loop over the trajectories in the message.
    for (int i=0; i<fwuMsg.part_status().size(); i++)
    {
        // Get the trajectory id.
        uint64_t trajId = fwuMsg.part_status(i).final_state().trajectory_id();

        //Make sure this trajectory is in the list.
        if (involvedTrajectories.count(trajId) == 0) THROW_EXCEPTION(ConsistencyException, "a trajectory in the work units finished message was not in the involved list: %d, %d", workUnitId, trajId);

        // Remove from the involved list.
        involvedTrajectories.erase(trajId);

        // Process the work unit part.
        processWorkUnitPartFinished(fwuMsg.part_status(i), trajId);
    }

    // Make sure we proccessed all of the trajectories.
    if (involvedTrajectories.size() != 0) THROW_EXCEPTION(ConsistencyException, "not all of the trajectories in the involved list were processed: %d remaining", involvedTrajectories.size());
}

void TrajectoryList::processWorkUnitPartFinished(const lm::message::WorkUnitStatus& wusBuf, uint64_t trajId)
{
    // Get the trajectory.
    Trajectory* traj = trajectories[trajId];

    // Make sure the state is consistent.
    if (traj->getStatus() != Trajectory::RUNNING) THROW_EXCEPTION(lm::ConsistencyException, "trajectory was not running: %d", trajId);
    if (runningTrajectories.count(trajId) == 0) THROW_EXCEPTION(lm::ConsistencyException, "trajectory was not in the running list: %d", trajId);

    // Update the state of the trajectory.
    traj->setState(wusBuf.final_state());
    traj->incrementWorkUnitsPerformed();

    // See why the trajectory stopped.
    if (wusBuf.status() == lm::message::WorkUnitStatus::STEPS_FINISHED)
    {
        // Update the status of the trajectory and move it to the waiting set.
        setTrajectoryStatus(trajId, Trajectory::WAITING);
    }
    else if (wusBuf.status() == lm::message::WorkUnitStatus::LIMIT_REACHED)
    {
        traj->setLimitReached(wusBuf.limit_reached());

        // Update the status of the trajectory and move it to the finished set.
        setTrajectoryStatus(trajId, Trajectory::FINISHED);
    }
    else if (wusBuf.status() == lm::message::WorkUnitStatus::ERROR)
    {
        if (wusBuf.has_error_message())
        {
            THROW_EXCEPTION(lm::RuntimeException, "error received in work unit status: %s", wusBuf.error_message().c_str());
        }
        else
        {
            THROW_EXCEPTION(lm::RuntimeException, "error received in work unit status: no message specified");
        }
    }
    else
    {
        THROW_EXCEPTION(lm::RuntimeException, "unknown work unit status: %d", wusBuf.status());
    }
}

std::set<uint64_t>& TrajectoryList::getTrajectorySet(Trajectory::Status status)
{
    // Return the correct set according to the status.
    switch (status)
    {
    case Trajectory::ABORTED: return abortedTrajectories;
    case Trajectory::FINISHED: return finishedTrajectories;
    case Trajectory::NOT_STARTED: return waitingTrajectories;
    case Trajectory::RUNNING: return runningTrajectories;
    case Trajectory::WAITING: return waitingTrajectories;
    }
    THROW_EXCEPTION(lm::RuntimeException, "unknown trajectory status: %d", status);
}

uint64_t TrajectoryList::pickNextTrajectoryToRun(const std::set<uint64_t>& possibleTrajectories) const
{
    return *(possibleTrajectories.begin());
}

void TrajectoryList::setTrajectoryStatus(uint64_t id, Trajectory::Status newStatus)
{
    // Get the current status.
    Trajectory::Status currentStatus = trajectories[id]->getStatus();

    // Remove it from the current list.
    getTrajectorySet(currentStatus).erase(id);

    // Add it to the new list.
    getTrajectorySet(newStatus).insert(id);

    // Update the trajectory's status.
    trajectories[id]->setStatus(newStatus);
}

void TrajectoryList::printTrajectoryStatistics() const
{
}

}
}
