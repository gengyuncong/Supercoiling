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

#ifndef LM_TRAJECTORY_TRAJECTORYLIST_H_
#define LM_TRAJECTORY_TRAJECTORYLIST_H_

#include <map>
#include <set>

#include "lm/Types.h"
#include "lm/message/RunWorkUnit.pb.h"
#include "lm/message/FinishedWorkUnit.pb.h"
#include "lm/message/WorkUnitStatus.pb.h"
#include "lm/trajectory/Trajectory.h"

namespace lm {
namespace trajectory {

class TrajectoryList
{
public:
    TrajectoryList();
    virtual ~TrajectoryList();

    // mutators
    virtual void addTrajectory(uint64_t trajectoryID, Trajectory* trajectory);
    virtual void removeTrajectory(uint64_t trajectoryID);
    virtual void removeAllTrajectories();

    // accessors
    virtual bool allFinished() const {return !(anyAborted() || anyRunning() || anyWaiting());}
    virtual bool anyAborted() const {return not abortedTrajectories.empty();}
    virtual bool anyRunning() const {return not runningTrajectories.empty();}
    virtual bool anyWaiting() const {return not waitingTrajectories.empty();}
    virtual size_t numberWaiting() const {return waitingTrajectories.size();}
    virtual uint64_t size() const {return uint64_t(trajectories.size());}
    virtual bool exists(uint64_t id) const {return trajectories.count(id) != 0;}
    virtual bool isTrajectoryAborted(uint64_t trajID) const {return abortedTrajectories.count(trajID) != 0;}
    virtual bool isTrajectoryFinished(uint64_t trajID) const {return finishedTrajectories.count(trajID) != 0;}
    virtual bool isTrajectoryRunning(uint64_t trajID) const {return runningTrajectories.count(trajID) != 0;}
    virtual bool isTrajectoryWaiting(uint64_t trajID) const {return waitingTrajectories.count(trajID) != 0;}
    virtual bool isTrajectoryAborted(Trajectory* traj) const {return abortedTrajectories.count(traj->getID()) != 0;}
    virtual bool isTrajectoryFinished(Trajectory* traj) const {return finishedTrajectories.count(traj->getID()) != 0;}
    virtual bool isTrajectoryRunning(Trajectory* traj) const {return runningTrajectories.count(traj->getID()) != 0;}
    virtual bool isTrajectoryWaiting(Trajectory* traj) const {return waitingTrajectories.count(traj->getID()) != 0;}
    virtual void restartFinishedTrajectories();

    // work unit functions.
    virtual uint64_t buildWorkUnitParts(uint64_t workUnitId, lm::message::RunWorkUnit* msg, uint64_t numberParts);
    virtual void processWorkUnitFinished(uint64_t workUnitId, const lm::message::FinishedWorkUnit& fwuMsg);
    virtual void processWorkUnitPartFinished(const lm::message::WorkUnitStatus& wusBuf, uint64_t trajId);

    virtual void printTrajectoryStatistics() const;

protected:
    virtual std::set<uint64_t>& getTrajectorySet(Trajectory::Status status);
    virtual uint64_t pickNextTrajectoryToRun(const std::set<uint64_t>& possibleTrajectories) const;
    virtual void setTrajectoryStatus(uint64_t id, Trajectory::Status newStatus);

protected:
    std::map<uint64_t,Trajectory*> trajectories;
    std::set<uint64_t> abortedTrajectories;
    std::set<uint64_t> finishedTrajectories;
    std::set<uint64_t> runningTrajectories;
    std::set<uint64_t> waitingTrajectories;
    std::map<uint64_t,std::set<uint64_t> > workUnitTrajectories;
};

}
}

#endif
