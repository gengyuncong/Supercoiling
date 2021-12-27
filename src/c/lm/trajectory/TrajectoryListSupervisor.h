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

#ifndef LM_SIMULATION_TRAJECTORYSIMULATIONSUPERVISOR_H
#define LM_SIMULATION_TRAJECTORYSIMULATIONSUPERVISOR_H

#include "lm/simulation/SimulationSupervisor.h"
#include "lm/trajectory/TrajectoryList.h"
#include "lm/types/TrajectoryBarrier.pb.h"
#include "lm/types/TrajectoryLimits.pb.h"

namespace lm {
namespace trajectory {

class TrajectoryListSupervisor : public lm::simulation::SimulationSupervisor
{
public:
    TrajectoryListSupervisor();
    virtual ~TrajectoryListSupervisor();
    virtual void printPerformanceStatistics();

protected:

    virtual void receivedFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg);
    virtual void startSimulation();
    virtual bool assignWork();
    virtual void buildTrajectoryList()=0;
    virtual void buildTrajectoryOptions()=0;
    virtual void finishSimulation();

    virtual void buildRunWorkUnitHeader(lm::message::RunWorkUnit* msg);
    virtual void buildRunWorkUnitParts(lm::message::RunWorkUnit* msg, uint minWorkUnits);

protected:
    lm::trajectory::TrajectoryList* trajectoryList;
    lm::types::TrajectoryLimits trajectoryLimits;
    lm::types::TrajectoryBarriers trajectoryBarriers;
};

}
}

#endif
