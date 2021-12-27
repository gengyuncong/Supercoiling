/*
 * Copyright 2016-2019 Johns Hopkins University
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

#ifndef LM_MICROENV_MICROENVIRONMENTSUPERVISOR_H_
#define LM_MICROENV_MICROENVIRONMENTSUPERVISOR_H_

#include <string>

#include "hrtime.h"
#include "lm/message/FinishedWorkUnit.pb.h"
#include "lm/microenv/METrajectoryList.h"
#include "lm/microenv/PDETrajectoryList.h"
#include "lm/simulation/SimulationSupervisor.h"
#include "robertslab/Types.h"

namespace lm {
namespace microenv {

using std::string;

class MicroenvironmentSupervisor : public lm::simulation::SimulationSupervisor
{
public:
    static bool registered;
    static bool registerClass();
    static void* allocateObject();

public:
    MicroenvironmentSupervisor();
    virtual ~MicroenvironmentSupervisor();
    virtual std::string getClassName();
    virtual void printPerformanceStatistics();

protected:
    virtual void receivedFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg);
    virtual void startSimulation();
    virtual void startReplicate();
    virtual void buildTrajectoryLists();
    virtual void continueReplicate();
    virtual bool assignWork();
    virtual void finishSimulation();

protected:
    virtual void buildRunWorkUnit(lm::message::RunWorkUnit* msg, bool me=true);

protected:
    uint numberReplicates;
    uint currentReplicateIndex;
    uint numberTimesteps;
    uint currentTimestep;
    double tau;
    double maxTime;
    lm::microenv::METrajectoryList* meTrajectoryList;
    lm::microenv::PDETrajectoryList* pdeTrajectoryList;
    double gridSpacing;
    double gridElementVolume;
    uint32_t numberCells;
    ndarray<double>* cellCoordinates;
    ndarray<uint32_t>* cellGridPoints;
    ndarray<double>* cellVolumes;
    ndarray<int32_t>* cellPreviousCounts;
    ndarray<int32_t>* cellCurrentCounts;
    ndarray<int32_t>* cellFlux;

private:
    long long stats_pdeWorkUnitsSteps;
    double stats_pdeWorkUnitsTime;
    long long stats_timesteps;
    hrtime stats_timestepStartTime;
    hrtime stats_timestepTotalTime;
    hrtime stats_timestepPDETime;
    hrtime stats_timestepMETime;
    hrtime stats_timestepReconcileTime;
};

}
}

#endif
