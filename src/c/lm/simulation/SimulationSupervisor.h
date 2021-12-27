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

#ifndef LM_SIMULATION_SUPERVISOR_H
#define LM_SIMULATION_SUPERVISOR_H

#include <string>
#include <google/protobuf/message.h>

#include "hrtime.h"
#include "lm/Exceptions.h"
#include "lm/input/Input.pb.h"
#include "lm/message/Communicator.h"
#include "lm/message/FinishedCheckpointing.pb.h"
#include "lm/message/FinishedWorkUnit.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/message/PerformCheckpointing.pb.h"
#include "lm/message/RunWorkUnit.pb.h"
#include "lm/message/StartedWorkUnit.pb.h"
#include "lm/message/WorkUnit.pb.h"
#include "lm/slot/SlotList.h"
#include "lm/thread/Thread.h"
#include "lm/thread/Worker.h"

namespace lm {
namespace simulation {

class SimulationSupervisor : public lm::thread::Worker
{
public:
    SimulationSupervisor();
    virtual ~SimulationSupervisor();
    virtual std::string getClassName()=0;

    virtual void setId(int64_t id);
    virtual void setInput(lm::input::Input* input);
    virtual void setSlotList(lm::slot::SlotList* slotList);
    virtual void setOutputWriterAddress(Endpoint outputWriterAddress);

    void wake();
    virtual void printPerformanceStatistics();

protected:
    virtual int run();

    virtual void receivedStartedWorkUnit(const lm::message::StartedWorkUnit& msg);
    virtual void receivedFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg);
    virtual void receivedPerformCheckpointing(const lm::message::PerformCheckpointing& msg);
    virtual void receivedFinishedCheckpointing(const lm::message::FinishedCheckpointing& msg);

    virtual void startSimulation();
    virtual bool assignWork()=0;
    virtual void finishSimulation();
    virtual bool processEvent(lm::message::Message& msg);

    virtual void buildRunWorkUnitHeader(lm::message::RunWorkUnit* msg);

protected:
    lm::message::Communicator* communicator;

    lm::input::Input* input;
    lm::slot::SlotList* slotList;
    Endpoint outputWriterAddress;
    int64_t id;

    bool simulationRunning;
    bool performingCheckpoint;

    long long nextWorkUnitID;

protected:
    hrtime simulationStartTime;

    long long stats_workUnits;
    long long stats_workUnitsParts;
    long long stats_minWorkUnitId;
    long long stats_maxWorkUnitId;
    long long stats_workUnitsSteps;
    double stats_workUnitTime;
};

}
}

#endif // LM_MAIN_SUPERVISOR_H
