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
#ifndef LM_MAIN_HYPERVISOR_H
#define LM_MAIN_HYPERVISOR_H

#include <map>
#include <string>
#include <vector>

//#include <google/protobuf/message.h>

//#include "hrtime.h"
//#include "lm/Exceptions.h"
#include "lm/input/Input.pb.h"
//#include "lm/types/BoundaryConditions.pb.h"
//#include "lm/input/DiffusionModel.pb.h"
//#include "lm/io/OrderParameters.pb.h"
//#include "lm/input/ReactionModel.pb.h"
//#include "lm/io/Tilings.pb.h"
#include "lm/message/Communicator.h"
#include "lm/message/FinishedCheckpointing.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/message/ResourcesAvailable.pb.h"
#include "lm/message/StartWorkUnitRunner.pb.h"
#include "lm/message/StartedCheckpointSignaler.pb.h"
#include "lm/message/StartedOutputWriter.pb.h"
#include "lm/message/StartedWorkUnitRunner.pb.h"
#include "lm/message/SupervisorFinished.pb.h"
#include "lm/message/SupervisorStarted.pb.h"
//#include "lm/message/WorkUnit.pb.h"
//#include "lm/oparam/OParams.h"
#include "lm/resource/ResourceMap.h"
#include "lm/slot/SlotList.h"
//#include "lm/trajectory/TrajectoryList.h"
#include "lm/thread/Thread.h"
#include "lm/thread/Worker.h"
//#include "lm/tiling/Tilings.h"

using std::map;
using std::string;
using std::vector;

namespace lm {
namespace simulation {

class JobHypervisor : public lm::thread::Worker
{
public:
    JobHypervisor();
    virtual ~JobHypervisor();
    virtual string getClassName()=0;

    virtual void init();
    void setOutputWriterClassName(string outputWriterClassName) {this->outputWriterClassName = outputWriterClassName;}
    void setSupervisorClassName(string supervisorClassName) {this->supervisorClassName = supervisorClassName;}
    void setResourceMap(lm::resource::ResourceMap& resourceMap) {this->resourceMap = resourceMap;}
    void setSimulationFilename(vector<string> simulationInputFilenames, string simulationOutputFilename) {this->simulationInputFilenames = simulationInputFilenames; this->simulationOutputFilename = simulationOutputFilename;}
    void setSolverClassName(string solverClassName) {this->solverClassName = solverClassName;}
    void setUseCPUAffinity(bool useCPUAffinity) {this->useCPUAffinity = useCPUAffinity;}
    void wake();
    virtual void printPerformanceStatistics()=0;

protected:
    virtual int run();

    virtual void receivedResourceAvailable(const lm::message::ResourcesAvailable& msg);
    virtual void allResourcesRegistered();
    virtual void startOutputWriter();
    virtual void startCheckpointSignaler();
    virtual void startOutputPerformanceSignaler();
    virtual void startWorkUnitRunners();

    virtual void receivedStartedOutputWriter(const lm::message::StartedOutputWriter& msg);
    virtual void receivedStartedCheckpointSignaler(const lm::message::StartedCheckpointSignaler& msg);
    virtual void receivedStartedOutputPerformanceSignaler(const lm::message::StartedOutputPerformanceSignaler& msg);
    virtual void receivedStartedWorkUnitRunner(const lm::message::StartedWorkUnitRunner & msg);
    virtual void receivedPerformCheckpointing(const lm::message::PerformCheckpointing& msg);
    virtual void receivedPerformOutputPerformance(const lm::message::PerformOutputPerformance& msg);
    virtual void receivedFinishedCheckpointing(const lm::message::FinishedCheckpointing& msg);
    virtual void supervisorStartedEvent(const lm::message::SupervisorStarted& msg);
    virtual void supervisorFinishedEvent(const lm::message::SupervisorFinished& msg);
    virtual bool processEvent(lm::message::Message& msg);

protected:
    virtual void startJobIfAllWorkersStarted();
    virtual void startJob();
    virtual bool startRemainingSimulations()=0;
    virtual void finishJob();

protected:
    lm::message::Communicator* communicator;
    bool hasCheckpointSignalerStarted;
    bool hasOutputPerformanceSignalerStarted;
    bool hasOutputWriterStarted;
    bool haveAllWorkUnitRunnersStarted;

    lm::input::Input input;
    std::string outputWriterClassName;
    Endpoint outputWriterAddress;
    bool performingCheckpoint;
    lm::resource::ResourceMap resourceMap;
    vector<string> simulationInputFilenames;
    string simulationOutputFilename;
    uint64_t simulationPhase;
    bool jobRunning;
    lm::slot::SlotList slots;
    std::string supervisorClassName;
    std::string solverClassName;
    bool useCPUAffinity;
    long long workUnitCount;
};

}
}

#endif // LM_MAIN_HYPERVISOR_H
