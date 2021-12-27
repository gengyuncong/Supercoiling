/*
 * Copyright 2019 Johns Hopkins University
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

#ifndef LM_PARAMSWEEP_PARAMETERSWEEPHYPERVISOR_H
#define LM_PARAMSWEEP_PARAMETERSWEEPHYPERVISOR_H

#include <map>
#include <string>
#include <vector>

#include "hrtime.h"
#include "lm/message/SupervisorFinished.pb.h"
#include "lm/message/SupervisorStarted.pb.h"
#include "lm/simulation/JobHypervisor.h"
#include "lm/simulation/SimulationSupervisor.h"
#include "lm/slot/SlotList.h"

using std::map;
using std::string;
using std::vector;

namespace lm {
namespace paramsweep {

class ParameterSimulation;
class SlotGroup;

class ParameterSweepHypervisor : public lm::simulation::JobHypervisor
{
public:
    static bool registered;
    static bool registerClass();
    static void* allocateObject();

public:
    ParameterSweepHypervisor();
    virtual ~ParameterSweepHypervisor();
    virtual std::string getClassName();
    virtual void init();
    virtual void printPerformanceStatistics();

protected:
    virtual void startJob();
    virtual bool startRemainingSimulations();
    virtual void supervisorStartedEvent(const lm::message::SupervisorStarted& msg);
    virtual void supervisorFinishedEvent(const lm::message::SupervisorFinished& msg);

protected:
    vector<SlotGroup> slotGroups;
    vector<string> parameterKeys;
    vector<ParameterSimulation*> parameterRuns;
    size_t parameterSimulationsRunning;
    size_t parameterSimulationsFinished;
};

class SlotGroup
{
public:
    enum Status {AVAILABLE, IN_USE};
    static const char* status_strings[];
    SlotGroup(lm::slot::SlotList slots);
    ~SlotGroup();
    Status state;
    lm::slot::SlotList slots;
    size_t parameterSimulationIndex;
};

class ParameterSimulation
{
public:
    enum Status {NOT_STARTED, STARTED, RUNNING, FINISHED};
    static const char* status_strings[];

public:
    ParameterSimulation(int index, const lm::input::Input& originalInput, const map<string,string>& parameterSet);
    ~ParameterSimulation();
    void releaseResources();
    Status state;
    lm::input::Input* input;
    lm::simulation::SimulationSupervisor* supervisor;
    hrtime startTime;
    hrtime stopTime;
    string name;
    lm::message::Endpoint supervisorEndpoint;
};

}
}

#endif /* LM_PARAMSWEEP_PARAMETERSWEEPHYPERVISOR_H */
