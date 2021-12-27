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

#ifndef LM_SIMULATION_SIMPLEHYPERVISOR_H
#define LM_SIMULATION_SIMPLEHYPERVISOR_H

#include <string>

#include "lm/message/SupervisorFinished.pb.h"
#include "lm/message/SupervisorStarted.pb.h"
#include "lm/simulation/JobHypervisor.h"
#include "lm/simulation/SimulationSupervisor.h"

namespace lm {
namespace simulation {

class SimpleHypervisor : public lm::simulation::JobHypervisor
{
public:
    static bool registered;
    static bool registerClass();
    static void* allocateObject();

public:
    SimpleHypervisor();
    virtual ~SimpleHypervisor();
    virtual std::string getClassName();
    virtual bool startRemainingSimulations();
    virtual void printPerformanceStatistics();

    virtual void supervisorStartedEvent(const lm::message::SupervisorStarted& msg);
    virtual void supervisorFinishedEvent(const lm::message::SupervisorFinished& msg);

protected:
    lm::simulation::SimulationSupervisor* simulationSupervisor;
    bool supervisorStarted;
    bool supervisorFinished;
};

}
}

#endif /* LM_SIMULATION_SIMPLEHYPERVISOR_H_ */
