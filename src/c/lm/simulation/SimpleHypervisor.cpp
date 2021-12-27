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


#include <map>
#include <string>

#include "hrtime.h"
#include "lm/ClassFactory.h"
#include "lm/Print.h"
#include "lm/main/Globals.h"
#include "lm/message/SupervisorFinished.pb.h"
#include "lm/message/SupervisorStarted.pb.h"
#include "lm/simulation/SimulationSupervisor.h"
#include "lm/simulation/SimpleHypervisor.h"

using std::map;
using std::string;
using lm::resource::ResourceMap;

namespace lm {
namespace simulation {

bool SimpleHypervisor::registered=SimpleHypervisor::registerClass();

bool SimpleHypervisor::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::simulation::JobHypervisor","lm::simulation::SimpleHypervisor",&SimpleHypervisor::allocateObject);
    return true;
}

void* SimpleHypervisor::allocateObject()
{
    return new SimpleHypervisor();
}

SimpleHypervisor::SimpleHypervisor()
:simulationSupervisor(NULL),supervisorStarted(false),supervisorFinished(false)
{
}

SimpleHypervisor::~SimpleHypervisor()
{
}

std::string SimpleHypervisor::getClassName()
{
    return "lm::simulation::SimpleHypervisor";
}

bool SimpleHypervisor::startRemainingSimulations()
{
    if (simulationSupervisor == NULL)
    {
        // Create the supervisor.
        simulationSupervisor = static_cast<lm::simulation::SimulationSupervisor*>(lm::ClassFactory::getInstance().allocateObjectOfClass("lm::simulation::SimulationSupervisor",supervisorClassName));
        simulationSupervisor->setId(0);
        simulationSupervisor->setInput(&input);
        simulationSupervisor->setSlotList(&slots);
        simulationSupervisor->setOutputWriterAddress(outputWriterAddress);

        // Start the supervisor.
        simulationSupervisor->start();
    }
    if (!supervisorStarted) return true;

    //return true;
    return !supervisorFinished;
}

void SimpleHypervisor::supervisorStartedEvent(const lm::message::SupervisorStarted& msg)
{
    // Mark that the supervisor has finished.
    supervisorStarted = true;

    // Call the base method now that we have processed the message.
    JobHypervisor::supervisorStartedEvent(msg);
}

void SimpleHypervisor::supervisorFinishedEvent(const lm::message::SupervisorFinished& msg)
{
    // Mark that the supervisor has finished.
    supervisorFinished = true;

    // Call the base method now that we have processed the message.
    JobHypervisor::supervisorFinishedEvent(msg);
}

void SimpleHypervisor::printPerformanceStatistics()
{
    if (simulationSupervisor != NULL)
    {
        simulationSupervisor->printPerformanceStatistics();
    }
}


}
}
