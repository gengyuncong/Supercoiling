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


#include <map>
#include <string>
#include <vector>

#include "hrtime.h"
#include "lm/ClassFactory.h"
#include "lm/Print.h"
#include "lm/io/props/PropertiesInputReader.h"
#include "lm/io/props/TabularPropertiesFile.h"
#include "lm/main/Globals.h"
#include "lm/message/SupervisorFinished.pb.h"
#include "lm/message/SupervisorStarted.pb.h"
#include "lm/paramsweep/ParameterSweepHypervisor.h"
#include "lm/simulation/SimulationSupervisor.h"

using std::map;
using std::string;
using std::vector;
using lm::resource::ResourceMap;

namespace lm {
namespace paramsweep {

bool ParameterSweepHypervisor::registered=ParameterSweepHypervisor::registerClass();

bool ParameterSweepHypervisor::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::simulation::JobHypervisor","lm::paramsweep::ParameterSweepHypervisor",&ParameterSweepHypervisor::allocateObject);
    return true;
}

void* ParameterSweepHypervisor::allocateObject()
{
    return new ParameterSweepHypervisor();
}

ParameterSweepHypervisor::ParameterSweepHypervisor()
:parameterSimulationsRunning(0),parameterSimulationsFinished(0)
{
}

ParameterSweepHypervisor::~ParameterSweepHypervisor()
{
    // Free the simulation runs.
    for (size_t i=0; i<parameterRuns.size(); i++) delete parameterRuns[i];
    parameterRuns.clear();
}

std::string ParameterSweepHypervisor::getClassName()
{
    return "lm::simulation::ParameterSweepHypervisor";
}

void ParameterSweepHypervisor::init()
{
    // Call the base class method.
    JobHypervisor::init();

    // Open the parameters data file.
    if (parameterSweepFilename == "") THROW_EXCEPTION(lm::RuntimeException, "No parameter-sweep filename was specified.");
    lm::io::props::TabularPropertiesFile f(parameterSweepFilename);
    if (!f.exists() || !f.isFile()) THROW_EXCEPTION(lm::RuntimeException, "Invalid parameter-sweep filename specified.");
    f.openRead();

    // Get the keys.
    parameterKeys = f.getKeys();

    // Read all of the properties.
    while (!f.isEof())
    {
        // Get the next set of properties.
        map<string,string> params = f.readNextProperties();
        if (params.size() != 0)
        {
            // Create an object for the parameter simulation.
            ParameterSimulation* s = new ParameterSimulation(static_cast<int>(parameterRuns.size()), input, params);

            // Add the simulation to the list.
            parameterRuns.push_back(s);
        }
    }

    Print::printf(Print::INFO, "ParameterSweepHypervisor loaded %d parameter sets with %d parameters in each from %s.", parameterRuns.size(), parameterKeys.size(), parameterSweepFilename.c_str());

    // Close the file.
    f.close();

}

void ParameterSweepHypervisor::startJob()
{
    // Split the slots into individual slot groups.
    vector<lm::slot::SlotList> splitSlots = slots.split(static_cast<size_t>(parameterSweepParallel));
    for (size_t i=0; i<splitSlots.size(); i++)
        slotGroups.push_back(SlotGroup(splitSlots[i]));

    // Call the base clas method.
    JobHypervisor::startJob();
}

bool ParameterSweepHypervisor::startRemainingSimulations()
{
    bool anyUnfinished = false;

    // Go through the simulations and see if any need to be started.
    for (size_t i=0; i<parameterRuns.size(); i++)
    {
        // See if we have any that have not been started.
        if (parameterRuns[i]->state == ParameterSimulation::NOT_STARTED)
        {
            // Mark that some are not finished.
            anyUnfinished = true;

            // Go through the slot groups and see if one is available.
            for (size_t slotGroupIndex=0; slotGroupIndex<slotGroups.size(); slotGroupIndex++)
            {
                if (slotGroups[slotGroupIndex].state == SlotGroup::AVAILABLE)
                {
                    // Increment the running simulation count.
                    parameterSimulationsRunning++;

                    // Mark that the slot group is in use.
                    slotGroups[slotGroupIndex].state = SlotGroup::IN_USE;
                    slotGroups[slotGroupIndex].parameterSimulationIndex = i;

                    // Mark that the simulation was started.
                    parameterRuns[i]->state = ParameterSimulation::STARTED;

                    // Start the parameter simulation.
                    Print::printf(Print::INFO, "Parameter sweep simulation supervisor started for parameter set %s.",parameterRuns[i]->name.c_str());
                    parameterRuns[i]->supervisor->setSlotList(&slotGroups[slotGroupIndex].slots);
                    parameterRuns[i]->supervisor->setOutputWriterAddress(outputWriterAddress);
                    parameterRuns[i]->supervisor->start();
                    break;
                }
            }
        }
        else if (parameterRuns[i]->state == ParameterSimulation::STARTED || parameterRuns[i]->state == ParameterSimulation::RUNNING)
        {
            // Mark that some are not finished.
            anyUnfinished = true;
        }
    }

    //return true;
    return anyUnfinished;
}

void ParameterSweepHypervisor::supervisorStartedEvent(const lm::message::SupervisorStarted& msg)
{
    // Mark that the supervisor has started.
    parameterRuns[static_cast<size_t>(msg.id())]->state = ParameterSimulation::RUNNING;
    parameterRuns[static_cast<size_t>(msg.id())]->startTime = getHrTime();
    parameterRuns[static_cast<size_t>(msg.id())]->stopTime = parameterRuns[static_cast<size_t>(msg.id())]->startTime;

    // Call the base method now that we have processed the message.
    JobHypervisor::supervisorStartedEvent(msg);
}

void ParameterSweepHypervisor::supervisorFinishedEvent(const lm::message::SupervisorFinished& msg)
{
    // Mark that the simulation has finished and free any resources.
    parameterRuns[static_cast<size_t>(msg.id())]->state = ParameterSimulation::FINISHED;
    parameterRuns[static_cast<size_t>(msg.id())]->stopTime = getHrTime();
    parameterRuns[static_cast<size_t>(msg.id())]->releaseResources();

    // Mark that the slot group is available.
    bool foundSlotGroup = false;
    for (size_t i=0; i<slotGroups.size(); i++)
    {
        if (slotGroups[i].state == SlotGroup::IN_USE && slotGroups[i].parameterSimulationIndex == static_cast<size_t>(msg.id()))
        {
            slotGroups[i].state = SlotGroup::AVAILABLE;
            slotGroups[i].parameterSimulationIndex = std::numeric_limits<size_t>::max();
            foundSlotGroup = true;
            break;
        }
    }
    if (!foundSlotGroup) THROW_EXCEPTION(lm::RuntimeException, "Could not find the slot group corresponding to parameter set: %d.", msg.id());

    // Decrement the simulation running count.
    parameterSimulationsRunning--;
    parameterSimulationsFinished++;

    // Call the base method now that we have processed the message.
    JobHypervisor::supervisorFinishedEvent(msg);
}

void ParameterSweepHypervisor::printPerformanceStatistics()
{
    Print::printf(Print::INFO, "Parameter sweep status: %d finished, %d running, %d waiting", parameterSimulationsFinished, parameterSimulationsRunning, (static_cast<int>(parameterRuns.size()-parameterSimulationsFinished-parameterSimulationsRunning)));
    Print::printf(Print::INFO, "%10s %20s %07s %11s %9s", "Id", "Name", "Address", "Status", "Time");
    Print::printf(Print::INFO, "-------------------------------------------");
    for (size_t i=0; i<parameterRuns.size(); i++)
    {
        Print::printf(Print::INFO, "%10d %20s %07s %11s %9.3e", i, parameterRuns[i]->name.c_str(), lm::message::Communicator::printableAddress(parameterRuns[i]->supervisorEndpoint).c_str(), ParameterSimulation::status_strings[parameterRuns[i]->state], convertHrToSeconds(parameterRuns[i]->stopTime-parameterRuns[i]->startTime));
    }

    for (size_t i=0; i<parameterRuns.size(); i++)
    {
        if (parameterRuns[i]->state == ParameterSimulation::STARTED || parameterRuns[i]->state == ParameterSimulation::RUNNING)
        {
            parameterRuns[i]->supervisor->printPerformanceStatistics();
        }
    }
}

const char* SlotGroup::status_strings[] = {"AVAILABLE", "IN_USE"};

SlotGroup::SlotGroup(lm::slot::SlotList slots)
:state(AVAILABLE),slots(slots),parameterSimulationIndex(std::numeric_limits<size_t>::max())
{
}

SlotGroup::~SlotGroup()
{
}

const char* ParameterSimulation::status_strings[] = {"NOT_STARTED", "STARTED", "RUNNING", "FINISHED"};

ParameterSimulation::ParameterSimulation(int index, const lm::input::Input& originalInput, const map<string,string>& parameterSet)
:state(NOT_STARTED),input(NULL),supervisor(NULL),startTime(0),stopTime(0),name("")
{
    // Make a copy of the input.
    input = new lm::input::Input(originalInput);

    // Update the input with the parameters.
    lm::io::props::PropertiesInputReader r;
    r.setPropertiesInto(parameterSet, input, true);

    // Create the supervisor.
    supervisor = static_cast<lm::simulation::SimulationSupervisor*>(lm::ClassFactory::getInstance().allocateObjectOfClass("lm::simulation::SimulationSupervisor", supervisorClassName));
    supervisor->setId(index);
    supervisor->setInput(input);

    // Get the name.
    if (parameterSet.count("name") != 0)
    {
        name = parameterSet.at("name");
    }
    else
    {
        char buffer[24];
        snprintf(buffer, 24, "%d", index);
        name = string(buffer);
    }
    input->mutable_output_options()->set_record_name_prefix("/Sweep/"+name);

    Print::printf(Print::DEBUG, "Parameter sweep simulation supervisor created for parameter set %s.",name.c_str());
    input->PrintDebugString(); //TODO: DEBUG
}

ParameterSimulation::~ParameterSimulation()
{
    releaseResources();
}

void ParameterSimulation::releaseResources()
{
    if (supervisor != NULL) delete supervisor; supervisor = NULL;
    if (input != NULL) delete input; input = NULL;
}

}
}
