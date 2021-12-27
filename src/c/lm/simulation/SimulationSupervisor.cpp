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

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "hrtime.h"
#include "lm/Exceptions.h"
#include "lm/Print.h"
#include "lm/input/Input.pb.h"
#include "lm/message/Communicator.h"
#include "lm/message/FinishedCheckpointing.pb.h"
#include "lm/message/FinishedWorkUnit.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/message/PerformCheckpointing.pb.h"
#include "lm/message/RunWorkUnit.pb.h"
#include "lm/message/StartedWorkUnit.pb.h"
#include "lm/message/WorkUnit.pb.h"
#include "lm/simulation/SimulationSupervisor.h"
#include "lm/slot/SlotList.h"
#include "lm/thread/Thread.h"
#include "lm/thread/Worker.h"
#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"

using lm::message::Communicator;
using lm::message::Endpoint;

namespace lm {
namespace simulation {

SimulationSupervisor::SimulationSupervisor()
:communicator(NULL),input(NULL),slotList(NULL),simulationRunning(false),performingCheckpoint(false),nextWorkUnitID(1),
stats_workUnits(0),stats_workUnitsParts(0),stats_minWorkUnitId(0),stats_maxWorkUnitId(0),stats_workUnitsSteps(0),stats_workUnitTime(0.0)
{
    // Create the communicator.
    communicator = lm::message::Communicator::createObjectOfDefaultSubclass(false);
}

SimulationSupervisor::~SimulationSupervisor()
{
    if (communicator != NULL) delete communicator; communicator = NULL;
}

void SimulationSupervisor::setId(int64_t id)
{
    this->id = id;
}

void SimulationSupervisor::setInput(lm::input::Input* input)
{
    this->input = input;
}

void SimulationSupervisor::setSlotList(lm::slot::SlotList* slotList)
{
    this->slotList = slotList;
}

void SimulationSupervisor::setOutputWriterAddress(Endpoint outputWriterAddress)
{
    this->outputWriterAddress = outputWriterAddress;
}

void SimulationSupervisor::wake()
{
    lm::message::Message msg;
    msg.mutable_ping_target()->set_id(0);
    communicator->sendMessage(communicator->getSourceAddress(), &msg);
}

int SimulationSupervisor::run()
{
    try
    {
        Print::printf(Print::INFO, "Supervisor %s of type %s started.", Communicator::printableAddress(communicator->getSourceAddress()).c_str(), getClassName().c_str());

        // Start the simulation.
        startSimulation();

        // Assign the initial work.
        if (!assignWork())
        {
            // If assignWork returned false, there was nothing to be done.
            Print::printf(Print::INFO, "No work to be performed.");
            finishSimulation();
        }

        // Loop reading messages.
        lm::message::Message message;
        while (running && simulationRunning)
        {
            // Read the next message.
            communicator->receiveMessage(&message);

            //message.PrintDebugString();

            // Do something with the message.
            if (message.has_started_work_unit())
            {
                receivedStartedWorkUnit(message.started_work_unit());
            }
            else if (message.has_finished_work_unit())
            {
                receivedFinishedWorkUnit(message.finished_work_unit());
            }
            else if (message.has_perform_checkpointing())
            {
                receivedPerformCheckpointing(message.perform_checkpointing());
            }
            else if (message.has_finished_checkpointing())
            {
                receivedFinishedCheckpointing(message.finished_checkpointing());
            }
            else if (message.has_ping_target())
            {
            }
            // Hook for adding generic behavior to this loop in children.
            else if (processEvent(message))
            {
            }
            else
            {
                Print::printf(Print::FATAL, "Supervisor received an unknown message: {\n%s}",message.DebugString().c_str());
            }

            // Clear the message object so it can be used again.
            message.Clear();
        }

        Print::printf(Print::INFO, "Supervisor %s finished.", Communicator::printableAddress(communicator->getSourceAddress()).c_str());

        return 0;
    }
    catch (lm::Exception& e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (lm::Exception* e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e->what(), __FILE__, __LINE__);
    }
    catch (std::exception& e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (...)
    {
        Print::printf(Print::FATAL, "Unknown Exception during execution (%s:%d)", __FILE__, __LINE__);
    }
    exit(-1);
}

bool SimulationSupervisor::processEvent(lm::message::Message& msg)
{
    return false;
}

void SimulationSupervisor::receivedStartedWorkUnit(const lm::message::StartedWorkUnit& msg)
{
    Print::printf(Print::VERBOSE_DEBUG, "Work unit %d started.",msg.work_unit_id());
}

void SimulationSupervisor::receivedFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg)
{
    Print::printf(Print::VERBOSE_DEBUG, "Work unit %d finished in %0.3f s.",msg.work_unit_id(),msg.run_time());

    // Collect global performance stats for printPerformanceStatistics
    stats_workUnits++;
    stats_minWorkUnitId = std::min(stats_minWorkUnitId,(long long)msg.work_unit_id());
    stats_maxWorkUnitId = std::max(stats_maxWorkUnitId,(long long)msg.work_unit_id());
    stats_workUnitsSteps += msg.steps();
    stats_workUnitTime += msg.run_time();
    for (int i=0; i<msg.part_status_size(); i++)
        stats_workUnitsParts++;

    PROF_BEGIN(PROF_SLOTS_WORK_UNIT_FINISHED);

    // Update the slots list.
    slotList->workUnitFinished(msg);

    PROF_END(PROF_SLOTS_WORK_UNIT_FINISHED);

    // If we are not performing a checkpoint, distribute more work.
    if (!performingCheckpoint)
    {
        // Try to assign more work.
        if (!assignWork())
        {
            finishSimulation();
        }
    }
}

void SimulationSupervisor::startSimulation()
{
    simulationStartTime = getHrTime();
    simulationRunning = true;

    // Send a message to the hypervisor that we are started.
    lm::message::Message msg;
    msg.mutable_supervisor_started()->mutable_address()->CopyFrom(communicator->getSourceAddress());
    msg.mutable_supervisor_started()->set_id(id);
    communicator->sendMessage(communicator->getHypervisorAddress(), &msg);
}


void SimulationSupervisor::buildRunWorkUnitHeader(lm::message::RunWorkUnit* msg)
{
    // Set the work unit id.
    msg->set_work_unit_id(nextWorkUnitID++);

    // Set the writer address.
    msg->mutable_output_address()->CopyFrom(outputWriterAddress);

    // Set the work unit's effective supervisor address.
    msg->mutable_supervisor_address()->CopyFrom(communicator->getSourceAddress());
}

void SimulationSupervisor::finishSimulation()
{
    // Mark that the simulation is finished so we exit our message loop.
    simulationRunning = false;

    // Send a message to the hypervisor that we are finished.
    lm::message::Message msg;
    msg.mutable_supervisor_finished()->mutable_address()->CopyFrom(communicator->getSourceAddress());
    msg.mutable_supervisor_finished()->set_id(id);
    communicator->sendMessage(communicator->getHypervisorAddress(), &msg);
}

void SimulationSupervisor::receivedPerformCheckpointing(const lm::message::PerformCheckpointing& msg)
{
    if (simulationRunning)
    {
        // Mark that we need to perform a checkpoint, so distribution of work units should pause until checkpointing is finished.
        performingCheckpoint = true;
    }
}

void SimulationSupervisor::receivedFinishedCheckpointing(const lm::message::FinishedCheckpointing& msg)
{
    // Mark that we are done checkpointing.
    performingCheckpoint = false;

    // Resume distribution of work.
    if (!assignWork())
    {
        finishSimulation();
    }
}

void SimulationSupervisor::printPerformanceStatistics()
{
    if (stats_workUnits > 0)
        Print::printf(Print::INFO, "Supervisor %s finished %lld work units (ids in range %lld to %lld) with %lld parts in the last interval. %lld steps in %0.3e seconds (%0.3e steps/second).",Communicator::printableAddress(communicator->getSourceAddress()).c_str(),stats_workUnits,stats_minWorkUnitId,stats_maxWorkUnitId,stats_workUnitsParts, stats_workUnitsSteps, stats_workUnitTime, double(stats_workUnitsSteps)/stats_workUnitTime);
    else
        Print::printf(Print::INFO, "Supervisor %s finished %lld work units in the last interval.",Communicator::printableAddress(communicator->getSourceAddress()).c_str(),stats_workUnits,stats_minWorkUnitId,stats_maxWorkUnitId,stats_workUnitsParts, stats_workUnitsSteps, stats_workUnitTime, double(stats_workUnitsSteps)/stats_workUnitTime);
    stats_workUnits = 0;
    stats_workUnitsParts = 0;
    stats_minWorkUnitId = std::numeric_limits<long long>::max();
    stats_maxWorkUnitId = 0;
    stats_workUnitsSteps = 0;
    stats_workUnitTime = 0.0;
}

}
}
