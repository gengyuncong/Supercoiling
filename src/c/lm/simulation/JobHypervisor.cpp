/*
 * University of Illinois Open Source License
 * Copyright 2012-2016 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to
 * do so, subject to the following conditions:
 *
 * - Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimers.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimers in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the names of the Roberts Group, Johns Hopkins University,
 * nor the names of its contributors may be used to endorse or
 * promote products derived from this Software without specific prior written
 * permission.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Elijah Roberts, Max Klein
 */

#include <algorithm>
#include <limits>
#include <string>
#include <vector>

#include "hrtime.h"
#include "lm/ClassFactory.h"
#include "lm/Exceptions.h"
#include "lm/Print.h"
#include "lm/input/DiffusionModel.pb.h"
#include "lm/input/Input.pb.h"
#include "lm/input/InputAggregator.h"
#include "lm/io/hdf5/HDF5.h"
#include "lm/io/hdf5/SimulationFile.h"
#include "lm/main/Globals.h"
#include "lm/simulation/JobHypervisor.h"
#include "lm/message/Communicator.h"
#include "lm/message/Endpoint.pb.h"
#include "lm/message/FinishedCheckpointing.pb.h"
#include "lm/message/FinishedWorkUnit.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/message/ResourcesAvailable.pb.h"
#include "lm/message/RunWorkUnit.pb.h"
#include "lm/message/StartWorkUnitRunner.pb.h"
#include "lm/message/StartedWorkUnit.pb.h"
#include "lm/message/StartedWorkUnitRunner.pb.h"
#include "lm/message/WorkUnit.pb.h"
#include "lm/resource/ComputeResources.h"
#include "lm/resource/ResourceMap.h"
#include "lm/slot/Slot.h"
#include "lm/slot/SlotList.h"

using lm::message::Communicator;
using lm::message::Endpoint;
using lm::resource::ComputeResources;
using lm::resource::ResourceMap;
using std::string;
using std::vector;

namespace lm {
namespace simulation {

JobHypervisor::JobHypervisor()
:communicator(NULL),hasCheckpointSignalerStarted(false),hasOutputPerformanceSignalerStarted(false),hasOutputWriterStarted(false),haveAllWorkUnitRunnersStarted(false),
 outputWriterClassName(""),performingCheckpoint(false),
 simulationOutputFilename(""),simulationPhase(0),jobRunning(true),slots(),
 solverClassName(""),useCPUAffinity(false),workUnitCount(0)
{
    // Create the communicator.
    communicator = lm::message::Communicator::createObjectOfDefaultSubclass(true);
}

JobHypervisor::~JobHypervisor()
{
    if (communicator != NULL) delete communicator; communicator = NULL;
}

void JobHypervisor::init()
{
    // Initialize the input object using the input filenames.
    lm::input::InputAggregator agg;
    agg.readFilesInto(simulationInputFilenames, &input);

    input.PrintDebugString(); //TODO: DEBUG
}

void JobHypervisor::wake()
{
    lm::message::Message msg;
    msg.mutable_ping_target()->set_id(0);
    communicator->sendMessage(communicator->getSourceAddress(), &msg);
}

int JobHypervisor::run()
{
    try
    {
        Print::printf(Print::INFO, "Hypervisor %s of type %s started.", Communicator::printableAddress(communicator->getSourceAddress()).c_str(), getClassName().c_str());

        // Loop reading messages.
        lm::message::Message message;
        while (running && jobRunning)
        {
            // Read the next message.
            communicator->receiveMessage(&message);

            //message.PrintDebugString();

            // Do something with the message.
            if (message.has_resources_available())
            {
                receivedResourceAvailable(message.resources_available());
            }
            else if (message.has_started_output_writer())
            {
                receivedStartedOutputWriter(message.started_output_writer());
            }
            else if (message.has_started_checkpoint_signaler())
            {
                receivedStartedCheckpointSignaler(message.started_checkpoint_signaler());
            }
            else if (message.has_started_output_performance_signaler())
            {
                receivedStartedOutputPerformanceSignaler(message.started_output_performance_signaler());
            }
            else if (message.has_started_work_unit_runner())
            {
                receivedStartedWorkUnitRunner(message.started_work_unit_runner());
            }
            else if (message.has_perform_checkpointing())
            {
                receivedPerformCheckpointing(message.perform_checkpointing());
            }
            else if (message.has_perform_output_performance())
            {
                receivedPerformOutputPerformance(message.perform_output_performance());
            }
            else if (message.has_finished_checkpointing())
            {
                receivedFinishedCheckpointing(message.finished_checkpointing());
            }
            else if (message.has_supervisor_started())
            {
                supervisorStartedEvent(message.supervisor_started());
            }
            else if (message.has_supervisor_finished())
            {
                supervisorFinishedEvent(message.supervisor_finished());
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
                Print::printf(Print::ERROR, "Hypervisor received an unknown message: {\n%s}",message.DebugString().c_str());
            }

            // Clear the message object so it can be used again.
            message.Clear();
        }

        Print::printf(Print::INFO, "Hypervisor %s finished.", Communicator::printableAddress(communicator->getSourceAddress()).c_str());

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

bool JobHypervisor::processEvent(lm::message::Message& msg)
{
    return false;
}

void JobHypervisor::receivedResourceAvailable(const lm::message::ResourcesAvailable& msg)
{
    Print::printf(Print::INFO, "Resource controller %s on %s reported %d cpu core(s) and %d gpu device(s). %d hosts registered, %d hosts remaining.", Communicator::printableAddress(msg.controller_address()).c_str(), msg.hostname().c_str(), msg.cpu_cores_size(), msg.gpu_devices_size(), resourceMap.numberHostsRegistered(), resourceMap.numberHostsAllocated());
    if (resourceMap.registerResources(msg))
    {
        allResourcesRegistered();
    }
}

void JobHypervisor::allResourcesRegistered()
{
    // Start the work unit runners.
    Print::printf(Print::INFO, "All resources registered with hypervisor, starting workers.");

    startOutputWriter();
    startCheckpointSignaler();
    startOutputPerformanceSignaler();
    startWorkUnitRunners();
}

void JobHypervisor::startOutputWriter()
{
    // Reserve a core for the output writer if the option is set.
    if (shouldReserveOutputCore)
    {
        // Start the output writer.
        ComputeResources resources = resourceMap.reserveCPUCores(communicator->getHostname(), 1);
        lm::message::Message msg;
        msg.mutable_start_output_writer()->set_use_cpu_affinity(useCPUAffinity);
        msg.mutable_start_output_writer()->set_cpu(resources.cpuCores[0]);
        msg.mutable_start_output_writer()->set_output_filename(simulationOutputFilename);
        msg.mutable_start_output_writer()->set_output_writer_class(outputWriterClassName);
        communicator->sendMessage(resources.controllerAddress, &msg);
        Print::printf(Print::INFO, "Reserved core %d on %s for the output writer.", resources.cpuCores[0], Communicator::printableAddress(resources.controllerAddress).c_str());
    }
    else
    {
        // Start the output writer.
        ComputeResources resources = resourceMap.reserveCPUCores(communicator->getHostname(), 0);
        lm::message::Message msg;
        msg.mutable_start_output_writer()->set_output_filename(simulationOutputFilename);
        msg.mutable_start_output_writer()->set_output_writer_class(outputWriterClassName);
        // thread 1 should be the resource controller
        communicator->sendMessage(resources.controllerAddress, &msg);
        Print::printf(Print::INFO, "Output writer is sharing resources on %s.", Communicator::printableAddress(resources.controllerAddress).c_str());
    }
}

void JobHypervisor::startCheckpointSignaler()
{
    //See if we need to start a checkpoint signaler.
    if (checkpointInterval > 0)
    {
        // Start the checkpoint signaler.
        ComputeResources resources = resourceMap.reserveCPUCores(communicator->getHostname(), 0);
        lm::message::Message msg;
        msg.mutable_start_checkpoint_signaler()->set_checkpoint_interval((int)checkpointInterval);
        communicator->sendMessage(resources.controllerAddress, &msg);
    }
    else
    {
        hasCheckpointSignalerStarted = true;
    }
}

void JobHypervisor::startOutputPerformanceSignaler()
{
    //Start an output performance signaler.
    ComputeResources resources = resourceMap.reserveCPUCores(communicator->getHostname(), 0);
    lm::message::Message msg;
    msg.mutable_start_output_performance_signaler()->set_output_interval((int)performanceOutputInterval);
    communicator->sendMessage(resources.controllerAddress, &msg);
}

void JobHypervisor::startWorkUnitRunners()
{
    map<string,ComputeResources> allResources = resourceMap.getAvailableResources();
    slots.createAllSlots(communicator, allResources, cpuCoresPerRunner, gpuDevicesPerRunner, useCPUAffinity, meSolverClassName, diffusionPDESolverClassName);
}

void JobHypervisor::receivedStartedOutputWriter(const lm::message::StartedOutputWriter& msg)
{
    Print::printf(Print::INFO, "Output writer started: %s.", Communicator::printableAddress(msg.address()).c_str());
    hasOutputWriterStarted = true;
    outputWriterAddress = msg.address();
    startJobIfAllWorkersStarted();
}

void JobHypervisor::receivedStartedCheckpointSignaler(const lm::message::StartedCheckpointSignaler& msg)
{
    Print::printf(Print::INFO, "Checkpoint signaller started: %s.", Communicator::printableAddress(msg.address()).c_str());
    hasCheckpointSignalerStarted = true;
    startJobIfAllWorkersStarted();
}

void JobHypervisor::receivedStartedOutputPerformanceSignaler(const lm::message::StartedOutputPerformanceSignaler& msg)
{
    Print::printf(Print::INFO, "Output performance signaller started: %s.", Communicator::printableAddress(msg.address()).c_str());
    hasOutputPerformanceSignalerStarted = true;
    startJobIfAllWorkersStarted();
}

void JobHypervisor::receivedStartedWorkUnitRunner(const lm::message::StartedWorkUnitRunner& msg)
{
    Print::printf(Print::INFO, "Work unit runner started: %s.", Communicator::printableAddress(msg.address()).c_str());

    slots.markSlotStarted(msg);
    if (!slots.hasUnstartedSlots())
    {
        haveAllWorkUnitRunnersStarted = true;
        startJobIfAllWorkersStarted();
    }
}

void JobHypervisor::startJobIfAllWorkersStarted()
{
    if (haveAllWorkUnitRunnersStarted && hasOutputWriterStarted && hasCheckpointSignalerStarted && hasOutputPerformanceSignalerStarted)
    {
        Print::printf(Print::INFO, "All hypervisor workers have started, beginning job.");
        startJob();
    }
}

void JobHypervisor::receivedPerformCheckpointing(const lm::message::PerformCheckpointing& msg)
{
    if (jobRunning)
    {
        Print::printf(Print::INFO, "Started creating a checkpoint, pausing work.");

        // Mark that we are performing a checkpoint.
        performingCheckpoint = true;

        // TODO
    }
}

void JobHypervisor::receivedPerformOutputPerformance(const lm::message::PerformOutputPerformance& msg)
{
    if (jobRunning)
    {
        Print::printf(Print::INFO, "========================================================================================================================");
        Print::printf(Print::INFO, "Performance statistics.");

        // Print the hypervisor stats.
        printPerformanceStatistics();

        // Print slot stats.
        Print::printf(Print::INFO, "Slot status");
        slots.printPerformanceStatistics();

        // Print the output writer stats.
        lm::message::Message msg;
        msg.mutable_perform_output_performance();
        communicator->sendMessage(outputWriterAddress, &msg);
    }
}

void JobHypervisor::receivedFinishedCheckpointing(const lm::message::FinishedCheckpointing& msg)
{
    Print::printf(Print::INFO, "Finished creating a checkpoint, resuming work.");

    // Mark that we are done checkpointing.
    performingCheckpoint = false;


    // TODO

    // Otherwise, see if all outstanding work units have finished.
    //else if (!slotList.hasBusySlots())
    //{
    //    Print::printf(Print::INFO, "Creating a checkpoint, pausing work.");

    //    // Send a message to the output writer to save a checkpoint. Calling .mutable_perform_checkpointing() initializes the submessage
    //    lm::message::Message msgp;
    //    msgp.mutable_perform_checkpointing();
    //    communicator->sendMessage(outputWriterAddress, &msgp);
    //}

//    // Resume distribution of work.
//    if (assignWork())
//    {
//        Print::printf(Print::INFO, "Simulation finished.");
//        finishSimulation();
//    }
}

void JobHypervisor::supervisorStartedEvent(const lm::message::SupervisorStarted& msg)
{
}

void JobHypervisor::supervisorFinishedEvent(const lm::message::SupervisorFinished& msg)
{
    // A supervisor finished, try to start some more.
    if (!startRemainingSimulations())
    {
        // If start supervisors return false, we are done.
        finishJob();
    }
}

void JobHypervisor::startJob()
{
    Print::printf(Print::INFO, "Job started.");

    if (!startRemainingSimulations())
    {
        // If start supervisors return false, there was nothing to be done.
        Print::printf(Print::INFO, "No simulations to start.");
        finishJob();
    }
}

void JobHypervisor::finishJob()
{
    Print::printf(Print::INFO, "Job finished.");

    // Mark that the simulation is finished so we exit our message loop.
    jobRunning = false;

    // Stop all of the resource controllers.
    map<string,ComputeResources> resources = resourceMap.getAvailableResources();
    for (map<string,ComputeResources>::iterator it=resources.begin(); it != resources.end(); it++)
    {
        // Send a message for the resource controller to stop.
        lm::message::Message msg;
        msg.mutable_stop_resource_controller()->set_abort(false);
        communicator->sendMessage(it->second.controllerAddress, &msg);
    }
}

}
}
