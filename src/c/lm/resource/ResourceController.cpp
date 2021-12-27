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
#include <list>

#if defined(MACOSX)
#include <sys/sysctl.h>
#elif defined(LINUX)
#include <sys/sysinfo.h>
#endif

#include "lm/Exceptions.h"
#include "lm/ClassFactory.h"
#ifdef OPT_CUDA
#include "lm/Cuda.h"
#endif
#include "lm/io/OutputWriter.h"
#include "lm/resource/ResourceController.h"
#include "lm/message/Communicator.h"
#include "lm/message/Endpoint.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/message/ResourcesAvailable.pb.h"
#include "lm/message/StartOutputWriter.pb.h"
#include "lm/message/StartWorkUnitRunner.pb.h"
#include "lm/message/StartedWorkUnitRunner.pb.h"
#include "lm/simulation/CheckpointSignaler.h"
#include "lm/simulation/OutputPerformanceSignaler.h"
#include "lm/simulation/WorkUnitRunner.h"
#include "lm/thread/WorkerManager.h"
#include "hrtime.h"

using lm::message::Communicator;
using lm::message::Endpoint;
using lm::thread::PthreadException;

namespace lm {
namespace main {

ResourceController::ResourceController()
:communicator(NULL)
{
    // Create the communicator.
    communicator = lm::message::Communicator::createObjectOfDefaultSubclass(false);
}

ResourceController::~ResourceController()
{
    for (std::list<lm::thread::Worker*>::iterator it=workers.begin(); it != workers.end(); it++)
        delete *it;
    workers.clear();
    if (communicator != NULL) delete communicator; communicator = NULL;
}

void ResourceController::wake()
{
    lm::message::Message msg;
    msg.mutable_ping_target()->set_id(0);
    communicator->sendMessage(communicator->getSourceAddress(), &msg);
}

/**
 * Gets the number of physical cpu cores on the system.
 */
std::vector<int> ResourceController::getPhysicalCPUCores()
{
    std::vector<int> cpus;

    // Get the number of processors.
    int numberCPUs=0;
    #if defined(MACOSX)
    uint physicalCpuCores;
    size_t  physicalCpuCoresSize=sizeof(physicalCpuCores);
    sysctlbyname("hw.activecpu",&physicalCpuCores,&physicalCpuCoresSize,NULL,0);
    numberCPUs=(int)physicalCpuCores;
    #elif defined(LINUX)
    numberCPUs=get_nprocs();
    #else
    #error "Unsupported architecture."
    #endif

    // Create a pid entry for each cpu.
    for (int i=0; i<numberCPUs; i++)
        cpus.push_back(i);

    return cpus;
}

std::vector<int> ResourceController::getPhysicalGPUs()
{
    std::vector<int> gpus;

    #ifdef OPT_CUDA
    for (int i=0; i<lm::CUDA::getNumberDevices(); i++)
        gpus.push_back(i);
    #endif

    return gpus;
}

int ResourceController::run()
{
    try
    {
        Print::printf(Print::INFO, "Resource controller %s started.", Communicator::printableAddress(communicator->getSourceAddress()).c_str());

        // Register our info with the supervisor.
        lm::message::Message msg;
        msg.mutable_resources_available()->set_hostname(communicator->getHostname());
        msg.mutable_resources_available()->mutable_controller_address()->CopyFrom(communicator->getSourceAddress());
        std::vector<int> cpus=getPhysicalCPUCores();
        for (std::vector<int>::iterator it = cpus.begin() ; it != cpus.end(); ++it)
            msg.mutable_resources_available()->add_cpu_cores(*it);
        std::vector<int> gpus=getPhysicalGPUs();
        for (std::vector<int>::iterator it = gpus.begin() ; it != gpus.end(); ++it)
            msg.mutable_resources_available()->add_gpu_devices(*it);
        communicator->sendMessage(communicator->getHypervisorAddress(), &msg);

        // Loop reading messages.
        lm::message::Message message;
        while (true)
        {
            // Read the next message.
            communicator->receiveMessage(&message);

            // Do something with the message.
            if (message.start_work_unit_runner_size() > 0)
            {
                for (int i=0; i<message.start_work_unit_runner_size(); i++)
                    startWorkUnitRunner(message.start_work_unit_runner(i));
            }
            else if (message.has_start_output_writer())
            {
                startOutputWriter(message.start_output_writer());
            }
            else if (message.has_start_checkpoint_signaler())
            {
                startCheckpointSignaler(message.start_checkpoint_signaler());
            }
            else if (message.has_start_output_performance_signaler())
            {
                startOutputPerformanceSignaler(message.start_output_performance_signaler());
            }
            else if (message.has_stop_resource_controller())
            {
                stopWorkers(message.stop_resource_controller().abort());
                break;
            }
            else
            {
                Print::printf(Print::ERROR, "Resource controller received an unknown message: {\n%s}",message.DebugString().c_str());
            }

            // Clear the message object so it can be used again.
            message.Clear();
        }

        Print::printf(Print::INFO, "Resource controller %s finished.", Communicator::printableAddress(communicator->getSourceAddress()).c_str());
        return 0;
    }
    catch (lm::Exception& e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (std::exception& e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (...)
    {
        Print::printf(Print::FATAL, "Unknown Exception during execution (%s:%d)", __FILE__, __LINE__);
    }
    return -1;
}

void ResourceController::startWorkUnitRunner(const lm::message::StartWorkUnitRunner& msg)
{
    // Start the work unit runner.
    WorkUnitRunner* runner = new WorkUnitRunner(msg);
    runner->start();
    workers.push_back(runner);
}

void ResourceController::startOutputWriter(const lm::message::StartOutputWriter& msg)
{
    lm::io::OutputWriter* writer = static_cast<lm::io::OutputWriter*>(lm::ClassFactory::getInstance().allocateObjectOfClass("lm::io::OutputWriter",msg.output_writer_class()));
    if (msg.use_cpu_affinity()) writer->setAffinity(msg.cpu());
    writer->setOutputFilename(msg.output_filename());
    writer->initialize();
    writer->start();
    workers.push_back(writer);
}

void ResourceController::startCheckpointSignaler(const lm::message::StartCheckpointSignaler& msg)
{
    lm::simulation::CheckpointSignaler* s = new lm::simulation::CheckpointSignaler(msg.checkpoint_interval());
    s->start();
    workers.push_back(s);
}

void ResourceController::startOutputPerformanceSignaler(const lm::message::StartOutputPerformanceSignaler& msg)
{
    lm::simulation::OutputPerformanceSignaler* s = new lm::simulation::OutputPerformanceSignaler(msg.output_interval());
    s->start();
    workers.push_back(s);
}

void ResourceController::stopWorkers(bool abort)
{
    Print::printf(Print::DEBUG, "Resource controller %s stopping workers.", Communicator::printableAddress(communicator->getSourceAddress()).c_str());
    for (std::list<lm::thread::Worker*>::reverse_iterator it=workers.rbegin(); it!=workers.rend(); it++)
	{
        if (abort)
            (*it)->abort();
        else
            (*it)->stop();
	}
}

}
}
