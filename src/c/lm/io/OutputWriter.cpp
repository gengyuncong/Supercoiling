/*
 * Copyright 2012-2019 Johns Hopkins University
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
#include <google/protobuf/message.h>
#include <queue>
#include <pthread.h>
#include <sstream>
#include <sys/time.h>
#include <time.h>

#include "hrtime.h"
#include "lm/ClassFactory.h"
#include "lm/Math.h"
#include "lm/Print.h"
#include "lm/String.h"
#include "lm/Types.h"
#include "lm/io/OutputWriter.h"
#include "lm/main/Globals.h"
#include "lm/simulation/SimulationSupervisor.h"
#include "lm/message/Communicator.h"
#include "lm/message/Endpoint.pb.h"
#include "lm/message/FinishedCheckpointing.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/message/ProcessAggregatedOutput.pb.h"
#include "lm/message/ProcessWorkUnitOutput.pb.h"
#include "lm/message/StartedOutputWriter.pb.h"
#include "lm/message/WorkUnitOutput.pb.h"
#include "lm/thread/Thread.h"
#include "lm/thread/Worker.h"
#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"

namespace lm {
namespace io {

using std::string;
using std::stringstream;
using std::vector;

using lm::message::Communicator;
using lm::message::Endpoint;

OutputWriter::OutputWriter()
:outputFilename(""),communicator(NULL),messageQueueSize(0)
{
    // Create the communicator.
    communicator = lm::message::Communicator::createObjectOfDefaultSubclass(false);

    // Create the queue mutex.
    pthread_mutexattr_t attr;
    PTHREAD_EXCEPTION_CHECK(pthread_mutexattr_init(&attr));
    PTHREAD_EXCEPTION_CHECK(pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_NORMAL));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_init(&messageQueueMutex, &attr));
    PTHREAD_EXCEPTION_CHECK(pthread_mutexattr_destroy(&attr));

    // Create the queue signal.
    PTHREAD_EXCEPTION_CHECK(pthread_cond_init(&messageQueueSignal, NULL));
}

OutputWriter::~OutputWriter()
{
    if (communicator != NULL) delete communicator; communicator = NULL;
    pthread_mutex_destroy(&messageQueueMutex);
    pthread_cond_destroy(&messageQueueSignal);
}

void OutputWriter::setOutputFilename(std::string outputFilename)
{
    this->outputFilename = outputFilename;
}

void OutputWriter::initialize()
{
}

void OutputWriter::finalize()
{
}

void OutputWriter::wake()
{
    lm::message::Message msg;
    msg.mutable_ping_target()->set_id(0);
    communicator->sendMessage(communicator->getSourceAddress(), &msg);
}

int OutputWriter::run()
{
    PROF_BEGIN(PROF_DATAOUTPUT_RUN);

    // Create a helper thread.
    HelperThread helperThread(this);

    try
    {
        Print::printf(Print::INFO, "OutputWriter %s started.", Communicator::printableAddress(communicator->getSourceAddress()).c_str());

        // Start the helper thread.
        if (cpuNumber >= 0) helperThread.setAffinity(cpuNumber);
        helperThread.start();

        // Register our info with the supervisor.
        lm::message::Message msgp;
        lm::message::StartedOutputWriter* msg = msgp.mutable_started_output_writer();
        msg->mutable_address()->CopyFrom(communicator->getSourceAddress());
        communicator->sendMessage(communicator->getHypervisorAddress(), &msgp);

        // Loop reading messages.
        while (true)
        {
            // Read the next message.
            lm::message::Message* message = new lm::message::Message();
            communicator->receiveMessage(message);

            if (message->has_process_work_unit_output() || message->has_process_aggregated_output())
            {
                //// BEGIN CRITICAL SECTION: messageQueueMutex
                PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&messageQueueMutex));

                // Add the message to the queue.
                messageQueue.push(message);

                // Track the size of the queue.
                messageQueueSize += message->ByteSize();
                int tmpMessageQueueSize = messageQueueSize;

                // Signal that data is aavailable.
                PTHREAD_EXCEPTION_CHECK(pthread_cond_signal(&messageQueueSignal));

                PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&messageQueueMutex));
                //// END CRITICAL SECTION: messageQueueMutex

                // If the queue is too full, wait until it empties before reading any more messages.
                while (tmpMessageQueueSize > MESSAGE_QUEUE_MAX_SIZE)
                {
                    Print::printf(Print::WARNING, "OutputWriter is receiving too much data, performance may be degraded. If this this message appear frequently, increase write intervals to increase performance. (%d bytes queued)",tmpMessageQueueSize);
                    sleep(5);

                    //// BEGIN CRITICAL SECTION: messageQueueMutex
                    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&messageQueueMutex));
                    tmpMessageQueueSize = messageQueueSize;
                    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&messageQueueMutex));
                    //// END CRITICAL SECTION: messageQueueMutex
                }
            }
            else if (message->has_perform_checkpointing())
            {
                // Wait for the queue to be empty.
                while (true)
                {
                    size_t tmpMessageQueueLength;

                    //// BEGIN CRITICAL SECTION: messageQueueMutex
                    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&messageQueueMutex));
                    tmpMessageQueueLength = messageQueue.size();
                    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&messageQueueMutex));
                    //// END CRITICAL SECTION: messageQueueMutex

                    if (tmpMessageQueueLength > 0)
                        sleep(1);
                    else
                        break;
                }

                // Perform the checkpointing.
                checkpoint();

                // Report back to the supervisor that the checkpoint is finished. Calling .mutable_finished_checkpointing() initializes the message
                msgp.Clear();
                msgp.mutable_finished_checkpointing();
                communicator->sendMessage(communicator->getHypervisorAddress(), &msgp);
            }
            else if (message->has_perform_output_performance())
            {
                helperThread.printPerformanceStatistics();
            }
            else if (message->has_ping_target())
            {
                // If we are done running, stop the loop.
                if (!running) break;
            }
            else
            {
                Print::printf(Print::ERROR, "OutputWriter received an unknown message: {\n%s}",message->DebugString().c_str());
            }
        }

        // Stop the helper thread.
        helperThread.stop();

        // Let the output writer close any resources.
        Print::printf(Print::INFO, "OutputWriter %s flushing and closing.", Communicator::printableAddress(communicator->getSourceAddress()).c_str());
        flush();
        finalize();

        Print::printf(Print::INFO, "OutputWriter %s finished.", Communicator::printableAddress(communicator->getSourceAddress()).c_str());
        PROF_END(PROF_DATAOUTPUT_RUN);
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
    PROF_END(PROF_DATAOUTPUT_RUN);

    // Stop the helper thread.
    helperThread.stop();

    return -1;
}

OutputWriter::HelperThread::HelperThread(OutputWriter* p)
:p(p),buffer(new char[MEBI+1]),
lastUpdateTime(getHrTime()),writingTime(0),bytesWritten(0),totalBytesWritten(0),messagesWritten(0),totalMessagesWritten(0)
{
}

OutputWriter::HelperThread::~HelperThread()
{
    if (buffer != NULL) delete[] buffer; buffer = NULL;
}

void OutputWriter::HelperThread::wake()
{
    //// BEGIN CRITICAL SECTION: messageQueueMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&p->messageQueueMutex));
    PTHREAD_EXCEPTION_CHECK(pthread_cond_signal(&p->messageQueueSignal));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&p->messageQueueMutex));
    //// END CRITICAL SECTION: messageQueueMutex
}

int OutputWriter::HelperThread::run()
{
    try
    {
        Print::printf(Print::INFO, "OutputWriter::HelperThread %s started.", Communicator::printableAddress(p->communicator->getSourceAddress()).c_str());

        // Performance stats.

        bool finished = false;
        while (!finished)
        {
            lm::message::Message* message=NULL;
            int messageSize=0;

            //// BEGIN CRITICAL SECTION: messageQueueMutex
            PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&p->messageQueueMutex));

            // See if there are any messages.
            if (!p->messageQueue.empty())
            {
                // Get the next message.
                message = p->messageQueue.front();
                p->messageQueue.pop();

                // Update the total message size in the queue.
                messageSize = message->ByteSize();
                p->messageQueueSize -= messageSize;
            }

            // If we are not running and the queue is empty, stop after this iteration of the loop.
            else if (!running)
            {
                finished = true;
            }

            // Otherwise, wait for more data.
            else
            {
                struct timeval tv;
                gettimeofday(&tv, NULL);
                struct timespec waitTime;
                waitTime.tv_sec = tv.tv_sec + 60;
                // have to set nanoseconds to avoid pthread_cond_timedwait throwing EINVAL under weirdly specific conditions
                waitTime.tv_nsec = 0;
                PTHREAD_TIMEOUT_EXCEPTION_CHECK(pthread_cond_timedwait(&p->messageQueueSignal, &p->messageQueueMutex, &waitTime));
            }

            PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&p->messageQueueMutex));
            //// END CRITICAL SECTION: messageQueueMutex

            // If we got a message off of the queue, process it.
            if (message != NULL)
            {
                // Loop over every output in the message.
                hrtime startWriting = getHrTime();

                if (message->has_process_work_unit_output())
                {
                    const lm::message::ProcessWorkUnitOutput& pwu = message->process_work_unit_output();
                    for (int i=0; i<pwu.part_output_size(); i++)
                    {
                        const lm::message::WorkUnitOutput& output = pwu.part_output(i);

                        // set the output options
                        const string& recordNamePrefix = output.record_name_prefix();

                        // process the actual output
                        if (output.has_concentrations_time_series())
                        {
                            p->processConcentrationsTimeSeries(recordNamePrefix, output.concentrations_time_series());
                        }
                        if (output.has_degree_advancement_time_series())
                        {
                            p->processDegreeAdvancementTimeSeries(recordNamePrefix, output.degree_advancement_time_series());
                        }
                        if (output.has_lattice_time_series())
                        {
                            p->processLatticeTimeSeries(recordNamePrefix, output.lattice_time_series());
                        }
                        if (output.has_order_parameter_time_series())
                        {
                            p->processOrderParameterTimeSeries(recordNamePrefix, output.order_parameter_time_series());
                        }
                        if (output.has_species_time_series())
                        {
                            p->processSpeciesTimeSeries(recordNamePrefix, output.species_time_series());
                        }
                        if (output.has_stochastic_process_time_series())
                        {
                            p->processStochasticProcessTimeSeries(recordNamePrefix, output.stochastic_process_time_series());
                        }
                        if (output.barrier_crossing_times().size() > 0)
                        {
                            for (int j=0; j<output.barrier_crossing_times().size(); j++)
                                p->processBarrierCrossingTimes(recordNamePrefix, output.barrier_crossing_times(j));
                        }
                        if (output.order_parameter_first_passage_times_size() > 0)
                        {
                            for (int j=0; j<output.order_parameter_first_passage_times_size(); j++)
                                p->processOrderParameterFirstPassageTimes(recordNamePrefix, output.order_parameter_first_passage_times(j));
                        }
                        if (output.first_passage_times_size() > 0)
                        {
                            for (int j=0; j<output.first_passage_times_size(); j++)
                                p->processSpeciesFirstPassageTimes(recordNamePrefix, output.first_passage_times(j));
                        }
                    }
                }
                if (message->has_process_aggregated_output())
                {
                    const lm::message::ProcessAggregatedOutput& output = message->process_aggregated_output();
                    const string& recordNamePrefix = output.record_name_prefix();
                    if (output.has_ffpilot_output())
                    {
                        p->processFFPilotOutput(recordNamePrefix, output.ffpilot_output());
                    }
                }

                // Track some performance statistics.
                writingTime += getHrTime()-startWriting;
                bytesWritten += messageSize;
                messagesWritten++;
                totalBytesWritten += messageSize;
                totalMessagesWritten++;

                // Delete the message.
                delete message;
                message = NULL;
            }

            //TODO flush the data if it has been a while
            //p->flush();

        }

        // Print the final statistics.
        Print::printf(Print::INFO, "OutputWriter wrote %lld messages and %lld bytes total.", totalMessagesWritten, totalBytesWritten);
        Print::printf(Print::INFO, "OutputWriter::HelperThread %s finished.", Communicator::printableAddress(p->communicator->getSourceAddress()).c_str());
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

    exit(-1);
}

void OutputWriter::HelperThread::printPerformanceStatistics()
{
    //// BEGIN CRITICAL SECTION: messageQueueMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&p->messageQueueMutex));

    // See if we should display some stats.
    hrtime currentTime = getHrTime();
    Print::printf(Print::INFO, "OutputWriter wrote %u messages (%lld bytes) in the last %0.1f seconds (%0.6f seconds writing). %u messages (%d bytes) queued.",messagesWritten,bytesWritten,convertHrToSeconds(currentTime-lastUpdateTime), convertHrToSeconds(writingTime), p->messageQueue.size(), p->messageQueueSize);

    lastUpdateTime = currentTime;
    writingTime = 0;
    messagesWritten = 0;
    bytesWritten = 0;

    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&p->messageQueueMutex));
    //// END CRITICAL SECTION: messageQueueMutex
}

void OutputWriter::HelperThread::processGenericMessage(const google::protobuf::Message& data)
{
    stringstream debugSS;
    debugSS << "--------------------------------------------------------------------------------\n";
    debugSS << data.DebugString();
    debugSS << "--------------------------------------------------------------------------------";
    Print::printf(Print::INFO, "OutputWriter received %s:\n%s", data.GetDescriptor()->name().c_str(), debugSS.str().c_str());
}

}
}
