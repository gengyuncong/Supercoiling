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

#include <cerrno>
#include <ctime>
#if defined(MACOSX)
#include <sys/time.h>
#endif
#include <pthread.h>
#include "lm/Print.h"
#include "lm/message/Communicator.h"
#include "lm/message/Endpoint.pb.h"
#include "lm/message/StartedOutputPerformanceSignaler.pb.h"
#include "lm/simulation/OutputPerformanceSignaler.h"
#include "lm/thread/Thread.h"
#include "lm/thread/Worker.h"

using lm::message::Communicator;
using lm::message::Endpoint;
using lm::thread::PthreadException;

namespace lm {
namespace simulation {

OutputPerformanceSignaler::OutputPerformanceSignaler(time_t checkpointInterval)
: outputInterval(checkpointInterval),communicator(NULL)
{
    // Create the communicator.
    communicator = lm::message::Communicator::createObjectOfDefaultSubclass(false);

    PTHREAD_EXCEPTION_CHECK(pthread_cond_init(&controlSignal, NULL));
    nextOutput.tv_sec = 0;
    nextOutput.tv_nsec = 0;
}

OutputPerformanceSignaler::~OutputPerformanceSignaler()
{
    if (communicator != NULL) delete communicator; communicator = NULL;
    pthread_cond_destroy(&controlSignal);
}

void OutputPerformanceSignaler::wake()
{
    //// BEGIN CRITICAL SECTION: controlMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));
    PTHREAD_EXCEPTION_CHECK(pthread_cond_signal(&controlSignal));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex
}

int OutputPerformanceSignaler::run()
{
    try
    {
        Print::printf(Print::INFO, "OutputPerformanceSignaler %s started, writing performance statistics every %d seconds.", Communicator::printableAddress(communicator->getSourceAddress()).c_str(), outputInterval);
        setNextOutputTime();

        // Register our info with the supervisor.
        lm::message::Message msgp;
        lm::message::StartedOutputPerformanceSignaler* msg = msgp.mutable_started_output_performance_signaler();
        msg->mutable_address()->CopyFrom(communicator->getSourceAddress());
        communicator->sendMessage(communicator->getHypervisorAddress(), &msgp);

        bool looping = true;
        while (looping)
        {
            bool doOutput = false;

            //// BEGIN CRITICAL SECTION: controlMutex
            PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));

            // See if we are outputting.
            if (outputInterval > 0)
            {
                // Sleep until our next output interval or we receive a control signal.
                int ret=pthread_cond_timedwait(&controlSignal, &controlMutex, &nextOutput);
				if (ret != 0 && ret != ETIMEDOUT) throw lm::thread::PthreadException(ret,__FILE__,__LINE__);
            }
            else
            {
                // Otherwise we are not checkpointing, so just sleep until our next control change.
                PTHREAD_EXCEPTION_CHECK(pthread_cond_wait(&controlSignal, &controlMutex));
            }

            // If we should abort or have been stopped, stop looping.
            if (aborted || !running)
            {
                looping=false;
            }

            // Otherwise, see if we need to perform a checkpoint.
            else if (outputInterval > 0)
            {
            	// Get the current time.
				#if defined(LINUX)
				struct timespec now;
				clock_gettime(CLOCK_REALTIME, &now);
				#elif defined(MACOSX)
				struct timeval now;
				gettimeofday(&now, NULL);
				#endif

            	// See if we are past the checkpoint time.
                if (now.tv_sec >= nextOutput.tv_sec)
				{
					// Mark that we should perform the checkpoint once we are out of the critical section.
                    doOutput = true;
				}
            }

			PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
			//// END CRITICAL SECTION: controlMutex

			// See if we should perform the checkpoint.
            if (doOutput)
			{
                // Send a output performance message to the hypervisor.
                msgp.Clear();
                msgp.mutable_perform_output_performance();
                communicator->sendMessage(communicator->getHypervisorAddress(), &msgp);

				// Update the next checkpoint time.
                setNextOutputTime();
			}
        }
        Print::printf(Print::INFO, "OutputPerformanceSignaler %s finished.", Communicator::printableAddress(communicator->getSourceAddress()).c_str());
        return 0;
    }
    catch (lm::thread::PthreadException & e)
    {
        Print::printf(Print::FATAL, "Pthread exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (lm::Exception & e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (std::exception & e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (...)
    {
        Print::printf(Print::FATAL, "Unknown Exception during execution (%s:%d)", __FILE__, __LINE__);
    }

    return -1;
}

void OutputPerformanceSignaler::setNextOutputTime()
{
	// Set the next checkpoint time.
	#if defined(LINUX)
	struct timespec now;
	clock_gettime(CLOCK_REALTIME, &now);
    nextOutput.tv_sec = now.tv_sec+outputInterval;
	#elif defined(MACOSX)
	struct timeval now;
	gettimeofday(&now, NULL);
    nextOutput.tv_sec = now.tv_sec+outputInterval;
	#endif
}


}
}
