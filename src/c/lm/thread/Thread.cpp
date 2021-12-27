/*
 * University of Illinois Open Source License
 * Copyright 2008-2012 Luthey-Schulten Group,
 * Copyright 2012-2016 Roberts Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 *               University of Illinois at Urbana-Champaign
 *               http://www.scs.uiuc.edu/~schulten
 * 
 * Developed by: Roberts Group
 *               Johns Hopkins University
 *               http://biophysics.jhu.edu/roberts/
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
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Elijah Roberts
 */

#include <pthread.h>
#include "lm/Print.h"
#include "lm/thread/Thread.h"
#include "lptf/Profile.h"

namespace lm {
namespace thread {

int Thread::nextThreadNumber=1;

Thread::Thread()
:threadNumber(0),threadId(0),running(false),cpuNumber(-1)
{
    // Assign the thread number.
    threadNumber = nextThreadNumber++; //TODO: add a global mutex lock.

    // Create the control mutex.
	pthread_mutexattr_t attr;
	PTHREAD_EXCEPTION_CHECK(pthread_mutexattr_init(&attr));
	PTHREAD_EXCEPTION_CHECK(pthread_mutexattr_settype(&attr, PTHREAD_MUTEX_NORMAL));
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_init(&controlMutex, &attr));
    PTHREAD_EXCEPTION_CHECK(pthread_mutexattr_destroy(&attr));
}

Thread::~Thread()
{
    pthread_mutex_destroy(&controlMutex);
}

void Thread::setAffinity(int cpuNumber)
{
    //// BEGIN CRITICAL SECTION: controlMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));

    // Set the assigned cpu.
	this->cpuNumber = cpuNumber;

	// If we are already running, reset the processor affinity.
    if (running)
    {
		#if defined(LINUX)
		cpu_set_t cpuset;
		CPU_ZERO(&cpuset);
		CPU_SET(cpuNumber, &cpuset);
		if (pthread_setaffinity_np(threadId, sizeof(cpu_set_t), &cpuset) != 0)
            Print::printf(Print::WARNING, "Could not bind thread %d to CPU core %d", threadNumber, cpuNumber);
		#endif
    }
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex
}

void Thread::start()
{
    //// BEGIN CRITICAL SECTION: controlMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));
    if (!running)
    {
        running=true;
        pthread_attr_t attr;
        PTHREAD_EXCEPTION_CHECK(pthread_attr_init(&attr));
        PTHREAD_EXCEPTION_CHECK(pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE));
        PTHREAD_EXCEPTION_CHECK(pthread_create(&threadId, &attr, &Thread::start_thread, this));
        PTHREAD_EXCEPTION_CHECK(pthread_attr_destroy(&attr));

    	// Set the processor affinity, if we have a cpu assigned.
        if (cpuNumber >= 0)
        {
			#if defined(LINUX)
			cpu_set_t cpuset;
			CPU_ZERO(&cpuset);
			CPU_SET(cpuNumber, &cpuset);
			if (pthread_setaffinity_np(threadId, sizeof(cpu_set_t), &cpuset) != 0)
                Print::printf(Print::WARNING, "Could not bind thread %d to CPU core %d", threadNumber, cpuNumber);
			#endif
        }
        Print::printf(Print::DEBUG, "Started thread %d.", threadNumber);
    }
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex
}

void * Thread::start_thread(void * obj)
{
    // Get the thread object.
    Thread* thread=(reinterpret_cast<Thread *>(obj));

    // Set the thread number in the profiler.
    PROF_SET_THREAD(thread->threadNumber);

    // Enter the thread run method.
    size_t ret = thread->run();

    // The run method has finished, so the thread should exit.
    pthread_exit((void *)ret);
}

void Thread::stop()
{
    bool waitForThread=false;

    //// BEGIN CRITICAL SECTION: controlMutex
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_lock(&controlMutex));
    if (running)
    {
        Print::printf(Print::DEBUG, "Stopping thread %d (0x%X).", threadNumber,threadId);
        running = false;
        waitForThread = true;
    }
    PTHREAD_EXCEPTION_CHECK(pthread_mutex_unlock(&controlMutex));
    //// END CRITICAL SECTION: controlMutex

    // Wait for the thread to stop, if we really stopped it.
    if (waitForThread)
    {
        // Try to wake up the thread in case it is sleeping.
        wake();

        // Join with the thread.
        void * ret;
        PTHREAD_EXCEPTION_CHECK(pthread_join(threadId, &ret));
        Print::printf(Print::DEBUG, "Thread %d stopped.", threadNumber);
    }
}

void Thread::wait()
{
    // Join with the thread.
    void * ret;
    PTHREAD_EXCEPTION_CHECK(pthread_join(threadId, &ret));
}

}
}
