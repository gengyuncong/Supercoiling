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
 * Author(s): Elijah Roberts, Max Klein
 */

#ifndef LM_THREAD_THREAD_H_
#define LM_THREAD_THREAD_H_

#include <errno.h>
#include <pthread.h>

#include "lm/Exceptions.h"

namespace lm {
namespace thread {

/**
 * PThread exception.
 */
class PthreadException : public Exception
{
public:
    virtual const char* preamble() const throw() {return "Error in pthread library";}

    PthreadException(const int line, const char* file, int err) {init(file, line, err);}
};

class Thread
{
public:
    static int nextThreadNumber;

public:
    Thread();
    virtual ~Thread();
    virtual void start();
    virtual void stop();
    virtual void wait();
    virtual void wake()=0;
    virtual pthread_t getId() {return threadId;}
    virtual int getThreadNumber() {return threadNumber;}
    virtual void setAffinity(int cpuNumber);

protected:
    virtual int run()=0;

private:
    static void * start_thread(void * obj);

protected:
    pthread_mutex_t controlMutex;
    int threadNumber;
    pthread_t threadId;
    volatile bool running;
    int cpuNumber;
};

}
}

/**
 * Exception wrapping for pthreads api.
 */
#define PTHREAD_EXCEPTION_CHECK(pthread_call) {int _pthread_ret_=pthread_call; if (_pthread_ret_ != 0) THROW_EXCEPTION(lm::thread::PthreadException, _pthread_ret_);}
#define PTHREAD_TIMEOUT_EXCEPTION_CHECK(pthread_call) {int _pthread_ret_=pthread_call; if (_pthread_ret_ != 0 && _pthread_ret_ != ETIMEDOUT) THROW_EXCEPTION(lm::thread::PthreadException, _pthread_ret_);}

#endif /* LM_THREAD_THREAD_H_ */
