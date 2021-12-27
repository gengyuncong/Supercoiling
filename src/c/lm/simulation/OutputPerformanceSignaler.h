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

#ifndef LM_SIMULATION_OUTPUTPERFORMANCESIGNALER
#define LM_SIMULATION_OUTPUTPERFORMANCESIGNALER

#include <ctime>
#include <pthread.h>
#include "lm/message/Communicator.h"
#include "lm/message/Endpoint.pb.h"
#include "lm/thread/Thread.h"
#include "lm/thread/Worker.h"

namespace lm {
namespace simulation {

class OutputPerformanceSignaler : public lm::thread::Worker
{
public:
    OutputPerformanceSignaler(time_t checkpointInterval);
    virtual ~OutputPerformanceSignaler();
    virtual void wake();

protected:
    virtual int run();
    void setNextOutputTime();

private:
    time_t outputInterval;
    lm::message::Communicator* communicator;
    pthread_cond_t controlSignal;
    struct timespec nextOutput;
};

}
}


#endif
