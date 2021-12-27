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

#ifndef LM_IO_OUTPUTWRITER
#define LM_IO_OUTPUTWRITER

#include <pthread.h>
#include <queue>
#include <string>
#include <vector>

#include "hrtime.h"
#include "lm/io/BarrierCrossingTimes.pb.h"
#include "lm/io/ConcentrationsTimeSeries.pb.h"
#include "lm/io/DegreeAdvancementTimeSeries.pb.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/io/LatticeTimeSeries.pb.h"
#include "lm/io/SpeciesTimeSeries.pb.h"
#include "lm/io/StochasticProcessTimeSeries.pb.h"
#include "lm/io/ffpilot/FFPilotOutput.pb.h"
#include "lm/message/Communicator.h"
#include "lm/message/Message.pb.h"
#include "lm/message/ProcessWorkUnitOutput.pb.h"
#include "lm/thread/Thread.h"
#include "lm/thread/Worker.h"

namespace lm {
namespace io {

class OutputWriter : public lm::thread::Worker
{
public:
    OutputWriter();
    virtual ~OutputWriter();
    virtual void setOutputFilename(std::string outputFilename);
    virtual void initialize();
    virtual void finalize();

    virtual void wake();

protected:
    virtual void checkpoint()=0;
    virtual void flush()=0;

    virtual void processBarrierCrossingTimes(const std::string& recordNamePrefix, const lm::io::BarrierCrossingTimes& data)=0;
    virtual void processConcentrationsTimeSeries(const std::string& recordNamePrefix, const lm::io::ConcentrationsTimeSeries& data)=0;
    virtual void processDegreeAdvancementTimeSeries(const std::string& recordNamePrefix, const lm::io::DegreeAdvancementTimeSeries& data)=0;
    virtual void processFFPilotOutput(const std::string& recordNamePrefix, const lm::io::ffpilot::FFPilotOutput& data)=0;
    virtual void processLatticeTimeSeries(const std::string& recordNamePrefix, const lm::io::LatticeTimeSeries& data)=0;
    virtual void processOrderParameterFirstPassageTimes(const std::string& recordNamePrefix, const lm::io::OrderParameterFirstPassageTimes& data)=0;
    virtual void processOrderParameterTimeSeries(const std::string& recordNamePrefix, const lm::io::OrderParameterTimeSeries& data)=0;
    virtual void processSpeciesFirstPassageTimes(const std::string& recordNamePrefix, const lm::io::FirstPassageTimes& data)=0;
    virtual void processSpeciesTimeSeries(const std::string& recordNamePrefix, const lm::io::SpeciesTimeSeries& data)=0;
    virtual void processStochasticProcessTimeSeries(const std::string& recordNamePrefix, const lm::io::StochasticProcessTimeSeries& data)=0;

    virtual int run();

private:
    static const int MESSAGE_QUEUE_MAX_SIZE=200*1024*1024;

protected:
    std::string outputFilename;

private:
    lm::message::Communicator* communicator;
    std::queue<lm::message::Message*> messageQueue;
    volatile int messageQueueSize;
    pthread_mutex_t messageQueueMutex;
    pthread_cond_t messageQueueSignal;

private:
    class HelperThread : public lm::thread::Thread
    {
    public:
        HelperThread(OutputWriter* p);
        virtual ~HelperThread();
        virtual void wake();
        virtual void printPerformanceStatistics();

    protected:
        virtual int run();

        virtual void processGenericMessage(const google::protobuf::Message& data);
    private:
        OutputWriter* p;
        char* buffer;

        hrtime lastUpdateTime;
        hrtime writingTime;
        long long bytesWritten;
        long long totalBytesWritten;
        long long messagesWritten;
        long long totalMessagesWritten;
    };
};

}
}


#endif
