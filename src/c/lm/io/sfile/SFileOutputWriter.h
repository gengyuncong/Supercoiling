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

#ifndef LM_IO_SFILEOutputWriter
#define LM_IO_SFILEOutputWriter

#include <google/protobuf/message.h>
#include <string>

#include "lm/io/DegreeAdvancementTimeSeries.pb.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/io/LatticeTimeSeries.pb.h"
#include "lm/io/OrderParameterFirstPassageTimes.pb.h"
#include "lm/io/OrderParameterTimeSeries.pb.h"
#include "lm/io/OutputWriter.h"
#include "lm/io/SpeciesCounts.pb.h"
#include "lm/io/SpeciesTimeSeries.pb.h"
#include "lm/io/ffpilot/FFPilotOutput.pb.h"
#include "lm/io/sfile/SFile.h"

namespace lm {
namespace io {
namespace sfile {

class SFileOutputWriter : public OutputWriter
{
public:
    static bool registered;
    static bool registerClass();
    static void* allocateObject();

private:
    static const size_t MAX_NAME_SIZE = 1024;
    static const size_t MAX_TYPE_SIZE = 1024;

public:
    SFileOutputWriter();
    virtual ~SFileOutputWriter();
    virtual void initialize();
    virtual void finalize();

protected:
    SFile* file;

protected:
    virtual void checkpoint();
    virtual void flush();

    virtual void processBarrierCrossingTimes(const std::string& recordNamePrefix, const lm::io::BarrierCrossingTimes& data);
    virtual void processConcentrationsTimeSeries(const std::string& recordNamePrefix, const lm::io::ConcentrationsTimeSeries& data);
    virtual void processDegreeAdvancementTimeSeries(const std::string& recordNamePrefix, const lm::io::DegreeAdvancementTimeSeries& data);
    virtual void processFFPilotOutput(const std::string& recordNamePrefix, const lm::io::ffpilot::FFPilotOutput& data);
    virtual void processLatticeTimeSeries(const std::string& recordNamePrefix, const lm::io::LatticeTimeSeries& data);
    virtual void processOrderParameterFirstPassageTimes(const std::string& recordNamePrefix, const lm::io::OrderParameterFirstPassageTimes& data);
    virtual void processOrderParameterTimeSeries(const std::string& recordNamePrefix, const lm::io::OrderParameterTimeSeries& data);
    virtual void processSpeciesFirstPassageTimes(const std::string& recordNamePrefix, const lm::io::FirstPassageTimes& data);
    virtual void processSpeciesTimeSeries(const std::string& recordNamePrefix, const lm::io::SpeciesTimeSeries& data);
    virtual void processStochasticProcessTimeSeries(const std::string& recordNamePrefix, const lm::io::StochasticProcessTimeSeries& data);

private:
    void processMessage(const string& name, const google::protobuf::Message& data);
};

}
}
}

#endif /* LM_IO_SFILEOutputWriter */
