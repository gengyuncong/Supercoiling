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

#include "unistd.h"
#include <string>

#include "lm/ClassFactory.h"
#include "lm/io/BarrierCrossingTimes.pb.h"
#include "lm/io/NullOutputWriter.h"
#include "lm/io/OutputWriter.h"


namespace lm {
namespace io {


bool NullOutputWriter::registered=NullOutputWriter::registerClass();

bool NullOutputWriter::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::io::OutputWriter","lm::io::NullOutputWriter",&NullOutputWriter::allocateObject);
    return true;
}

void* NullOutputWriter::allocateObject()
{
    return new NullOutputWriter();
}

NullOutputWriter::NullOutputWriter()
:secondsToDelay(0)
{
}

NullOutputWriter::~NullOutputWriter()
{
}

void NullOutputWriter::checkpoint()
{
}

void NullOutputWriter::flush()
{
}

void NullOutputWriter::processBarrierCrossingTimes(const std::string& recordNamePrefix, const lm::io::BarrierCrossingTimes& data)
{
    processMessage();
}

void NullOutputWriter::processConcentrationsTimeSeries(const std::string& recordNamePrefix, const lm::io::ConcentrationsTimeSeries& data)
{
    processMessage();
}

void NullOutputWriter::processDegreeAdvancementTimeSeries(const std::string& recordNamePrefix, const lm::io::DegreeAdvancementTimeSeries& data)
{
    processMessage();
}

void NullOutputWriter::processFFPilotOutput(const std::string& recordNamePrefix, const lm::io::ffpilot::FFPilotOutput& data)
{
    processMessage();
}

void NullOutputWriter::processLatticeTimeSeries(const std::string& recordNamePrefix, const lm::io::LatticeTimeSeries& data)
{
    processMessage();
}

void NullOutputWriter::processOrderParameterFirstPassageTimes(const std::string& recordNamePrefix, const lm::io::OrderParameterFirstPassageTimes& data)
{
    processMessage();
}

void NullOutputWriter::processOrderParameterTimeSeries(const std::string& recordNamePrefix, const lm::io::OrderParameterTimeSeries& data)
{
    processMessage();
}

void NullOutputWriter::processSpeciesFirstPassageTimes(const std::string& recordNamePrefix, const lm::io::FirstPassageTimes& data)
{
    processMessage();
}

void NullOutputWriter::processSpeciesTimeSeries(const std::string& recordNamePrefix, const lm::io::SpeciesTimeSeries& data)
{
    processMessage();
}

void NullOutputWriter::processStochasticProcessTimeSeries(const std::string& recordNamePrefix, const lm::io::StochasticProcessTimeSeries& data)
{
    processMessage();
}

void NullOutputWriter::processMessage()
{
    if (secondsToDelay > 0)
        sleep(secondsToDelay);
}

}
}
