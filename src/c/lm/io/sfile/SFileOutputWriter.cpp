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

#include <cstdio>
#include <string>
#include <google/protobuf/message.h>

#include "lm/ClassFactory.h"
#include "lm/io/BarrierCrossingTimes.pb.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/io/OrderParameterFirstPassageTimes.pb.h"
#include "lm/io/OrderParameterTimeSeries.pb.h"
#include "lm/io/OutputWriter.h"
#include "lm/io/SpeciesTimeSeries.pb.h"
#include "lm/io/ffpilot/FFPilotOutput.pb.h"
#include "lm/io/sfile/LocalSFile.h"
#include "lm/io/sfile/SFileOutputWriter.h"
#include "lm/io/sfile/SFile.h"

using std::string;

namespace lm {
namespace io {
namespace sfile {

bool SFileOutputWriter::registered=SFileOutputWriter::registerClass();

bool SFileOutputWriter::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::io::OutputWriter","lm::io::sfile::SFileOutputWriter",&SFileOutputWriter::allocateObject);
    return true;
}

void* SFileOutputWriter::allocateObject()
{
    return new SFileOutputWriter();
}

SFileOutputWriter::SFileOutputWriter(): file(NULL) {}

SFileOutputWriter::~SFileOutputWriter()
{
    if (file != NULL) delete file; file = NULL;
}

void SFileOutputWriter::initialize()
{
    OutputWriter::initialize();

    // Make sure we have an output filename.
    if (outputFilename == "") throw Exception("Invalid output filename",outputFilename.c_str());

    // Open the file.
    file = new LocalSFile(outputFilename);
    file->openAppend();
}

void SFileOutputWriter::finalize()
{
    OutputWriter::finalize();

    file->close();
    delete file;
    file = NULL;
}

void SFileOutputWriter::checkpoint()
{
}

void SFileOutputWriter::flush()
{
    file->flush();
}

void SFileOutputWriter::processBarrierCrossingTimes(const std::string& recordNamePrefix, const lm::io::BarrierCrossingTimes& data)
{
    char name[MAX_NAME_SIZE];
    snprintf(name, MAX_NAME_SIZE, "%s/Simulations/%lu/BarrierCrossingTimes/%u", recordNamePrefix.c_str(), data.trajectory_id(), data.barrier_index());
    processMessage(name, data);
}

void SFileOutputWriter::processConcentrationsTimeSeries(const std::string& recordNamePrefix, const lm::io::ConcentrationsTimeSeries& data)
{
    char name[MAX_NAME_SIZE];
    snprintf(name, MAX_NAME_SIZE, "%s/Simulations/%lu/ConcentrationsTimeSeries", recordNamePrefix.c_str(), data.trajectory_id());
    processMessage(name, data);
}

void SFileOutputWriter::processDegreeAdvancementTimeSeries(const std::string& recordNamePrefix, const lm::io::DegreeAdvancementTimeSeries& data)
{
    char name[MAX_NAME_SIZE];
    snprintf(name, MAX_NAME_SIZE, "%s/Simulations/%lu/DegreeAdvancementTimeSeries", recordNamePrefix.c_str(), data.trajectory_id());
    processMessage(name, data);
}

void SFileOutputWriter::processFFPilotOutput(const std::string& recordNamePrefix, const lm::io::ffpilot::FFPilotOutput& data)
{
    char name[MAX_NAME_SIZE];
    snprintf(name, MAX_NAME_SIZE, "%s/FFPilotOutput", recordNamePrefix.c_str());
    processMessage(name, data);
}

void SFileOutputWriter::processLatticeTimeSeries(const std::string& recordNamePrefix, const lm::io::LatticeTimeSeries& data)
{
    char name[MAX_NAME_SIZE];
    snprintf(name, MAX_NAME_SIZE, "%s/Simulations/%lu/LatticeTimeSeries", recordNamePrefix.c_str(), data.trajectory_id());
    processMessage(name, data);
}

void SFileOutputWriter::processOrderParameterFirstPassageTimes(const std::string& recordNamePrefix, const lm::io::OrderParameterFirstPassageTimes& data)
{
    char name[MAX_NAME_SIZE];
    snprintf(name, MAX_NAME_SIZE, "%s/Simulations/%lu/OrderParameterFirstPassageTimes", recordNamePrefix.c_str(), data.trajectory_id());
    processMessage(name, data);
}

void SFileOutputWriter::processOrderParameterTimeSeries(const std::string& recordNamePrefix, const lm::io::OrderParameterTimeSeries& data)
{
    char name[MAX_NAME_SIZE];
    snprintf(name, MAX_NAME_SIZE, "%s/Simulations/%lu/OrderParameterTimeSeries", recordNamePrefix.c_str(), data.trajectory_id());
    processMessage(name, data);
}

void SFileOutputWriter::processSpeciesFirstPassageTimes(const std::string& recordNamePrefix, const lm::io::FirstPassageTimes& data)
{
    char name[MAX_NAME_SIZE];
    snprintf(name, MAX_NAME_SIZE, "%s/Simulations/%lu/FirstPassageTimes/%d", recordNamePrefix.c_str(), data.trajectory_id(), data.species());
    processMessage(name, data);
}

void SFileOutputWriter::processSpeciesTimeSeries(const std::string& recordNamePrefix, const lm::io::SpeciesTimeSeries& data)
{
    char name[MAX_NAME_SIZE];
    snprintf(name, MAX_NAME_SIZE, "%s/Simulations/%lu/SpeciesTimeSeries", recordNamePrefix.c_str(), data.trajectory_id());
    processMessage(name, data);
}

void SFileOutputWriter::processStochasticProcessTimeSeries(const std::string& recordNamePrefix, const lm::io::StochasticProcessTimeSeries& data)
{
    char name[MAX_NAME_SIZE];
    snprintf(name, MAX_NAME_SIZE, "%s/Simulations/%lu/StochasticProcessTimeSeries", recordNamePrefix.c_str(), data.trajectory_id());
    processMessage(name, data);
}

void SFileOutputWriter::processMessage(const string& name, const google::protobuf::Message& data)
{
    // Create the type string.
    char type[MAX_TYPE_SIZE];
    snprintf(type, MAX_TYPE_SIZE, "protobuf:%s", data.GetDescriptor()->full_name().c_str());

    // Write the record.
    SFileRecord record(name, type, uint64_t(data.ByteSize()));
    file->writeSFileRecord(record);
    file->writeMessage(data);
}

}
}
}
