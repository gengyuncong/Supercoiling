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

#include <iomanip>
#include <sstream>
#include <string>
#include <sys/stat.h>

#include "lm/ClassFactory.h"
#include "lm/io/OutputWriter.h"
#include "lm/io/hdf5/Hdf5OutputWriter.h"

namespace lm {
namespace io {
namespace hdf5 {

using std::stringstream;
using std::string;

bool Hdf5OutputWriter::registered=Hdf5OutputWriter::registerClass();

bool Hdf5OutputWriter::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::io::OutputWriter","lm::io::hdf5::Hdf5OutputWriter",&Hdf5OutputWriter::allocateObject);
    return true;
}

void* Hdf5OutputWriter::allocateObject()
{
    return new Hdf5OutputWriter();
}

Hdf5OutputWriter::Hdf5OutputWriter()
:file(NULL)
{
}

Hdf5OutputWriter::~Hdf5OutputWriter()
{
    if (file != NULL) delete file; file = NULL;
}

void Hdf5OutputWriter::initialize()
{
    OutputWriter::initialize();

    // Make sure we have an output filename.
    if (outputFilename == "") throw Exception("Invalid output filename",outputFilename.c_str());

    // If the file doesn't exists, create it.
    struct stat fileStats;
    if (stat(outputFilename.c_str(), &fileStats) != 0)
    {
        Hdf5File::create(outputFilename);
    }

    // Open the file.
    file = new Hdf5File(outputFilename);
}

void Hdf5OutputWriter::finalize()
{
    OutputWriter::finalize();

    file->close();
    delete file;
    file = NULL;
}

void Hdf5OutputWriter::checkpoint()
{
    file->checkpoint();
}

void Hdf5OutputWriter::flush()
{
    file->flush();
}

void Hdf5OutputWriter::processBarrierCrossingTimes(const std::string& recordNamePrefix, const lm::io::BarrierCrossingTimes& data)
{
    // Sall the record name method in hdf5 file.
    file->setRecordNamePrefix(recordNamePrefix);

    if (data.times().shape(0) != data.counts().shape(0))
        THROW_EXCEPTION(InvalidArgException,
            "In Hdf5OutputWriter::processBarrierCrossingTimes, "
            "number of rows in time array inconsistent with number of rows in counts array.\n"
            "data.times().shape(0): %d, data.counts().shape(0): %d", data.times().shape(0), data.counts().shape(0));

    char groupName[64];
    snprintf(groupName, 64, "BarrierCrossingTimes/%d", data.barrier_index());
    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), groupName, "SpeciesCounts", data.counts(), false);
    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), groupName, "Times", data.times(), false);
}

void Hdf5OutputWriter::processConcentrationsTimeSeries(const std::string& recordNamePrefix, const lm::io::ConcentrationsTimeSeries& data)
{
    // Sall the record name method in hdf5 file.
    file->setRecordNamePrefix(recordNamePrefix);

    THROW_EXCEPTION(lm::RuntimeException, "not implemented");
}

void Hdf5OutputWriter::processDegreeAdvancementTimeSeries(const std::string& recordNamePrefix, const lm::io::DegreeAdvancementTimeSeries& data)
{
    // Sall the record name method in hdf5 file.
    file->setRecordNamePrefix(recordNamePrefix);

    if (data.times().shape(0) != data.counts().shape(0))
        THROW_EXCEPTION(InvalidArgException,
            "In Hdf5OutputWriter::processDegreeAdvancementTimeSeries, "
            "number of rows in time array inconsistent with number of rows in counts array.\n"
            "data.times().shape(0): %d, data.counts().shape(0): %d", data.times().shape(0), data.counts().shape(0));

    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), "", "DegreeAdvancementCounts", data.counts(), false);
    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), "", "DegreeAdvancementTimes", data.times(), false);
}

void Hdf5OutputWriter::processFFPilotOutput(const std::string& recordNamePrefix, const lm::io::ffpilot::FFPilotOutput& data)
{
    // Sall the record name method in hdf5 file.
    file->setRecordNamePrefix(recordNamePrefix);

    file->setFFPilotOutput(data);
}

void Hdf5OutputWriter::processLatticeTimeSeries(const std::string& recordNamePrefix, const lm::io::LatticeTimeSeries& data)
{
    // Sall the record name method in hdf5 file.
    file->setRecordNamePrefix(recordNamePrefix);

    file->appendLatticeTimeSeries(data.trajectory_id(), data);
}

void Hdf5OutputWriter::processOrderParameterFirstPassageTimes(const std::string& recordNamePrefix, const lm::io::OrderParameterFirstPassageTimes& data)
{
    // Sall the record name method in hdf5 file.
    file->setRecordNamePrefix(recordNamePrefix);

    // construct the relative path to the group we're storing the opfpt datasets in
    std::stringstream ss;
    ss << "OrderParameterFirstPassageTime" << "/";
    ss << std::setfill('0') << std::setw(2) << data.order_parameter_id();

    std::string groupRelativePath(ss.str()), valuesDatasetName("Values"), timesDatasetName("Times");

    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), groupRelativePath, valuesDatasetName, data.order_parameter_value(), false);
    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), groupRelativePath, timesDatasetName, data.first_passage_time(), false);
}

void Hdf5OutputWriter::processOrderParameterTimeSeries(const std::string& recordNamePrefix, const lm::io::OrderParameterTimeSeries& data)
{
    // Sall the record name method in hdf5 file.
    file->setRecordNamePrefix(recordNamePrefix);

    if (data.times().shape(0) != data.values().shape(0))
        THROW_EXCEPTION(InvalidArgException,
            "In Hdf5OutputWriter::processOrderParameterTimeSeries, "
            "number of rows in time array inconsistent with number of rows in value array.\n"
            "data.times().shape(0): %d, data.values().shape(0): %d", data.times().shape(0), data.values().shape(0));

    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), "", "OrderParameterValues", data.values(), false);
    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), "", "OrderParameterTimes", data.times(), false);
}

void Hdf5OutputWriter::processSpeciesFirstPassageTimes(const std::string& recordNamePrefix, const lm::io::FirstPassageTimes& data)
{
    // Sall the record name method in hdf5 file.
    file->setRecordNamePrefix(recordNamePrefix);

    file->setFirstPassageTimes(data.trajectory_id(), data);
}

void Hdf5OutputWriter::processSpeciesTimeSeries(const std::string& recordNamePrefix, const lm::io::SpeciesTimeSeries& data)
{
    // Sall the record name method in hdf5 file.
    file->setRecordNamePrefix(recordNamePrefix);

    if (data.times().shape(0) != data.counts().shape(0))
        THROW_EXCEPTION(InvalidArgException,
            "In Hdf5OutputWriter::processSpeciesTimeSeries, "
            "number of rows in time array inconsistent with number of rows in counts array.\n"
            "data.times().shape(0): %d, data.counts().shape(0): %d", data.times().shape(0), data.counts().shape(0));

    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), "", "SpeciesCounts", data.counts(), false);
    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), "", "SpeciesCountTimes", data.times(), false);
}

void Hdf5OutputWriter::processStochasticProcessTimeSeries(const std::string& recordNamePrefix, const lm::io::StochasticProcessTimeSeries& data)
{
    // Sall the record name method in hdf5 file.
    file->setRecordNamePrefix(recordNamePrefix);

    if (data.times().shape(0) != data.values().shape(0))
        THROW_EXCEPTION(InvalidArgException,
            "In Hdf5OutputWriter::processStochasticProcessTimeSeries, "
            "number of rows in time array inconsistent with number of rows in counts array.\n"
            "data.times().shape(0): %d, data.counts().shape(0): %d", data.times().shape(0), data.values().shape(0));

    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), "", "StochasticProcessValues", data.values(), false);
    file->setDatasetFromNDArrayReplicate(data.trajectory_id(), "", "StochasticProcessTimes", data.times(), false);
}

}
}
}
