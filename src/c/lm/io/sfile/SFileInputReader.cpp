/*
 * Copyright 2019 Johns Hopkins University
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

#include <string>

#include "lm/String.h"
#include "lm/ClassFactory.h"
#include "lm/input/Input.pb.h"
#include "lm/io/sfile/SFileInputReader.h"
#include "lm/io/sfile/SFile.h"
#include "lm/io/sfile/LocalSFile.h"
#include "lm/types/BoundaryConditions.pb.h"


namespace lm {
namespace io {
namespace sfile {

bool SFileInputReader::registered=SFileInputReader::registerClass();

bool SFileInputReader::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::input::InputReader","lm::io::sfile::SFileInputReader",&SFileInputReader::allocateObject);
    return true;
}

void* SFileInputReader::allocateObject()
{
    return new SFileInputReader();
}

SFileInputReader::SFileInputReader()
{
}

SFileInputReader::~SFileInputReader()
{
}

bool SFileInputReader::canReadFile(std::string filename)
{
    LocalSFile file(filename);
    return (file.exists() && file.isFile() && file.isSFile());
}

void SFileInputReader::readFileInto(std::string filename, lm::input::Input* input)
{
    LocalSFile file(filename);
    file.openRead();
    readInput(file, input);
    file.close();
}

void SFileInputReader::readInput(lm::io::sfile::SFile& file, lm::input::Input* input)
{
    // Read all of the records.
    while (!file.isEof())
    {
        // Read the next record.
        lm::io::sfile::SFileRecord record = file.readNextSFileRecord();

        // See if this is an SimulationInput record.
        if (record.type == "protobuf:lm.input.Input" || record.type == "protobuf:lm.input.SimulationInput")
        {
            // Read the record.
            unsigned char* data = new unsigned char[record.dataSize];
            file.readFully(data, record.dataSize);

            // Parse the record.
            lm::input::Input newInput;
            if (!newInput.ParseFromArray(data, int(record.dataSize))) THROW_EXCEPTION(lm::IOException, "unable to deserialize lm::input::Input record");

            // Copy the parts of the input that are present.
            if (newInput.has_simulation_options()) input->mutable_simulation_options()->CopyFrom(newInput.simulation_options());
            if (newInput.has_output_options()) input->mutable_output_options()->CopyFrom(newInput.output_options());
            if (newInput.has_reaction_model()) input->mutable_reaction_model()->CopyFrom(newInput.reaction_model());
            if (newInput.has_diffusion_model()) input->mutable_diffusion_model()->CopyFrom(newInput.diffusion_model());
            if (newInput.has_microenv_model()) input->mutable_microenv_model()->CopyFrom(newInput.microenv_model());
            if (newInput.has_order_parameters()) input->mutable_order_parameters()->CopyFrom(newInput.order_parameters());
            if (newInput.has_tilings()) input->mutable_tilings()->CopyFrom(newInput.tilings());
            if (newInput.has_ffpilot_options()) input->mutable_ffpilot_options()->CopyFrom(newInput.ffpilot_options());
            if (newInput.has_cme_restart()) input->mutable_cme_restart()->CopyFrom(newInput.cme_restart());
            if (newInput.has_rdme_restart()) input->mutable_rdme_restart()->CopyFrom(newInput.rdme_restart());
        }
        else
        {
            // Skip over the record's data.
            file.skip(record.dataSize);
        }
    }
}

}
}
}
