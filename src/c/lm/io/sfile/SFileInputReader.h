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
#ifndef LM_IO_HDF5_SFileInputReader_H
#define LM_IO_HDF5_SFileInputReader_H

#include <string>

#include "lm/input/Input.pb.h"
#include "lm/input/InputReader.h"
#include "lm/io/sfile/SFile.h"
#include "lm/types/OrderParameters.pb.h"

namespace lm {
namespace io {
namespace sfile {

class SFileInputReader : public lm::input::InputReader
{
public:
    static bool registered;
    static bool registerClass();
    static void* allocateObject();

public:
    SFileInputReader();
    virtual ~SFileInputReader();
    virtual bool canReadFile(std::string filename);
    virtual void readFileInto(std::string filename, lm::input::Input* input);

protected:
    virtual void readInput(lm::io::sfile::SFile& file, lm::input::Input* input);
};

}
}
}
#endif // LM_IO_HDF5_SFileInputReader_H
