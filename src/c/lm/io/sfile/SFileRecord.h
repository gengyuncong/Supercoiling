/*
 * Copyright 2014-2019 Johns Hopkins University
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

#ifndef LM_IO_SFILE_SFILERECORD_H
#define LM_IO_SFILE_SFILERECORD_H

#include <string>

namespace lm {
namespace io {
namespace sfile {

class SFileRecord
{
public:
    static const unsigned char RECORD_SEPARATOR[16];

public:
    SFileRecord();
    SFileRecord(const std::string& name, const std::string& type, uint64_t dataSize);

    std::string name;
    std::string type;
    uint64_t dataSize;
};

}
}
}

#endif
