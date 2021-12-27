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
#ifndef LM_IO_PROPS_PROPERTIESINPUTREADER_H
#define LM_IO_PROPS_PROPERTIESINPUTREADER_H

#include <string>
#include <vector>
#include <google/protobuf/message.h>

#include "lm/input/Input.pb.h"
#include "lm/input/InputReader.h"
#include "lm/io/props/PropertiesFile.h"
#include "robertslab/Types.h"

namespace lm {
namespace io {
namespace props {

class PropertiesInputReader : public lm::input::InputReader
{
public:
    static bool registered;
    static bool registerClass();
    static void* allocateObject();

public:
    PropertiesInputReader();
    virtual ~PropertiesInputReader();
    virtual bool canReadFile(std::string filename);
    virtual void readFileInto(std::string filename, lm::input::Input* input);
    virtual void setPropertiesInto(const std::map<std::string,std::string>& props, google::protobuf::Message* msg, bool ignoreUnknown=false);

protected:
    virtual void readInput(lm::io::props::PropertiesFile& file, lm::input::Input* input);
    virtual bool setInputProperty(std::vector<std::string>& keyParts, string value, google::protobuf::Message* msg);
    virtual bool setInputProperty(const google::protobuf::FieldDescriptor* fd, std::vector<utuple> indices, string value, google::protobuf::Message* msg);
    virtual bool setNDArrayInputProperty(google::protobuf::Message* msg, std::vector<utuple>& indices, std::vector<string>& values);
};

}
}
}

#endif // PROPERTIESINPUTREADER_H
