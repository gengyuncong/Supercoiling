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

#include <cmath>
#include <map>
#include <string>
#include <vector>
#include <google/protobuf/descriptor.h>
#include <google/protobuf/message.h>

#include "lm/Print.h"
#include "lm/String.h"
#include "lm/ClassFactory.h"
#include "lm/input/Input.pb.h"
#include "lm/io/props/PropertiesFile.h"
#include "lm/io/props/PropertiesInputReader.h"
#include "robertslab/Types.h"
#include "robertslab/pbuf/NDArray.pb.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using std::map;
using std::string;
using std::vector;

namespace lm {
namespace io {
namespace props {

bool PropertiesInputReader::registered=PropertiesInputReader::registerClass();

bool PropertiesInputReader::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::input::InputReader", "lm::io::props::PropertiesInputReader", &PropertiesInputReader::allocateObject);
    return true;
}

void* PropertiesInputReader::allocateObject()
{
    return new PropertiesInputReader();
}

PropertiesInputReader::PropertiesInputReader()
{
}

PropertiesInputReader::~PropertiesInputReader()
{
}

bool PropertiesInputReader::canReadFile(std::string filename)
{
    PropertiesFile file(filename);
    return (file.exists() && file.isFile() && (endsWith(filename, ".prop") || endsWith(filename, ".props") || endsWith(filename, ".properties") || endsWith(filename, ".txt")));
}

void PropertiesInputReader::readFileInto(std::string filename, lm::input::Input* input)
{
    PropertiesFile file(filename);
    file.openRead();
    readInput(file, input);
    file.close();
}

void PropertiesInputReader::readInput(lm::io::props::PropertiesFile& file, lm::input::Input* input)
{
    // Read all of the records.
    map<string,string> props = file.readProperties();

    // Set the properties on the input.
    setPropertiesInto(props, input);
}

void PropertiesInputReader::setPropertiesInto(const map<string,string>& props, google::protobuf::Message* msg, bool ignoreUnknown)
{
    for (map<string,string>::const_iterator it=props.begin(); it != props.end(); it++)
    {
        string key = it->first;
        string value = it->second;

        // Split the key into its parts.
        vector<string> keyParts = parseStringList(key, ".");

        // See what to do with the key.
        if (keyParts.size() >= 2 && keyParts[0] == "input")
        {
            keyParts.erase(keyParts.begin());
            if (!setInputProperty(keyParts, value, msg))
                THROW_EXCEPTION(lm::InputException, "Could not parse input property: %s", key.c_str());
        }
        else if (!ignoreUnknown)
        {
            THROW_EXCEPTION(lm::InputException, "Unknown property: %s", key.c_str());
        }
    }
}

bool PropertiesInputReader::setInputProperty(vector<string>& keyParts, string value, google::protobuf::Message* msg)
{
    // Get the next key from the list.
    string key = keyParts[0];

    // Remove the key from the list.
    keyParts.erase(keyParts.begin());

    // See if there are any indices in this key.
    vector<utuple> indices;
    size_t indicesPos = key.find("(");
    if (indicesPos != string::npos)
    {
        indices = parseUintIndices(key.substr(indicesPos));
        key = key.substr(0, indicesPos);
    }

    // Find the field in the message.
    const google::protobuf::FieldDescriptor* fd = NULL;
    const google::protobuf::Descriptor* d = msg->GetDescriptor();
    const google::protobuf::Reflection* r = msg->GetReflection();
    for (int i=0; i<d->field_count(); i++)
    {
        if (d->field(i)->name() == key)
            fd = d->field(i);
    }
    if (fd == NULL) return false;

    // If this is the terminal key, set it.
    if (keyParts.size() == 0)
    {
        return setInputProperty(fd, indices, value, msg);
    }

    // Otherwise, process the next value down the chain.
    else
    {
        // See if this is a singular field.
        if (!fd->is_repeated() && indices.size() == 0)
        {
            return setInputProperty(keyParts, value, r->MutableMessage(msg, fd));
        }

        // See if this field is a repeated field and we have a single index.
        else if (fd->is_repeated() && indices.size() == 1 && indices[0].len == 2 && indices[0][0] == indices[0][1])
        {
            return setInputProperty(keyParts, value, r->MutableRepeatedMessage(msg, fd, static_cast<int>(indices[0][0])));
        }

        // Otherwise, we don't know how to process the field.
        else
        {
            THROW_EXCEPTION(lm::InputException, "Invalid key/field combination: %d %d", fd->is_repeated(), indices.size());
        }
    }
}

bool PropertiesInputReader::setInputProperty(const google::protobuf::FieldDescriptor* fd, vector<utuple> indices, string valuesString, google::protobuf::Message* msg)
{
    // Parse the values.
    vector<string> values = parseStringValues(valuesString);

    const google::protobuf::Reflection* r = msg->GetReflection();

    Print::printf(Print::DEBUG, "Setting field %s of type %d (%s)", fd->name().c_str(), fd->type(), (fd->type() == google::protobuf::FieldDescriptor::Type::TYPE_MESSAGE)?(fd->message_type()->full_name().c_str()):"");

    // See if this is a protobuf field.
    if (fd->type() == google::protobuf::FieldDescriptor::Type::TYPE_MESSAGE && fd->message_type()->full_name() == "robertslab.pbuf.NDArray")
    {
        return setNDArrayInputProperty(r->MutableMessage(msg, fd), indices, values);
    }

    // If this is a singular field.
    else if (!fd->is_repeated() && indices.size() == 0 && values.size() == 1)
    {
        switch (fd->type())
        {
        case google::protobuf::FieldDescriptor::Type::TYPE_FLOAT:
            r->SetFloat(msg, fd, static_cast<float>(atof(values[0].c_str())));
            break;
        case google::protobuf::FieldDescriptor::Type::TYPE_DOUBLE:
            r->SetDouble(msg, fd, static_cast<double>(atof(values[0].c_str())));
            break;
        case google::protobuf::FieldDescriptor::Type::TYPE_INT32:
        case google::protobuf::FieldDescriptor::Type::TYPE_SFIXED32:
        case google::protobuf::FieldDescriptor::Type::TYPE_SINT32:
            r->SetInt32(msg, fd, static_cast<int32_t>(lround(atof(values[0].c_str()))));
            break;
        case google::protobuf::FieldDescriptor::Type::TYPE_UINT32:
        case google::protobuf::FieldDescriptor::Type::TYPE_FIXED32:
            r->SetUInt32(msg, fd, static_cast<uint32_t>(lround(atof(values[0].c_str()))));
            break;
        case google::protobuf::FieldDescriptor::Type::TYPE_INT64:
        case google::protobuf::FieldDescriptor::Type::TYPE_SFIXED64:
        case google::protobuf::FieldDescriptor::Type::TYPE_SINT64:
            r->SetInt64(msg, fd, static_cast<int64_t>(lround(atof(values[0].c_str()))));
            break;
        case google::protobuf::FieldDescriptor::Type::TYPE_UINT64:
        case google::protobuf::FieldDescriptor::Type::TYPE_FIXED64:
            r->SetUInt64(msg, fd, static_cast<uint64_t>(lround(atof(values[0].c_str()))));
            break;
        case google::protobuf::FieldDescriptor::Type::TYPE_BOOL:
            r->SetBool(msg, fd, static_cast<bool>(lround(atof(values[0].c_str()))));
            break;
        case google::protobuf::FieldDescriptor::Type::TYPE_STRING:
            r->SetString(msg, fd, values[0].c_str());
            break;
        default:
            THROW_EXCEPTION(lm::RuntimeException, "Unsupported protobuf field type: %d", fd->type());
        }
        return true;
    }

    // See if this field is a repeated field and we have the same number of indices and values.
    else if (fd->is_repeated() && indices.size() == 1 && indices[0].len == 2 && (indices[0][1]-indices[0][0]+1) == values.size())
    {
        size_t i=0;
        for (int index=static_cast<int>(indices[0][0]); index<=static_cast<int>(indices[0][1]); i++,index++)
        {
            switch (fd->type())
            {
            case google::protobuf::FieldDescriptor::Type::TYPE_FLOAT:
                while (r->FieldSize(*msg, fd) <= index) r->AddDouble(msg, fd, 0.0);
                r->SetRepeatedFloat(msg, fd, index, static_cast<float>(atof(values[i].c_str())));
                break;
            case google::protobuf::FieldDescriptor::Type::TYPE_DOUBLE:
                while (r->FieldSize(*msg, fd) <= index) r->AddDouble(msg, fd, 0.0);
                r->SetRepeatedDouble(msg, fd, index, static_cast<double>(atof(values[i].c_str())));
                break;
            case google::protobuf::FieldDescriptor::Type::TYPE_INT32:
            case google::protobuf::FieldDescriptor::Type::TYPE_SFIXED32:
            case google::protobuf::FieldDescriptor::Type::TYPE_SINT32:
                while (r->FieldSize(*msg, fd) <= index) r->AddInt32(msg, fd, 0);
                r->SetRepeatedInt32(msg, fd, index, static_cast<int32_t>(lround(atof(values[i].c_str()))));
                break;
            case google::protobuf::FieldDescriptor::Type::TYPE_UINT32:
            case google::protobuf::FieldDescriptor::Type::TYPE_FIXED32:
                while (r->FieldSize(*msg, fd) <= index) r->AddUInt32(msg, fd, 0U);
                r->SetRepeatedUInt32(msg, fd, index, static_cast<uint32_t>(lround(atof(values[i].c_str()))));
                break;
            case google::protobuf::FieldDescriptor::Type::TYPE_INT64:
            case google::protobuf::FieldDescriptor::Type::TYPE_SFIXED64:
            case google::protobuf::FieldDescriptor::Type::TYPE_SINT64:
                while (r->FieldSize(*msg, fd) <= index) r->AddInt64(msg, fd, 0);
                r->SetRepeatedInt64(msg, fd, index, static_cast<int64_t>(lround(atof(values[i].c_str()))));
                break;
            case google::protobuf::FieldDescriptor::Type::TYPE_UINT64:
            case google::protobuf::FieldDescriptor::Type::TYPE_FIXED64:
                while (r->FieldSize(*msg, fd) <= index) r->AddUInt64(msg, fd, 0U);
                r->SetRepeatedUInt64(msg, fd, index, static_cast<uint64_t>(lround(atof(values[i].c_str()))));
                break;
            case google::protobuf::FieldDescriptor::Type::TYPE_BOOL:
                while (r->FieldSize(*msg, fd) <= index) r->AddBool(msg, fd, false);
                r->SetRepeatedBool(msg, fd, index, static_cast<bool>(lround(atof(values[i].c_str()))));
                break;
            case google::protobuf::FieldDescriptor::Type::TYPE_STRING:
                while (r->FieldSize(*msg, fd) <= index) r->AddString(msg, fd, "");
                r->SetRepeatedString(msg, fd, index, values[i].c_str());
                break;
            default:
                THROW_EXCEPTION(lm::RuntimeException, "Unsupported protobuf field type: %d", fd->type());
            }
        }
        return true;
    }

    // Otherwise, we don't know how to process the field.
    else
    {
        THROW_EXCEPTION(lm::InputException, "Invalid key/field combination: %d %d", fd->is_repeated(), indices.size());
    }
}

bool PropertiesInputReader::setNDArrayInputProperty(google::protobuf::Message* msg, vector<utuple>& indices, vector<string>& values)
{
    // Copy the message into an NDArray object.
    robertslab::pbuf::NDArray arrayMsg;
    arrayMsg.CopyFrom(*msg);

    // Process the data according to its type.
    if (arrayMsg.data_type() == robertslab::pbuf::NDArray::uint64)
    {
        ndarray<uint64_t> data = robertslab::pbuf::NDArraySerializer::deserialize<uint64_t>(arrayMsg);
        if (indices.size() == 0 && values.size() == 1)
        {
            for (size_t i=0; i<data.size; i++)
                data.values[i] = static_cast<uint64_t>(lround(atof(values[0].c_str())));
        }
        else
        {
            THROW_EXCEPTION(lm::RuntimeException, "reading ndarray properties in that combination currently unimplemented: %d,%d", indices.size(), values.size());
        }
        robertslab::pbuf::NDArraySerializer::serializeInto(&arrayMsg, data);
    }
    else
    {
        THROW_EXCEPTION(lm::RuntimeException, "reading ndarray properties of type %d currently unimplemented: %d", arrayMsg.data_type());
    }

    // Copy the new data back into the message.
    msg->CopyFrom(arrayMsg);

    return true;
}


}
}
}
