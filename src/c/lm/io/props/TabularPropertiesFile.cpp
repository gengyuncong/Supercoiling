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

#include <cstdio>
#include <stdexcept>
#include <string>
#include <vector>
#include <sys/stat.h>

#include "lm/Exceptions.h"
#include "lm/String.h"
#include "lm/io/props/TabularPropertiesFile.h"

using std::string;
using std::vector;

namespace lm {
namespace io {
namespace props {

TabularPropertiesFile::TabularPropertiesFile(const string filename)
:File(filename),readHeader(false)
{
}

TabularPropertiesFile::~TabularPropertiesFile()
{
}

void TabularPropertiesFile::openRead()
{
    // Call the base class method.
    File::openRead();

    // Read the header line.
    keys = readColumns();
    readHeader = true;
}


std::map<std::string,std::string> TabularPropertiesFile::readNextProperties()
{
    // Read the next set of columns from the file.
    vector<string> values = readColumns();

    // Make sure we have the right number of values.
    if (values.size() != keys.size()) THROW_EXCEPTION(lm::InputException, "Invalid tabular properties format: %d %d", values.size(), keys.size());

    // Create the key-value map.
    std::map<std::string,std::string> props;
    for (size_t i=0; i<keys.size(); i++)
    {
        props[keys[i]] = values[i];
    }

    return props;
}

vector<string> TabularPropertiesFile::getKeys()
{
    return keys;
}

vector<string> TabularPropertiesFile::readColumns()
{
    // Read until we have processed a non-blank line.
    while (!isEof())
    {
        // Read the next line.
        string line = trim(readLine());

        // See if this is a non-blank line.
        if (line.size() > 0)
        {
            // Split the columns by whitespace and return.
            vector<string> cols = parseStringList(line, " \t");

            // TODO: add support for quotation marks surrounding whitespace.

            return cols;
        }
    }

    return vector<string>();
}


}
}
}
