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
#include "lm/io/props/PropertiesFile.h"

using std::string;
using std::vector;

namespace lm {
namespace io {
namespace props {

PropertiesFile::PropertiesFile(const string filename)
:File(filename)
{
}

PropertiesFile::~PropertiesFile()
{
}

std::map<std::string,std::string> PropertiesFile::readProperties()
{
    std::map<std::string,std::string> props;

    // Line saved from a continuation.
    string previousLine = "";

    // Loop through the file.
    while (!isEof())
    {
        // Read the line.
        string line = trim(readLine());

        // Remove any comments.
        size_t commentPos = line.find("#");
        if (commentPos != string::npos) line = trim(line.substr(0, commentPos));
        commentPos = line.find("//");
        if (commentPos != string::npos) line = trim(line.substr(0, commentPos));
        commentPos = line.find("!");
        if (commentPos != string::npos) line = trim(line.substr(0, commentPos));

        if (line.length() > 0)
        {
            // Add any previous saved line.
            if (previousLine != "")
            {
                line = previousLine + line;
                previousLine = "";
            }

            // If this line is continued, save it and continue with the next line.
            if (line.at(line.size()-1) == '\\')
            {
                line = line.substr(0, line.size()-1);
                previousLine = line;
                continue;
            }

            // Parse the line.
            vector<string> cols = parseStringList(line, "=");
            if (cols.size() != 2) THROW_EXCEPTION(lm::RuntimeException, "Invalid properties line: %s", line.c_str());

            // Add the key-value pair to the map.
            string key = trim(cols[0]);
            string value = trim(cols[1]);
            props[key] = value;
        }
    }

    return props;
}

}
}
}
