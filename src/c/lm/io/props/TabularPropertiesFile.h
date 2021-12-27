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


#ifndef LM_IO_PROPS_TABULARPROPERTIESFILE_H
#define LM_IO_PROPS_TABULARPROPERTIESFILE_H

#include <cstdio>
#include <map>
#include <string>
#include <vector>
#include <sys/stat.h>

#include "lm/io/File.h"

using std::string;

namespace lm {
namespace io {
namespace props {

class TabularPropertiesFile : public File
{
public:
    TabularPropertiesFile(const string filename);
    virtual ~TabularPropertiesFile();
    virtual void openRead();
    virtual std::map<std::string,std::string> readNextProperties();
    virtual std::vector<std::string> getKeys();

protected:
    virtual std::vector<std::string> readColumns();

protected:
    bool readHeader;
    std::vector<std::string> keys;
};

}
}
}

#endif
