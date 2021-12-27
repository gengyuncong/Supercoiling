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

#include <list>
#include <string>
#include <vector>

#include "lm/ClassFactory.h"
#include "lm/Print.h"
#include "lm/input/Input.pb.h"
#include "lm/input/InputAggregator.h"
#include "lm/input/InputReader.h"

using std::list;
using std::string;
using std::vector;

namespace lm {
namespace input {

InputAggregator::InputAggregator()
{
}

InputAggregator::~InputAggregator()
{
}

void InputAggregator::readFilesInto(vector<string> inputFilenames, lm::input::Input* input)
{
    for (size_t i = 0; i < inputFilenames.size(); i++)
    {
        readFileInto(inputFilenames[i], input);
    }
}

void InputAggregator::readFileInto(string inputFilename, lm::input::Input* input)
{
    // Go through all of the input readers and see if one can read the file.
    bool processed=false;
    list<string> inputReaderClasses = lm::ClassFactory::getInstance().getAllSubclasses("lm::input::InputReader");
    for (list<string>::iterator it=inputReaderClasses.begin(); it !=inputReaderClasses.end(); it++)
    {
        InputReader* reader = static_cast<lm::input::InputReader*>(lm::ClassFactory::getInstance().allocateObjectOfClass("lm::input::InputReader", *it));
        if (reader->canReadFile(inputFilename))
        {
            Print::printf(Print::INFO, "Loading simulation input from %s using %s.", inputFilename.c_str(), (*it).c_str());
            reader->readFileInto(inputFilename, input);
            processed = true;
            break;
        }
    }

    if (!processed) THROW_EXCEPTION(lm::RuntimeException, "Could not find an input reader to process the file: %s", inputFilename.c_str());
}

}
}
