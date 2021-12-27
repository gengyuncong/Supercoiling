/*
 * University of Illinois Open Source License
 * Copyright 2012-2014 Roberts Group,
 * All rights reserved.
 * 
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with 
 * the Software without restriction, including without limitation the rights to 
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies 
 * of the Software, and to permit persons to whom the Software is furnished to 
 * do so, subject to the following conditions:
 * 
 * - Redistributions of source code must retain the above copyright notice, 
 * this list of conditions and the following disclaimers.
 * 
 * - Redistributions in binary form must reproduce the above copyright notice, 
 * this list of conditions and the following disclaimers in the documentation 
 * and/or other materials provided with the distribution.
 * 
 * - Neither the names of the Roberts Group, Johns Hopkins University
 * nor the names of its contributors may be used to endorse or promote products
 * derived from this Software without specific prior written permission.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL 
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR 
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR 
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Elijah Roberts
 */

#include <sstream>
#include <string>
#include <vector>
#include "lm/message/Communicator.h"
#include "lm/message/Endpoint.pb.h"
#include "lm/resource/ComputeResources.h"

using std::stringstream;
using std::string;
using std::vector;
using lm::message::Communicator;
using lm::message::Endpoint;

namespace lm {
namespace resource {

string ComputeResources::toString()
{
    stringstream ss;
    ss << "";
    ss << hostname;
    if (cpuCores.size() > 0)
    {
        ss << ":c=";
        for (vector<int>::const_iterator it = cpuCores.begin(); it != cpuCores.end(); ++it)
        {
            if (it !=  cpuCores.begin()) ss << ",";
            ss << *it;
        }
    }

    if (gpuDevices.size() > 0)
    {
        ss << ":g=";
        for (vector<int>::const_iterator it = gpuDevices.begin(); it != gpuDevices.end(); ++it)
        {
            if (it !=  gpuDevices.begin()) ss << ",";
            ss << *it;
        }
    }

    return ss.str();
}
}
}
