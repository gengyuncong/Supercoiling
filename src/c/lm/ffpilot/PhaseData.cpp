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


#include <map>
#include <limits>

#include "lm/Math.h"
#include "lm/ffpilot/TrajectoryData.h"
#include "lm/ffpilot/PhaseData.h"

using std::map;


namespace lm {
namespace ffpilot {

PhaseData::PhaseData()
:phaseZero(false),forward(true),phaseLimit(0),goalEdge(0.0)
{
}

PhaseData::~PhaseData()
{
}

void PhaseData::initialize(bool phaseZero, bool forward, size_t phaseLimit)
{
    this->phaseZero = phaseZero;
    this->forward = forward;
    this->phaseLimit = phaseLimit;
}

ndarray<uint64_t> PhaseData::getCosts()
{
    ndarray<uint64_t> costs(utuple(static_cast<uint>(trajectoryData.size())));
    uint i=0;
    for (map<uint64_t,TrajectoryData>::iterator it = trajectoryData.begin(); it != trajectoryData.end(); it++, i++)
    {
        costs[i] = it->second.steps;
    }

    return costs;
}

ndarray<double> PhaseData::getValues()
{
    ndarray<double> values(utuple(static_cast<uint>(trajectoryData.size())));
    uint i=0;
    if (phaseZero)
    {
        for (map<uint64_t,TrajectoryData>::iterator it = trajectoryData.begin(); it != trajectoryData.end(); it++, i++)
        {
            values[i] = it->second.time;
        }
    }
    else
    {
        for (map<uint64_t,TrajectoryData>::iterator it = trajectoryData.begin(); it != trajectoryData.end(); it++, i++)
        {
            if (forward)
                values[i] = it->second.endingValue >= goalEdge;
            else
                values[i] = it->second.endingValue <= goalEdge;
        }
    }

    return values;
}

}
}
