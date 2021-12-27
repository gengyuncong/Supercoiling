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


#include <limits>

#include "lm/ffpilot/TrajectoryData.h"

namespace lm {
namespace ffpilot {

TrajectoryData::TrajectoryData()
:parentId(std::numeric_limits<uint64_t>::max()),id(std::numeric_limits<uint64_t>::max()),startingValue(0.0),endingValue(0.0),time(0.0),steps(0)
{
}

TrajectoryData::TrajectoryData(uint64_t parentId, uint64_t id, double startingValue, ndarray<int32_t> startingState)
:parentId(parentId),id(id),startingValue(startingValue),startingState(startingState),endingValue(0.0),time(0.0),steps(0)
{
}

TrajectoryData::TrajectoryData(uint64_t id, double startingValue, ndarray<int32_t> startingState)
:parentId(std::numeric_limits<uint64_t>::max()),id(id),startingValue(startingValue),startingState(startingState),endingValue(0.0),time(0.0),steps(0)
{
}

TrajectoryData::~TrajectoryData()
{
}

}
}
