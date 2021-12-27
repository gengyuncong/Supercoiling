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

#ifndef LM_FFPILOT_TRAJECTORY_DATA_H_
#define LM_FFPILOT_TRAJECTORY_DATA_H_

#include "lm/Types.h"

namespace lm {
namespace ffpilot {

class TrajectoryData
{
public:
    TrajectoryData();
    TrajectoryData(uint64_t parentId, uint64_t id, double startingValue, ndarray<int32_t> startingState);
    TrajectoryData(uint64_t id, double startingValue, ndarray<int32_t> startingState);
    virtual ~TrajectoryData();

public:
    uint64_t parentId;
    uint64_t id;
    double startingValue;
    ndarray<int32_t> startingState;
    double endingValue;
    ndarray<int32_t> endingState;
    double time;
    uint64_t steps;
};

}
}

#endif
