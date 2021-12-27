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

#ifndef LM_FFPILOT_PHASE_DATA_H_
#define LM_FFPILOT_PHASE_DATA_H_

#include <map>
#include <vector>

#include "lm/ffpilot/TrajectoryData.h"

using std::map;
using std::vector;

namespace lm {
namespace ffpilot {

class PhaseData
{
public:
    PhaseData();
    virtual ~PhaseData();

public:
    void initialize(bool phaseZero, bool forward, size_t phaseLimit);
    ndarray<uint64_t> getCosts();
    ndarray<double> getValues();

public:
    bool phaseZero;
    bool forward;

    // The limit for this phase.
    size_t phaseLimit;

    double goalEdge;

    // The data for the phase's trajectories.
    map<uint64_t,TrajectoryData> trajectoryData;
    vector<size_t> successfulCrossings;
    vector<size_t> failedCrossings;

    // Data used by the fallback tiling method.
    map<uint64_t,int> phase0NumberCrossings;
    map<uint64_t,size_t> phase0LastCrossingDirection;
    map<uint64_t,double> phase0LastCrossingTime;
    map<uint64_t,uint64_t> phase0LastCrossingStep;
    map<uint64_t,ndarray<int32_t>> phase0LastCrossingState;

};

}
}

#endif
