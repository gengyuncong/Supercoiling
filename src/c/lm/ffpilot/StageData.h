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

#ifndef LM_FFPILOT_STAGE_DATA_H_
#define LM_FFPILOT_STAGE_DATA_H_

#include <limits>

#include "lm/ffpilot/PhaseData.h"

namespace lm {
namespace ffpilot {

class StageData
{
public:
    StageData();
    virtual ~StageData();

public:
    void initialize(size_t numberPhases, bool pilot, bool forward, size_t linkedProductionStage=std::numeric_limits<size_t>::max());
    size_t getNumberPhases();
    PhaseData& getPhaseData(size_t phaseIndex);
    PhaseData& getPhaseData(int phaseIndex);

    void calculateStageStatistics();

protected:
    size_t numberPhases;
    PhaseData* phaseData;

public:
    bool pilot;
    bool forward;
    size_t linkedProductionStage;

    ndarray<double> phaseCosts;
    ndarray<double> phaseWeights;
    ndarray<double> phaseWeightVariances;
};

}
}

#endif
