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


#include "lm/ffpilot/FFPilotData.h"
#include "lm/ffpilot/StageData.h"

namespace lm {
namespace ffpilot {

FFPilotData::FFPilotData()
:numberStages(0),stageData(NULL)
{
}

FFPilotData::~FFPilotData()
{
    if (stageData != NULL) delete[] stageData; stageData = NULL;
    numberStages = 0;
}

void FFPilotData::initialize(size_t numberStages)
{
    this->numberStages = numberStages;
    if (stageData != NULL) THROW_EXCEPTION(lm::RuntimeException, "attempted to reinitialize an ffpilot data");
    stageData = new StageData[numberStages];
}

size_t FFPilotData::getNumberStages()
{
    return numberStages;
}

StageData& FFPilotData::getStageData(size_t stageIndex)
{
    if (stageData == NULL) THROW_EXCEPTION(lm::RuntimeException, "stage data not initialized");
    if (stageIndex >= numberStages) THROW_EXCEPTION(lm::RuntimeException, "stage index %d greater than the number of phases %d", stageIndex, numberStages);

    return stageData[stageIndex];
}

StageData& FFPilotData::getStageData(int stageIndex)
{
    return getStageData(static_cast<size_t>(stageIndex));
}

}
}
