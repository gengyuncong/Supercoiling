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

#ifndef LM_FFPILOT_FFPILOT_DATA_H_
#define LM_FFPILOT_FFPILOT_DATA_H_

#include "lm/ffpilot/StageData.h"

namespace lm {
namespace ffpilot {

class FFPilotData
{
public:
    FFPilotData();
    virtual ~FFPilotData();

public:
    void initialize(size_t numberStages);
    size_t getNumberStages();
    StageData& getStageData(size_t stageIndex);
    StageData& getStageData(int stageIndex);

protected:
    size_t numberStages;
    StageData* stageData;
};

}
}

#endif
