/*
 * Copyright 2012-2019 Johns Hopkins University
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
 * Author(s): Elijah Roberts, Max Klein
 */

#ifndef LM_ME_MESOLVER_H
#define LM_ME_MESOLVER_H

#include <map>
#include <string>
#include <vector>

#include "lm/Types.h"
#include "lm/input/DiffusionModel.pb.h"
#include "lm/main/Solver.h"
#include "lm/types/OrderParameters.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/types/Tilings.pb.h"

using std::map;
using std::string;
using std::vector;

namespace lm {
namespace me {

class MESolver : public lm::main::Solver
{
public:
    MESolver();
    virtual ~MESolver();
    virtual bool needsDiffusionModel()=0;
    virtual bool needsReactionModel()=0;
    virtual void setReactionModel(const lm::input::ReactionModel& rm)=0;
    virtual void setDiffusionModel(const lm::input::DiffusionModel& dm)=0;
    virtual void setOrderParameters(const lm::types::OrderParameters& ops)=0;
    virtual void setTilings(const lm::types::Tilings& ops)=0;
};

}
}

#endif
