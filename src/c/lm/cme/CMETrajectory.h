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
 * Author(s): Elijah Roberts, Max Klein
 */

#ifndef LM_CME_CMETRAJECTORY_H_
#define LM_CME_CMETRAJECTORY_H_

#include "lm/types/OrderParameters.pb.h"
#include "lm/input/CMERestart.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/io/OrderParameterFirstPassageTimes.pb.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/trajectory/Trajectory.h"

namespace lm {
namespace cme {

class CMETrajectory : public lm::trajectory::Trajectory
{
public:
    CMETrajectory(uint64_t id);
    virtual ~CMETrajectory();

public:
    virtual void initializeSpeciesCountsState(const lm::input::ReactionModel& rm);
    virtual void initializeSpeciesCountsState(const ndarray<int32_t>& initialCounts, double time);
    virtual void initializeDegreeAdvancementsState(const lm::input::ReactionModel& rm);
    virtual void initializeOrderParametersState(const lm::types::OrderParameters& op);
    virtual void initializeSpeciesFirstPassageTimeState(uint32_t speciesToTrack);
    virtual void initializeOrderParameterFirstPassageTimeState(uint32_t orderParameterToTrack);

protected:
};

}
}

#endif
