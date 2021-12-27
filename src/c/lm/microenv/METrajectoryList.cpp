/*
 * Copyright 2016-2019 Johns Hopkins University
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

#include <cstdint>
#include <limits>
#include <list>
#include <map>
#include <string>
#include "hrtime.h"
#include "lm/Print.h"
#include "lm/Types.h"
#include "lm/input/Input.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/io/SpeciesCounts.pb.h"
#include "lm/io/TrajectoryState.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/microenv/METrajectoryList.h"
#include "lm/trajectory/Trajectory.h"
#include "robertslab/Types.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using std::map;
using std::string;

namespace lm {
namespace microenv {

METrajectoryList::METrajectoryList()
{
}

METrajectoryList::~METrajectoryList()
{
}

void METrajectoryList::copySpeciesCountInto(ndarray<int32_t>* counts, uint32_t column, uint32_t speciesId)
{
    for (map<uint64_t,lm::trajectory::Trajectory*>::iterator it=trajectories.begin(); it!=trajectories.end(); it++)
    {
        (*counts)[utuple(uint(it->first),column)] = it->second->getState().cme_state().species_counts().species_count(int(speciesId));
    }
}

void METrajectoryList::copySpeciesCountFrom(const ndarray<int32_t>& counts, uint32_t column, uint32_t speciesId)
{
    for (map<uint64_t,lm::trajectory::Trajectory*>::iterator it=trajectories.begin(); it!=trajectories.end(); it++)
    {
        it->second->getMutableState()->mutable_cme_state()->mutable_species_counts()->set_species_count(int(speciesId), counts[utuple(uint(it->first),column)]);
    }
}

}
}
