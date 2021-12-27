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

#ifndef LM_RDME_RDMETRAJECTORY_H_
#define LM_RDME_RDMETRAJECTORY_H_

#include "lm/cme/CMETrajectory.h"
#include "lm/input/DiffusionModel.pb.h"
#include "lm/types/Lattice.pb.h"

namespace lm {
namespace rdme {

class RDMETrajectory : public lm::cme::CMETrajectory
{
public:
    RDMETrajectory(uint64_t id);
    virtual ~RDMETrajectory();

public:
    virtual void initializeLatticeState(const lm::input::DiffusionModel& dm);
    virtual void initializeLatticeParticlesState(const ndarray<uint8_t>& initialParticles);
    virtual void initializeLatticeSitesState(const ndarray<uint8_t>& initialSites);
};

}
}

#endif
