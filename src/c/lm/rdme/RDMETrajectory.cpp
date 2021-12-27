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

#include "lm/rdme/RDMETrajectory.h"
#include "lm/types/Lattice.pb.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using std::map;
using std::set;
using robertslab::pbuf::NDArraySerializer;

namespace lm {
namespace rdme {

RDMETrajectory::RDMETrajectory(uint64_t id)
:CMETrajectory(id)
{
}

RDMETrajectory::~RDMETrajectory()
{
}

void RDMETrajectory::initializeLatticeState(const lm::input::DiffusionModel& dm)
{
    state.mutable_rdme_state()->mutable_lattice()->CopyFrom(dm.initial_lattice());
}

void RDMETrajectory::initializeLatticeParticlesState(const ndarray<uint8_t>& initialParticles)
{
    NDArraySerializer::serializeInto(state.mutable_rdme_state()->mutable_lattice()->mutable_particles(), initialParticles);
}

void RDMETrajectory::initializeLatticeSitesState(const ndarray<uint8_t>& initialSites)
{
    NDArraySerializer::serializeInto(state.mutable_rdme_state()->mutable_lattice()->mutable_sites(), initialSites);
}

}
}
