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

#include <cmath>
#include <map>
#include <set>

#include "lm/Exceptions.h"
#include "lm/Types.h"
#include "lm/pde/DiffusionPDETrajectory.h"
#include "lm/trajectory/Trajectory.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using std::map;
using std::set;
using robertslab::pbuf::NDArraySerializer;

namespace lm {
namespace pde {

DiffusionPDETrajectory::DiffusionPDETrajectory(uint64_t id)
:Trajectory(id)
{
}

DiffusionPDETrajectory::~DiffusionPDETrajectory()
{
}

void DiffusionPDETrajectory::initializeConcentrations(const robertslab::pbuf::NDArray& initialConcentrations)
{
    state.mutable_diffusion_pde_state()->set_time(0.0);
    state.mutable_diffusion_pde_state()->add_concentrations()->CopyFrom(initialConcentrations);
}

void DiffusionPDETrajectory::initializeConcentrations(const ndarray<double>& initialConcentrations)
{
    THROW_EXCEPTION(lm::RuntimeException, "unimplemented");
}

}
}
