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

#ifndef LM_PDE_DIFFUSIONPDETRAJECTORY_H_
#define LM_PDE_DIFFUSIONPDETRAJECTORY_H_

#include "lm/Types.h"
#include "lm/trajectory/Trajectory.h"
#include "robertslab/pbuf/NDArray.pb.h"

namespace lm {
namespace pde {

class DiffusionPDETrajectory : public lm::trajectory::Trajectory
{
public:
    DiffusionPDETrajectory(uint64_t id);
    virtual ~DiffusionPDETrajectory();

public:
    virtual void initializeConcentrations(const robertslab::pbuf::NDArray& initialConcentrations);
    virtual void initializeConcentrations(const ndarray<double>& initialConcentrations);

protected:
};

}
}

#endif
