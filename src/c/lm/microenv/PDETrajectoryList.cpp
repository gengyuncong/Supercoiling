/*
 * Copyright 2016 Johns Hopkins University
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
#include <limits>
#include <list>
#include <map>
#include <set>
#include <string>

#include "hrtime.h"
#include "lm/Math.h"
#include "lm/Print.h"
#include "lm/Types.h"
#include "lm/input/Input.pb.h"
#include "lm/io/DiffusionPDEState.pb.h"
#include "lm/io/TrajectoryState.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/microenv/PDETrajectoryList.h"
#include "lm/trajectory/Trajectory.h"
#include "lm/trajectory/TrajectoryList.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using std::map;
using std::set;
using std::string;

using lm::trajectory::Trajectory;

namespace lm {
namespace microenv {

PDETrajectoryList::PDETrajectoryList()
{
}

PDETrajectoryList::~PDETrajectoryList()
{
}

void PDETrajectoryList::reconcileDiffusionGrid(double gridElementVolume, ndarray<uint32_t>* cellGridPoints, ndarray<double>* cellVolumes, ndarray<int32_t>* cellCurrentCounts, ndarray<int32_t>* cellFlux, int column)
{
    // Get the trajectory state.
    lm::io::DiffusionPDEState* state = trajectories.begin()->second->getMutableState()->mutable_diffusion_pde_state();

    // Get the current diffusion grid.
    ndarray<double>* concentrations = robertslab::pbuf::NDArraySerializer::deserializeAllocate<double>(state->concentrations(column));

    // Go through each cell.
    for (uint i=0; i<cellGridPoints->shape[0]; i++)
    {
        // Get the grid index of the cell.
        utuple gridIndex((*cellGridPoints)[utuple(i,0U)],(*cellGridPoints)[utuple(i,1U)],(*cellGridPoints)[utuple(i,2U)]);

        // Update the grid with the flux from the cells.
        double newC = concentrations->get(gridIndex) + double((*cellFlux)[utuple(i,uint(column))])/(gridElementVolume*NA);
        if (newC < 0.0 && fabs(newC*gridElementVolume*NA) > 0.01)
            Print::printf(Print::WARNING, "Had a negative concentation after cell %d flux: %e M, %e particles.", i, newC, newC*gridElementVolume*NA);
        concentrations->get(gridIndex) = newC>0.0?newC:0.0;

        // Reset the cell counts with the concentration from the grid.
        (*cellCurrentCounts)[utuple(i,uint(column))] = int32_t(floor(concentrations->get(gridIndex)*((*cellVolumes)[utuple(i)]*NA)));
    }

    // Update the state with the new diffusion grid.
    robertslab::pbuf::NDArraySerializer::serializeInto<double>(state->mutable_concentrations(column), *concentrations);

    // Free the grid.
    delete concentrations;
}

}
}
