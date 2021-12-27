/*
 * University of Illinois Open Source License
 * Copyright 2012-2016 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy of
 * this software and associated documentation files (the Software), to deal with
 * the Software without restriction, including without limitation the rights to
 * use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
 * of the Software, and to permit persons to whom the Software is furnished to
 * do so, subject to the following conditions:
 *
 * - Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimers.
 *
 * - Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimers in the documentation
 * and/or other materials provided with the distribution.
 *
 * - Neither the names of the Roberts Group, Johns Hopkins University,
 * nor the names of its contributors may be used to endorse or
 * promote products derived from this Software without specific prior written
 * permission.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
 * THE CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR
 * OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS WITH THE SOFTWARE.
 *
 * Author(s): Elijah Roberts, Max Klein
 */
#include <cmath>
#include <limits>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "lm/Print.h"
#include "lm/Types.h"
#include "lm/input/Input.pb.h"
#include "lm/io/OrderParametersValues.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/io/SpeciesCounts.pb.h"
#include "lm/io/TrajectoryState.pb.h"
#include "lm/tiling/Tilings.h"
#include "lm/trajectory/Trajectory.h"

using std::list;
using std::map;
using std::string;
using std::vector;

using lm::input::DiffusionModel;
using lm::input::ReactionModel;
using lm::io::TrajectoryState;

namespace lm {
namespace trajectory {

const std::string Trajectory::status_strings[] = {"ABORTED",
                                                  "FINISHED",
                                                  "NOT_STARTED",
                                                  "RUNNING",
                                                  "WAITING"};

Trajectory::Trajectory(uint64_t id)
:id(id),workUnitsPerformed(0),state(),limitReached(),status(NOT_STARTED)
{
    state.set_trajectory_id(id);
    state.set_trajectory_started(false);
}

Trajectory::~Trajectory()
{
}

uint64_t Trajectory::getID() const
{
    return id;
}

Trajectory::Status Trajectory::getStatus() const
{
    return status;
}

const lm::types::TrajectoryLimit& Trajectory::getLimitReached() const
{
    return limitReached;
}

const lm::io::TrajectoryState& Trajectory::getState() const
{
    return state;
}

lm::io::TrajectoryState* Trajectory::getMutableState()
{
    return &state;
}

int64_t Trajectory::getWorkUnitsPerformed() const
{
    return workUnitsPerformed;
}

void Trajectory::setStatus(Status newStatus)
{
    status = newStatus;
}

void Trajectory::setLimitReached(const lm::types::TrajectoryLimit& newLimitReached)
{
    limitReached.CopyFrom(newLimitReached);
}

void Trajectory::setState(const lm::io::TrajectoryState& newState)
{
    state.CopyFrom(newState);
}

void Trajectory::incrementWorkUnitsPerformed()
{
    workUnitsPerformed++;
}

}
}
