/*
 * Copyright 2012-2018 Johns Hopkins University
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

#include <cstdio>
#include <iomanip>
#include <map>
#include <sstream>
#include <string>

#include "lm/Print.h"
#include "lm/message/Communicator.h"
#include "lm/message/FinishedWorkUnit.pb.h"
#include "lm/slot/Slot.h"

using std::fixed;
using std::map;
using std::ostringstream;
using std::setprecision;
using std::scientific;
using std::string;
using std::setfill;
using std::setw;

namespace lm {
namespace slot {

Slot::Slot()
:id(-1),status(NOT_STARTED),resources(),simultaneousWorkUnits(0),
stats_workUnitsME(0),stats_workUnitsStepsME(0),stats_workUnitsTimeME(0.0),stats_workUnitsPDE(0),stats_workUnitsStepsPDE(0),stats_workUnitsTimePDE(0.0)
{
}

Slot::Slot(int id, lm::resource::ComputeResources resources)
:id(id),status(NOT_STARTED),resources(resources),simultaneousWorkUnits(0),
stats_workUnitsME(0),stats_workUnitsStepsME(0),stats_workUnitsTimeME(0.0),stats_workUnitsPDE(0),stats_workUnitsStepsPDE(0),stats_workUnitsTimePDE(0.0)
{
}

Slot::~Slot()
{
}

void Slot::getStatsFromFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg)
{
    if (msg.solver_type() == lm::types::SolverType::ME)
    {
        stats_workUnitsME++;
        stats_workUnitsStepsME += msg.steps();
        stats_workUnitsTimeME += msg.run_time();
    }
    else if (msg.solver_type() == lm::types::SolverType::DIFFUSION_PDE)
    {
        stats_workUnitsPDE++;
        stats_workUnitsStepsPDE += msg.steps();
        stats_workUnitsTimePDE += msg.run_time();
    }
}

void Slot::resetStatistics()
{
    stats_workUnitsME = 0;
    stats_workUnitsStepsME = 0;
    stats_workUnitsTimeME = 0.0;
    stats_workUnitsPDE = 0;
    stats_workUnitsStepsPDE = 0;
    stats_workUnitsTimePDE = 0.0;
}

}
}
