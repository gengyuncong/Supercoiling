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

#ifndef SLOT_SLOT_H
#define SLOT_SLOT_H

#include <string>

#include "lm/Types.h"
#include "lm/message/Endpoint.pb.h"
#include "lm/message/FinishedWorkUnit.pb.h"
#include "lm/resource/ComputeResources.h"

namespace lm {
namespace slot {

class SlotList;

class Slot
{
public:
    enum Status {NOT_STARTED, FREE, BUSY, DEAD};

public:
    Slot();
    Slot(int id, lm::resource::ComputeResources resources);
    ~Slot();

    int getId() const {return id;}
    uint getSimultaneousWorkUnits() const {return simultaneousWorkUnits;}
    void getStatsFromFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg);
    void resetStatistics();

protected:
    int id;
    Status status;
    lm::resource::ComputeResources resources;
    lm::message::Endpoint workUnitRunnerAddress;

	// the value in this attribute ultimately comes from solver->getSimultaneousTrajectories(). See WorkUnitRunner::run() for details (look for "msg->set_simultaneous_work_units(solver->getSimultaneousTrajectories());")
    uint simultaneousWorkUnits;

    long long stats_workUnitsME;
    long long stats_workUnitsStepsME;
    double stats_workUnitsTimeME;
    long long stats_workUnitsPDE;
    long long stats_workUnitsStepsPDE;
    double stats_workUnitsTimePDE;

    friend class SlotList;

};

}
}

#endif /* SLOT_SLOT_H */
