/*
 * Copyright 2012-2019 Johns Hopkins University
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

#ifndef LM_SLOT_SLOTLIST_H
#define LM_SLOT_SLOTLIST_H

#include <map>
#include <string>
#include <vector>

#include "hrtime.h"
#include "lm/Types.h"
#include "lm/message/Communicator.h"
#include "lm/message/Message.pb.h"
#include "lm/resource/ComputeResources.h"
#include "lm/slot/Slot.h"
#include "lm/thread/Thread.h"

using lm::resource::ComputeResources;
using lm::slot::Slot;
using lm::thread::PthreadException;
using std::map;
using std::string;
using std::vector;

namespace lm {
namespace slot {

class SlotList
{
public:
    SlotList();
    ~SlotList();

    //create slot methods
    void createAllSlots(lm::message::Communicator* communicator, map<string,ComputeResources> & allResources, double cpusPerSlot, double gpusPerSlot, bool useCPUAffinity, string meSolver, string pdeSolver);
    bool isManagingSlot(int slotId) const;
    bool isRunningWorkUnit(uint64_t workUnitId) const;
    size_t getNumberSlots() const {return slotMap.size();}
    uint getSimultaneousWorkUnits() const {uint count=0; for (map<int,Slot*>::const_iterator it=slotMap.begin(); it!=slotMap.end(); it++) count+=it->second->simultaneousWorkUnits; return count;}
    void markSlotStarted(const lm::message::StartedWorkUnitRunner & msg);
    bool hasUnstartedSlots();
    bool hasFreeSlots();
    size_t numberFreeSlots();
    const Slot& getFreeSlot();
    bool hasBusySlots();
    void runWorkUnit(lm::message::Communicator* communicator, lm::message::Message* runWorkUnitMsg);
    void workUnitFinished(const lm::message::FinishedWorkUnit& msg);
    vector<SlotList> split(size_t groups);

    void printPerformanceStatistics();

protected:
    void createHostSlots(lm::message::Communicator* communicator, ComputeResources resources, double cpusPerSlot, double gpusPerSlot, bool useCPUAffinity, string meSolver, string pdeSolver);
    void createSlot(int slotId, ComputeResources resources, bool useCPUAffinity, lm::message::Message* msg, string meSolver, string pdeSolver);
    void addSlotReference(Slot* slot);

protected:
    static int nextSlotId;

private:
    bool ownSlots;
    map<int,Slot*> slotMap;
    map<uint64_t,int> workUnitToSlotMap;
};

}
}

#endif /* LM_SLOT_SLOTLIST_H */
