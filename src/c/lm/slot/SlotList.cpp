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

#include <climits>
#include <cstdio>
#include <iostream>
#include <pthread.h>
#include <sstream>
#include <string>
#include <vector>

#include "hrtime.h"
#include "lm/Exceptions.h"
#include "lm/Math.h"
#include "lm/Print.h"
#include "lm/Types.h"
#include "lm/message/Communicator.h"
#include "lm/resource/ComputeResources.h"
#include "lm/slot/Slot.h"
#include "lm/slot/SlotList.h"
#include "lm/thread/Thread.h"

using lm::resource::ComputeResources;
using lm::thread::PthreadException;
using lm::slot::Slot;
using std::string;
using std::vector;

namespace lm {
namespace slot {

int SlotList::nextSlotId = 0;

SlotList::SlotList()
:ownSlots(false)
{
}

SlotList::~SlotList()
{
    // If we own the slots, delete them.
    if (ownSlots)
    {
        for (map<int,Slot*>::const_iterator it=slotMap.begin(); it!=slotMap.end(); it++)
            delete it->second;
        slotMap.clear();
    }
}

void SlotList::createAllSlots(lm::message::Communicator* communicator, map<string,ComputeResources> & allResources, double cpusPerSlot, double gpusPerSlot, bool useCPUAffinity, string meSolver, string pdeSolver)
{
    // Mark that we own the slots, since we created them.
    ownSlots = true;

    // Create the slots.
    for (map<string,ComputeResources>::iterator it=allResources.begin(); it != allResources.end(); it++)
	{
        createHostSlots(communicator, it->second, cpusPerSlot, gpusPerSlot, useCPUAffinity, meSolver, pdeSolver);
	}
}

void SlotList::createHostSlots(lm::message::Communicator* communicator, ComputeResources resources, double cpusPerSlot, double gpusPerSlot, bool useCPUAffinity, string meSolver, string pdeSolver)
{
    // Make sure that a least one compute resource was being requested, otherwise we can create infinite slots.
    if (cpusPerSlot == 0.0 && gpusPerSlot == 0.0)
        throw Exception("No compute resources were requested for each work unit runner.");

    // Figure out the constraints on the number of slots.
    int cpuSlotsConstraint = (cpusPerSlot > 0.0)?(int(floor(double(resources.cpuCores.size())/cpusPerSlot))):(INT_MAX);
    int gpuSlotsConstraint = (gpusPerSlot > 0.0)?(int(floor(double(resources.gpuDevices.size())/gpusPerSlot))):(INT_MAX);

    // Loop while we have enough resources to create another slot.
    int i;
    for (i=0; i<min(cpuSlotsConstraint,gpuSlotsConstraint); i++)
    {
        // Create the resources object for the slot.
        ComputeResources slotResources;
        slotResources.hostname = resources.hostname;
        slotResources.controllerAddress = resources.controllerAddress;

        // Assign the cpu resources.
        if (cpusPerSlot >= 1.0)
        {
            // Assign sequential resources to the same slot.
            int cpuResourcesToAssign=int(floor(cpusPerSlot));
            for (int j=0; j<cpuResourcesToAssign; j++)
            {
                slotResources.cpuCores.push_back(resources.cpuCores[static_cast<size_t>(i*cpuResourcesToAssign+j)]);
            }
        }
        else if (cpusPerSlot > 0.0)
        {
            // Assign the same resource to multiple slots using a round-robin approach.
            slotResources.cpuCores.push_back(resources.cpuCores[static_cast<size_t>(i)%resources.cpuCores.size()]);
        }

        // Assign the gpu resources.
        if (gpusPerSlot >= 1.0)
        {
            // Assign sequential resources to the same slot.
            int gpuResourcesToAssign=int(floor(gpusPerSlot));
            for (int j=0; j<gpuResourcesToAssign; j++)
            {
                slotResources.gpuDevices.push_back(resources.gpuDevices[static_cast<size_t>(i*gpuResourcesToAssign+j)]);
            }
        }
        else if (gpusPerSlot > 0.0)
        {
            // Assign the same resource to multiple slots using a round-robin approach.
            slotResources.gpuDevices.push_back(resources.gpuDevices[static_cast<size_t>(i)%resources.gpuDevices.size()]);
        }

        // Create the message to start the workers.
        lm::message::Message msg;

        // Create the slot.
        createSlot(nextSlotId++, slotResources, useCPUAffinity, &msg, meSolver, pdeSolver);

        // Send the message to create all of the work units runners for this process.
        communicator->sendMessage(resources.controllerAddress, &msg);

    }
}

void SlotList::createSlot(int slotId, ComputeResources resources, bool useCPUAffinity, lm::message::Message* msg, string meSolver, string pdeSolver)
{
    Print::printf(Print::INFO, "Creating slot %d on %s using resources: %s.", slotId, lm::message::Communicator::printableAddress(resources.controllerAddress).c_str(), resources.toString().c_str());

    // Add the slots to our list.
    slotMap[slotId] = new Slot(slotId, resources);

    // Add the start work unit runner message for the slot.
    lm::message::StartWorkUnitRunner* s = msg->add_start_work_unit_runner();
    s->set_work_unit_runner_id(slotId);
    s->set_use_cpu_affinity(useCPUAffinity);
    for (vector<int>::iterator it=resources.cpuCores.begin(); it != resources.cpuCores.end(); it++)
        s->add_cpu(*it);
    for (vector<int>::iterator it=resources.gpuDevices.begin(); it != resources.gpuDevices.end(); it++)
        s->add_gpu(*it);
    s->set_me_solver(meSolver);
    s->set_diffusion_pde_solver(pdeSolver);
}

bool SlotList::isManagingSlot(int slotId) const
{
    return (slotMap.count(slotId) != 0);
}

bool SlotList::isRunningWorkUnit(uint64_t workUnitId) const
{
    return (workUnitToSlotMap.count(workUnitId) != 0);
}

void SlotList::markSlotStarted(const lm::message::StartedWorkUnitRunner& msg)
{
    // Mark the work unit runner as started.
    if (slotMap.count(msg.work_unit_runner_id()) == 0) throw Exception("Invalid work unit runner id received in started work unit runner message", msg.work_unit_runner_id());
    if (slotMap[msg.work_unit_runner_id()]->status != Slot::NOT_STARTED) throw Exception("Work unit runner was previosuly started", msg.work_unit_runner_id());
    slotMap[msg.work_unit_runner_id()]->status = Slot::FREE;
    slotMap[msg.work_unit_runner_id()]->workUnitRunnerAddress = msg.address();
    slotMap[msg.work_unit_runner_id()]->simultaneousWorkUnits = msg.simultaneous_work_units();
}

bool SlotList::hasUnstartedSlots()
{
    // Check all of the slots to see if anything is unstarted.
    for (map<int,Slot*>::const_iterator it=slotMap.begin(); it!=slotMap.end(); it++)
    {
        if (it->second->status == Slot::NOT_STARTED)
            return true;
    }
    return false;
}

bool SlotList::hasFreeSlots()
{
    // Check all of the slots to see if anything is free.
    for (map<int,Slot*>::const_iterator it=slotMap.begin(); it!=slotMap.end(); it++)
    {
        if (it->second->status == Slot::FREE)
            return true;
    }
    return false;
}

size_t SlotList::numberFreeSlots()
{
    // Check all of the slots to see if anything is free.
    size_t count=0;
    for (map<int,Slot*>::const_iterator it=slotMap.begin(); it!=slotMap.end(); it++)
    {
        if (it->second->status == Slot::FREE)
            count++;
    }
    return count;
}

const Slot& SlotList::getFreeSlot()
{
    // Find the first free slot.
    for (map<int,Slot*>::iterator it=slotMap.begin(); it!=slotMap.end(); it++)
    {
        if (it->second->status == Slot::FREE)
            return *it->second;
    }
    throw Exception("No free slots available.");
}

bool SlotList::hasBusySlots()
{
    // Check all of the slots to see if anything is free.
    for (map<int,Slot*>::const_iterator it=slotMap.begin(); it!=slotMap.end(); it++)
    {
        if (it->second->status == Slot::BUSY)
            return true;
    }
    return false;
}

void SlotList::runWorkUnit(lm::message::Communicator* communicator, lm::message::Message* runWorkUnitMsg)
{
    uint64_t workUnitId = runWorkUnitMsg->run_work_unit().work_unit_id();

    // Make sure we are not processing this work unit.
    if (workUnitToSlotMap.count(workUnitId) > 0)
        throw new Exception("The work unit is already assigned to a work unit runner",static_cast<int>(workUnitId),workUnitToSlotMap[workUnitId]);

    // Find a free slot.
    for (map<int,Slot*>::iterator it=slotMap.begin(); it!=slotMap.end(); it++)
    {
        if (it->second->status == Slot::FREE)
        {
            // Mark the slot as busy.
            it->second->status = Slot::BUSY;

            // Track the association between the work unit and the slot.
            workUnitToSlotMap[workUnitId] = it->first;

            // Send the message to start the work unit.
            communicator->sendMessage(it->second->workUnitRunnerAddress, runWorkUnitMsg);
            return;
        }
    }
    throw Exception("Could not find a free slot to run the work unit",static_cast<int>(workUnitId));
}

void SlotList::workUnitFinished(const lm::message::FinishedWorkUnit& msg)
{
    uint64_t workUnitId = msg.work_unit_id();

    // Make sure we were processing this work unit.
    if (workUnitToSlotMap.count(workUnitId) == 0) throw new Exception("The completed work unit was not found",static_cast<int>(workUnitId));

    // Get the slot that was running the work unit.
    int slotId = workUnitToSlotMap[workUnitId];

    // Make sure the slot was correctly marked as busy.
    if (slotMap[slotId]->status != Slot::BUSY) throw Exception("Work unit runner was not marked as busy while running work unit",slotId,static_cast<int>(workUnitId));

    // store some info for later use by printSlotsStatistics
    slotMap[slotId]->getStatsFromFinishedWorkUnit(msg);

    // Erase the work unit from our map.
    workUnitToSlotMap.erase(workUnitId);

    // Mark the slot as free.
    slotMap[slotId]->status = Slot::FREE;
}

vector<SlotList> SlotList::split(size_t numberGroups)
{
    // Make sure the number of groups is valid.
    if (numberGroups == 0 || numberGroups > slotMap.size()) THROW_EXCEPTION(lm::RuntimeException, "Invalid number of groups to split slots into: %d, %lld.", numberGroups, slotMap.size());

    // Figure out how many slots go into each of the first N-1 groups.
    size_t slotsPerGroup = slotMap.size()/numberGroups;

    // Create the groups.
    vector<SlotList> slotGroups;
    for (size_t i=0; i<numberGroups; i++)
        slotGroups.push_back(SlotList());

    // Go through the slots and distribute them into the groups.
    size_t groupIndex=0;
    for (map<int,Slot*>::const_iterator it=slotMap.begin(); it!=slotMap.end(); it++)
    {
        // Add the slot to the next group.
        slotGroups[groupIndex].addSlotReference(it->second);

        // If the group is full, move to the next group, unless this is the last group.
        if (slotGroups[groupIndex].getNumberSlots() >= slotsPerGroup && groupIndex < numberGroups-1)
            groupIndex++;
    }

    return slotGroups;
}

void SlotList::addSlotReference(Slot* slot)
{
    slotMap[slot->id] = slot;
}

void SlotList::printPerformanceStatistics()
{
    // Print some performance statistics.
    Print::printf(Print::INFO, "%5s %40s %7s %14s %15s %9s %9s %14s %15s %9s %9s", "Slot", "Resources", "Address", "ME_Work_Units", "Steps", "Time", "Steps/Sec", "PDE_Work_Units", "Steps", "Time", "Steps/Sec");
    Print::printf(Print::INFO, "------------------------------------------------------------------------------------------------------------------------------------------------------------");
    for (map<int,Slot*>::const_iterator it=slotMap.begin(); it!=slotMap.end(); it++)
    {
        Slot& s = *it->second;
        Print::printf(Print::INFO, "% 5d %40s %7s % 14lld % 15lld %9.3e %9.3e % 14lld % 15lld %9.3e %9.3e", s.id, s.resources.toString().c_str(), lm::message::Communicator::printableAddress(s.workUnitRunnerAddress).c_str(), s.stats_workUnitsME, s.stats_workUnitsStepsME, s.stats_workUnitsTimeME, double(s.stats_workUnitsStepsME)/s.stats_workUnitsTimeME, s.stats_workUnitsPDE, s.stats_workUnitsStepsPDE, s.stats_workUnitsTimePDE, double(s.stats_workUnitsStepsPDE)/s.stats_workUnitsTimePDE);
    }
}

}
}
