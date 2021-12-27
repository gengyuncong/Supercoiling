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
#include <set>
#include <string>

#include "hrtime.h"
#include "lm/Print.h"
#include "lm/Types.h"
#include "lm/input/Input.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/io/SpeciesCounts.pb.h"
#include "lm/io/TrajectoryState.pb.h"
#include "lm/main/Globals.h"
#include "lm/message/Message.pb.h"
#include "lm/replicates/ReplicateTrajectoryList.h"
#include "lm/trajectory/Trajectory.h"
#include "lm/trajectory/TrajectoryList.h"

#ifndef UINT64_MAX
#define UINT64_MAX        18446744073709551615ULL
#endif

using std::map;
using std::set;
using std::string;

using lm::trajectory::Trajectory;

namespace lm {
namespace replicates {

ReplicateTrajectoryList::ReplicateTrajectoryList()
{
}

ReplicateTrajectoryList::~ReplicateTrajectoryList()
{
}

void ReplicateTrajectoryList::processWorkUnitFinished(uint64_t workUnitId, const lm::message::FinishedWorkUnit& msg)
{
    // Call the base class method.
    TrajectoryList::processWorkUnitFinished(workUnitId, msg);

    // Print out a message for any trajectories that finished.
    if (replicatePrintMessages)
    {
        for (int i=0; i<msg.part_status_size(); i++)
        {
            uint64_t id = msg.part_status(i).final_state().trajectory_id();
            Trajectory* t = trajectories[id];
            if (t->getStatus() == Trajectory::FINISHED)
            {
                Print::printf(Print::INFO, "Replicate %lld completed with %8.2e of simulation time using %d work units.",
                              t->getID(), t->getState().cme_state().species_counts().time(0), t->getWorkUnitsPerformed());
            }
        }
    }
}

uint64_t ReplicateTrajectoryList::pickNextTrajectoryToRun(const set<uint64_t>& possibleTrajectories) const
{
    if (possibleTrajectories.size() <= 1000)
    {
        uint64_t minId=UINT64_MAX;
        double minTime=std::numeric_limits<double>::infinity();
        for (set<uint64_t>::const_iterator it=waitingTrajectories.begin(); it!=waitingTrajectories.end(); it++)
        {
            double time = trajectories.at(*it)->getState().cme_state().species_counts().time(0);
            if (time < minTime)
            {
                minTime = time;
                minId = *it;
            }
        }
        return minId;
    }
    else
    {
        return TrajectoryList::pickNextTrajectoryToRun(possibleTrajectories);
    }
}

void ReplicateTrajectoryList::printTrajectoryStatistics() const
{
    Print::printf(Print::INFO, "Replicate status: %d aborted, %d finished, %d running, %d waiting", abortedTrajectories.size(), finishedTrajectories.size(), runningTrajectories.size(), waitingTrajectories.size());
    if (replicatePrintMessages)
    {
        Print::printf(Print::INFO, "%10s %11s %9s %10s", "Id", "State", "Time", "Work_Units");
        Print::printf(Print::INFO, "-------------------------------------------");
        for (map<uint64_t,Trajectory*>::const_iterator it=trajectories.begin(); it!=trajectories.end(); it++)
        {
            uint64_t id = it->first;
            Trajectory* t = it->second;
            Print::printf(Print::INFO, "%10lld %11s %0.3e %10lld", id, Trajectory::status_strings[t->getStatus()].c_str(), t->getState().cme_state().species_counts().time(0), t->getWorkUnitsPerformed());
        }
    }
}

}
}
