/*
 * University of Illinois Open Source License
 * Copyright 2012-2014 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Roberts Group
 *                  Johns Hopkins University
 *                  http://biophysics.jhu.edu/roberts/
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
#include <csignal>
#include <functional>
#include <list>
#include <map>
#include <numeric>
#include <string>
#include <vector>

#include "lm/Exceptions.h"
#include "lm/neus/NeusTrajectoryList.h"
#include "lm/io/CMEState.pb.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/io/SpeciesCounts.pb.h"
#include "lm/io/TrajectoryState.pb.h"
#include "lm/main/Main.h"
#include "lm/Print.h"
#include "lm/trajectory/Trajectory.h"
#include "lm/tiling/Tilings.h"

using lm::input::DiffusionModel;
using lm::input::ReactionModel;
using std::map;
using std::string;
using std::vector;

namespace lm {
namespace neus {

// these typedefs have to be here in .cpp because Direction is part of the FFPilotTrajectoryList definition
typedef map<lm::ffpilot::FFPilotTrajectoryList::Direction, CrossingsMap> CrossingsMapMap;
typedef map<lm::ffpilot::FFPilotTrajectoryList::Direction, DwellTimeMap> DwellTimeMapMap;
typedef map<lm::ffpilot::FFPilotTrajectoryList::Direction, FinishedTrajectoriesCountMap> FinishedTrajectoriesCountMapMap;

NeusTrajectoryList::NeusTrajectoryList(uint64_t simultaneousTrajectoryCount,lm::input::Input& input)
:TrajectoryList(input),
 crossingsPerPhase(atof(input.simulationParameters["crossingsPerPhase"].c_str())),
 direction(FORWARD),
 dwellTimes(),
 ffpilotPhase(0),
 finishedTrajectoriesCounts(),
 maxFFPilotPhase(input.tilings[0]->getEdgesCount()),
 maxPhaseZeroTime(atof(input.simulationParameters["maxPhaseZeroTime"].c_str())),
 simultaneousTrajectoryCount(simultaneousTrajectoryCount),
 xorShift(0,0)  //the rng object xorShift uses the current time as a seed when given 0,0 as constructor arguments
{
    init();
}

NeusTrajectoryList::~NeusTrajectoryList()
{
    // TODO: for the sake of this damn destructor, if for nothing else, I'm going to tear down the CrossingsMapMap stuff and replace it with something less obstinate
    // free all of the memory used by the savedCrossings member
    for (CrossingsMapMap::iterator mvit=savedCrossings.begin();mvit!=savedCrossings.end();++mvit)
    {
        for (CrossingsMap::iterator mit=mvit->second.begin();mit!=mvit->second.end();++mit)
        {
            for (CrossingVector::iterator vit=mit->second.begin();vit!=mit->second.end();++vit)
            {
                if (*vit!=NULL)
                {
                    delete *vit;
                    *vit=NULL;
                }
            }
        }
    }
}

void NeusTrajectoryList::init()
{
    initTrajectories(simultaneousTrajectoryCount);
    averageTilingHist.set_tiling_id(input.tilings[0]->getID());
    for (lm::tiling::EdgeIterator e_it=input.tilings[0]->begin();e_it!=input.tilings[0]->end();e_it++)
    {
        averageTilingHist.add_tile_vals(0);
    }
    dwellTimes[0] = 0;
    finishedTrajectoriesCounts[0] = 0;

}

void NeusTrajectoryList::initReversed() // TODO: need to verify that reactionModel has a reversed_initial_species_count field before running this method
{
    initTrajectories(simultaneousTrajectoryCount, true);
}

void NeusTrajectoryList::initTrajectories(uint64_t trajectoriesToStart,bool reversed)
{
    for (long long i=0; i<trajectoriesToStart; i++)
    {
        trajectories[trajectoryCount] = new lm::ffpilot::FFPilotTrajectory(trajectoryCount,ffpilotPhase,input,reversed);
        trajectoryCount++;
    }
}

void NeusTrajectoryList::initTrajectories(uint64_t trajectoriesToStart, lm::io::TrajectoryState* zerothTraj)
{
    for (long long i=0; i<trajectoriesToStart; i++)
    {
        trajectories[trajectoryCount] = new lm::ffpilot::FFPilotTrajectory(trajectoryCount,ffpilotPhase,input,zerothTraj);
        trajectoryCount++;
    }
}

void NeusTrajectoryList::initPhaseNTrajectories(uint64_t trajectoriesToStart)
{
    for (long long i=0; i<trajectoriesToStart; i++)
    {
        // Randomly choose a crossing state collected in the last round of ffpilot sampling, and use as the starting state for a new trajectory
        lm::io::TrajectoryState* randomCrossing = getRandomCrossing(ffpilotPhase - 1);
        initTrajectories(1, randomCrossing);
    }
}

lm::ffpilot::FFPilotTrajectory* NeusTrajectoryList::workUnitFinished(const lm::message::FinishedWorkUnit & finishedWorkUnitMsg)
{
    double simTime;
    // Call the base class method.
    lm::ffpilot::FFPilotTrajectory* traj = static_cast<lm::ffpilot::FFPilotTrajectory*>(TrajectoryList::workUnitFinished(finishedWorkUnitMsg));
//    Print::printf(Print::DEBUG, "finishedTrajectoryCount is: %d",finishedTrajectoriesCounts[ffpilotPhaseID]);
    // If the work unit stopped because it detected a crossing event...]
    if (finishedWorkUnitMsg.status()==lm::message::FinishedWorkUnit::LIMIT_REACHED)
    {
        // ...and if the crossing event was a forward flux pilot...
        if (traj->fluxedForward())
        {
            // ...add the work unit's final state to the appropriate list of crossings
            Print::printf(Print::INFO,"Crossing %d added to phase %d list", crossings[ffpilotPhase].size(), ffpilotPhase);
            addCrossing(finishedWorkUnitMsg);
        }
        // Regardless of whether this crossing was a forward or backwards flux, increment this phase's finished trajectories counter and dwell time, and delete the finished trajectory
        simTime = traj->getSimTime();
        dwellTimes[ffpilotPhase] += simTime;
        ++finishedTrajectoriesCounts[ffpilotPhase];
        deleteTrajectory(finishedWorkUnitMsg.final_state().trajectory_id());
        //Print::printf(Print::INFO, "ffpilotPhaseID: %d, crossings[ffpilot].size(): %d, finishedTrajectoriesCount %d, time: %f, oparam: %f", ffpilotPhaseID, crossings[ffpilotPhaseID].size(), finishedTrajectoriesCounts[ffpilotPhaseID], crossings[ffpilotPhaseID].back()->cme_state().species_counts().time(crossings[ffpilotPhaseID].back()->cme_state().species_counts().number_entries() - 1), calcTestCaseOParam(finishedWorkUnitMsg.final_state()));
        // If the forward flux pilot sampling is still in its 0th (ie initial) phase...
        if (isZerothPhase())
        {
            // ...and if enough time has passed for phase zero to be complete...
            if (isZerothPhaseDone(simTime))
            {
                if (crossings.find(0)==crossings.end()) Print::printf(Print::ERROR, "No crossings were recorded during forward flux pilot phase zero. Try increasing maxPhaseZeroTime");
                Print::printf(Print::INFO,"By the end of forward flux pilot phase zero, %d forward crossings were recorded", crossings[ffpilotPhase].size());
                // ...delete the currently running set of trajectories.
                deleteAllTrajectories();
                // Next, increment the ffpilot phase counter. If there are still more phases to run...
                ++ffpilotPhase;
                if (!isFFPilotDone())
                {
                    // ...start up a new set of trajectories and make room to store their data
                    dwellTimes[ffpilotPhase] = 0;
                    finishedTrajectoriesCounts[ffpilotPhase] = 0;
                    initPhaseNTrajectories(simultaneousTrajectoryCount);
                }
            }
            // ...otherwise we still have more time to go in phase zero...
            else
            {
                // ...so start one phase zero trajectory.
                initTrajectories(1, crossings[ffpilotPhase].back());
            }
        }
        // ...otherwise if ffpilotPhaseID > 0...
        else
        {
            // ...and if enough crossing events have been detected for this phase of forward flux pilot sampling...
            if (isPhaseDone())
            {
                // ...delete the currently running set of trajectories
                deleteAllTrajectories();
                // Next, increment the ffpilot phase counter. If there are still more phases to run...
                ++ffpilotPhase;
                if (!isFFPilotDone())
                {
                    // ...increment the interface position (by altering the increasing/decreasing limits)...
//                    ratchetInterfaces();
                    // ...and start up a new set of trajectories
                    dwellTimes[ffpilotPhase] = 0;
                    finishedTrajectoriesCounts[ffpilotPhase] = 0;
                    initPhaseNTrajectories(simultaneousTrajectoryCount);
                }
                // ...otherwise if the whole simulation is complete, output some data.
                else
                {
                    saveCrossings();
                    saveDwellTimes();
                    saveFinishedTrajectoriesCounts();
                    Print::printf(Print::INFO, "Phase 0 probability flux: %.10f", (double)crossings[0].size()/(maxPhaseZeroTime*simultaneousTrajectoryCount));
                    for (int i=1;i<maxFFPilotPhase;i++)
                    {
                        Print::printf(Print::INFO, "Crossing probability for interface at %f: %.10f", input.tilings[0]->getEdge(i), (double)crossings[i].size()/finishedTrajectoriesCounts[i]);
                    }
                    double Kab = (double)crossings[0].size()/(maxPhaseZeroTime*simultaneousTrajectoryCount);
                    for (int i=1;i<maxFFPilotPhase;i++)
                    {
                        Kab *= (double)crossings[i].size()/finishedTrajectoriesCounts[i];
                    }
                    Print::printf(Print::INFO, "Pseudo first order rate constant: %.10f", Kab);

                    // If we have to run ffpilot sampling in both directions, check if we're on the forward phase...
                    if (direction==FORWARD) // if (direction==FORWARD && bothDirections==TRUE)
                    {
                        // ...and if we are, reverse the arrangement of the edges and restart the simulation
                        reverse();
                        restart();
                    }
                    // If we're completely done with sampling in both directions, do the probability calculations
                    else
                    {
                        // This version of the probability calculation is taken from Dinner, 2010
                        double phaseZeroFluxA = (double)savedCrossings[FORWARD][0].size()/(maxPhaseZeroTime*simultaneousTrajectoryCount);
                        double phaseZeroFluxB = (double)savedCrossings[BACKWARD][0].size()/(maxPhaseZeroTime*simultaneousTrajectoryCount);
                        vector<double> paiaiplusone, pbibiplusone, pa0ai, pb0bi, pa0aiNormed, pb0biNormed;
                        for (int i=0;i<maxFFPilotPhase;i++)
                        {
                            paiaiplusone.push_back((double)savedCrossings[FORWARD][i].size()/savedFinishedTrajectoriesCounts[FORWARD][i]);
                            pbibiplusone.push_back((double)savedCrossings[BACKWARD][i].size()/savedFinishedTrajectoriesCounts[BACKWARD][i]);
                            pa0ai.push_back(paiaiplusone[1]);
                            pb0bi.push_back(pbibiplusone[1]);
                            for (int j=1;j<i;j++)
                            {
                                pa0ai[i]*=paiaiplusone[j];
                                pb0bi[i]*=pbibiplusone[j];
                            }
                        }
//                        std::partial_sum(paiaiplusone.begin(), paiaiplusone.end(), &pa0ai, std::multiplies<double>());
                        double Kab = pa0ai.back();
                        double Kba = pb0bi.back();
                        double Pa = Kba/(Kab + Kba);
                        double Pb = Kab/(Kab + Kba);
                        double totalWeight = 0;
//                        double totalProb = 0;
                        for (int i=1;i<maxFFPilotPhase;i++)
                        {
                            Print::printf(Print::INFO, "Forward tile %d probability: %.10f", i, Pa*phaseZeroFluxA*pa0ai[i]*savedDwellTimes[FORWARD][i]);
                            totalWeight+=Pa*phaseZeroFluxA*pa0ai[i]*savedDwellTimes[FORWARD][i];
                        }
                        for (int i=1;i<maxFFPilotPhase;i++)
                        {
                            Print::printf(Print::INFO, "Backward tile %d probability: %.10f", i, Pb*phaseZeroFluxB*pb0bi[i]*savedDwellTimes[BACKWARD][i]);
                            totalWeight+=Pb*phaseZeroFluxB*pb0bi[i]*savedDwellTimes[BACKWARD][i];
                        }
                        for (int i=0;i<maxFFPilotPhase;i++)
                        {
                            pa0aiNormed.push_back((Pa*phaseZeroFluxA*pa0ai[i]*savedDwellTimes[FORWARD][i])/totalWeight);
                            pb0biNormed.push_back((Pb*phaseZeroFluxB*pb0bi[i]*savedDwellTimes[BACKWARD][i])/totalWeight);
                        }
                        Print::printf(Print::INFO, "Pa: %.10f", Pa);
                        for (int i=1;i<maxFFPilotPhase;i++)
                        {
//                            Print::printf(Print::INFO, "Forward bit %d probability: %.10f", i, Pa*pa0aiNormed[i]);
//                            Print::printf(Print::INFO, "Backward bit %d probability: %.10f", maxFFPilotPhase-i, Pb*pb0biNormed[maxFFPilotPhase-i]);
                            Print::printf(Print::INFO, "Normalized tile %d probability: %.10f", i, pa0aiNormed[i]+pb0biNormed[maxFFPilotPhase-i]);
                        }
                        Print::printf(Print::INFO, "Pb: %.10f", Pb);
                    }
                }
            }
            // ...otherwise we still need to collect more crossing events for this phase of forward flux pilot sampling...
            else
            {
                // ...so start one phase N trajectory.
                initPhaseNTrajectories(1);
            }
        }
    }
    return traj;
        //        *run.mutable_initial_state() = trajectories->getTrajectoryState(nextTrajectory);
        //        Print::printf(Print::INFO, "Sending message to start work unit %d with trajectory %d on slot %d:%d.", run.work_unit_id(), nextTrajectory, workSlot->getSlotKey()[0], workSlot->getSlotKey()[1]);
        //        communicator.sendMessage(workSlot->getSlotKey()[0], workSlot->getSlotKey()[1], &msg);
        //        trajectories->updateTrajectoryStatus(nextTrajectory, FFPilotTrajectoryList::RUNNING);
}

// getters
CrossingVector NeusTrajectoryList::getCrossings(long long ffpilotPhase)
{
    return crossings[ffpilotPhase];
}

uint NeusTrajectoryList::getCrossingsPerPhase()
{
    return crossingsPerPhase;
}

long long NeusTrajectoryList::getFFPilotPhase()
{
    return ffpilotPhase;
}

double NeusTrajectoryList::getMaxPhaseZeroTime()
{
    return maxPhaseZeroTime;
}

lm::io::TrajectoryState* NeusTrajectoryList::getRandomCrossing(long long ffpilotPhase)
{
    unsigned i = floor(xorShift.getRandomDouble()*crossings[ffpilotPhase].size());
    return crossings[ffpilotPhase][i];
}


CrossingsMap NeusTrajectoryList::getSavedCrossings(lm::ffpilot::FFPilotTrajectoryList::Direction dir)
{
    return savedCrossings[dir];
}

// encapsulated inner loop functions
void NeusTrajectoryList::addCrossing(const lm::message::FinishedWorkUnit& finishedWorkUnitMsg)
{
    lm::io::TrajectoryState* newCrossing = new lm::io::TrajectoryState(finishedWorkUnitMsg.final_state());
    crossings[ffpilotPhase].push_back(newCrossing);
}

uint NeusTrajectoryList::incrFFPilotPhase()
{
    return ++ffpilotPhase;
}

bool NeusTrajectoryList::isFFPilotDone()
{
    return (ffpilotPhase>=maxFFPilotPhase);
}

bool NeusTrajectoryList::isPhaseDone()
{
    return (crossings[ffpilotPhase].size()>=crossingsPerPhase);
}

bool NeusTrajectoryList::isZerothPhase()
{
    return (ffpilotPhase==0);
}

bool NeusTrajectoryList::isZerothPhaseDone(double simTime)
{
    return simTime>=maxPhaseZeroTime;
}

void NeusTrajectoryList::reduceTilingHist(const lm::io::TilingHist& tHist)
{
}

void NeusTrajectoryList::restart()
{
    deleteAllTrajectories();
    crossings.clear();
    dwellTimes.clear(); dwellTimes[0] = 0;
    finishedTrajectoriesCounts.clear(); finishedTrajectoriesCounts[0] = 0;
    //std::fill(finishedTrajectoriesCounts.begin(), finishedTrajectoriesCounts.end(), 0);
    std::fill(averageTilingHist.mutable_tile_vals()->begin(), averageTilingHist.mutable_tile_vals()->end(), 0);
    ffpilotPhase = 0;
    this->initReversed();
}

void NeusTrajectoryList::reverse()
{
    direction = direction==FORWARD ? BACKWARD : FORWARD;
    input.tilings.reverse();
}

void NeusTrajectoryList::saveCrossings()
{

    savedCrossings.insert(CrossingsMapMap::value_type(direction, crossings));
}

void NeusTrajectoryList::saveDwellTimes()
{
    savedDwellTimes.insert(DwellTimeMapMap::value_type(direction, dwellTimes));
}

void NeusTrajectoryList::saveFinishedTrajectoriesCounts()
{
    savedFinishedTrajectoriesCounts.insert(FinishedTrajectoriesCountMapMap::value_type(direction, finishedTrajectoriesCounts));
}

}
}
