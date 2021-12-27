/*
 * University of Illinois Open Source License
 * Copyright 2012-2014 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Roberts Group
 *               Johns Hopkins University
 *               http://biophysics.jhu.edu/roberts/
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
#ifndef LM_NEUS_NEUSTRAJECTORYLIST_H_
#define LM_NEUS_NEUSTRAJECTORYLIST_H_

#include <google/protobuf/repeated_field.h>
#include <map>
#include <string>
#include <vector>

#include "lm/neus/NeusTrajectory.h"
#include "lm/input/Input.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/io/TrajectoryState.pb.h"
#include "lm/message/FinishedWorkUnit.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/trajectory/TrajectoryList.h"
#include "lm/rng/XORShift.h"
#include "lm/tiling/Tilings.h"
#include "lm/Types.h"

namespace lm {
namespace neus {

typedef std::vector<lm::io::TrajectoryState*> CrossingVector;
typedef std::map<long long, CrossingVector> CrossingsMap;
typedef std::map<long long, double> DwellTimeMap;
typedef std::map<long long, long long> FinishedTrajectoriesCountMap;
typedef google::protobuf::RepeatedPtrField<lm::input::TrajectoryLimits::DecreasingOrderParameterLimit>::iterator decrLimitIterator;
typedef google::protobuf::RepeatedPtrField<lm::input::TrajectoryLimits::IncreasingOrderParameterLimit>::iterator incrLimitIterator;
typedef std::vector<lm::io::TilingHist*> TilingVector;

class NeusTrajectoryList : public lm::trajectory::TrajectoryList
{
public:
    // enumerated type used for describing the direction of the current ffpilot simulation relative to the arrangements (low-to-high or high-to-low) of the individual interfaces
    enum Direction {FORWARD, BACKWARD};

//    FFPilotTrajectoryList(uint64_t simultaneousTrajectoryCount,const lm::input::ReactionModel& reactionModel,const lm::input::DiffusionModel& diffusionModel,std::map<std::string,std::string>& simulationParameters, lm::tiling::Tilings& tilings);
    NeusTrajectoryList(uint64_t simultaneousTrajectoryCount,lm::input::Input& input);
    virtual ~NeusTrajectoryList();
    virtual void init();
    virtual void initReversed();
    virtual void initTrajectories(uint64_t toStartCount,bool reversed=false);
    virtual void initTrajectories(uint64_t toStartCount, lm::io::TrajectoryState* zerothTraj);
    virtual void initPhaseNTrajectories(uint64_t trajectoriesToStart);

    virtual lm::ffpilot::FFPilotTrajectory* workUnitFinished(const lm::message::FinishedWorkUnit & finishedWorkUnitMsg);

    // getters
    virtual CrossingVector getCrossings(long long ffpilotPhase);
    virtual uint getCrossingsPerPhase();
    virtual long long getFFPilotPhase();
    virtual double getMaxPhaseZeroTime();
    // Returns a randomly chosen crossing event (in the form of a TrajectoryState) collected durring forward flux pilot phase ffpilotPhase
    virtual lm::io::TrajectoryState* getRandomCrossing(long long ffpilotPhase);
    virtual CrossingsMap getSavedCrossings(lm::ffpilot::FFPilotTrajectoryList::Direction dir);

    // methods that encapsulate workUnitFinished inner loop tasks
    virtual void addCrossing(const lm::message::FinishedWorkUnit& finishedWorkUnitMsg);
    virtual uint incrFFPilotPhase();
    virtual bool isFFPilotDone();
    virtual bool isPhaseDone();
    virtual bool isZerothPhase();
    virtual bool isZerothPhaseDone(double);
    virtual void reduceTilingHist(const lm::io::TilingHist& tHist);
    virtual void restart();
    virtual void reverse();
    virtual void saveCrossings();
    virtual void saveDwellTimes();
    virtual void saveFinishedTrajectoriesCounts();

protected:
    Direction direction;
    lm::io::FFPilotOutput ffpilotOutput;
    long long ffpilotPhase;
    long long maxFFPilotPhase;
    uint64_t simultaneousTrajectoryCount;
    lm::rng::XORShift xorShift; //RNG used for randomly choosing a crossing in a crossing vector

    // members that hold the trajectory data used for the calculations at the end of ffpilot
    lm::io::TilingHist averageTilingHist;
    CrossingsMap crossings;
    DwellTimeMap dwellTimes;
    FinishedTrajectoriesCountMap finishedTrajectoriesCounts;
    std::map<lm::ffpilot::FFPilotTrajectoryList::Direction, CrossingsMap> savedCrossings;
    std::map<lm::ffpilot::FFPilotTrajectoryList::Direction, DwellTimeMap> savedDwellTimes;
    std::map<lm::ffpilot::FFPilotTrajectoryList::Direction, FinishedTrajectoriesCountMap> savedFinishedTrajectoriesCounts;
    std::map<lm::ffpilot::FFPilotTrajectoryList::Direction, lm::io::TilingHist*> savedHists;

    // user defined parameters that determine how the forward flux pilot sampling is carried out
    unsigned crossingsPerPhase; //the count of crossing events that should be collected for every ffpilot sampling phase
    double maxPhaseZeroTime;
};

}
}

#endif
