/*
 * University of Illinois Open Source License
 * Copyright 2008-2012 Luthey-Schulten Group,
 * Copyright 2012-2019 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
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
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
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
#ifndef LM_CME_CMESOLVER_H_
#define LM_CME_CMESOLVER_H_

#include <algorithm>
#include <cstdio>
#include <list>
#include <map>
#include <pthread.h>
#include <string>
#include <utility>
#include <vector>

#include "lm/Math.h"
#include "lm/Types.h"
#include "lm/cme/ReactionModel.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/io/OrderParameterFirstPassageTimes.pb.h"
#include "lm/io/ParameterValues.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/io/TrajectoryState.pb.h"
#include "lm/me/FPTDeque.h"
#include "lm/me/MESolver.h"
#include "lm/me/PropensityFunction.h"
#include "lm/me/TimeSeriesList.h"
#include "lm/me/TrajectoryHistogram.h"
#include "lm/message/WorkUnitOutput.pb.h"
#include "lm/message/WorkUnitStatus.pb.h"
#include "lm/oparam/OrderParameterFunction.h"
#include "lm/rng/RandomGenerator.h"
#include "lm/thread/Thread.h"
#include "lm/types/Tilings.pb.h"
#include "lm/types/TrajectoryBarrier.pb.h"
#include "lm/types/TrajectoryLimits.pb.h"

using std::list;
using std::map;
using std::pair;
using std::string;
using std::vector;
using lm::me::MESolver;
using lm::rng::RandomGenerator;

namespace lm {

namespace ioFirstOrderPropensity {
class ReactionModel;
}

namespace cme {

class CMESolver : public MESolver
{
public:
    CMESolver(RandomGenerator::Distributions neededDists);
    virtual ~CMESolver();
    virtual void setComputeResources(vector<int> cpus, vector<int> gpus);
    virtual bool needsDiffusionModel() {return false;}
    virtual void setDiffusionModel(const lm::input::DiffusionModel& dm) {}
    virtual void setLimits(const lm::types::TrajectoryLimits& limits);
    virtual void setBarriers(const lm::types::TrajectoryBarriers& barriers);
    virtual void setOrderParameters(const lm::types::OrderParameters& opsBuf);
    virtual void setTilings(const lm::types::Tilings& t);
    virtual void setOutputOptions(const lm::input::OutputOptions& outputOptions);
    virtual bool needsReactionModel() {return true;}
    virtual void setReactionModel(const lm::input::ReactionModel& rm);
    virtual void reset();
    virtual void getState(lm::io::TrajectoryState* state, uint trajectoryNumber=0);
    virtual void setState(const lm::io::TrajectoryState& state, uint trajectoryNumber=0);
    virtual lm::message::WorkUnitOutput* getOutput(uint trajectoryNumber=0);
    virtual lm::message::WorkUnitStatus::Status getStatus(uint trajectoryNumber=0);
    virtual const lm::types::TrajectoryLimit& getLimitReached(uint trajectoryNumber=0);

protected:

    //
    // Trajectory helper functions.
    //
    virtual void trajectoryStarted();
    virtual bool performTimeIncrement(double deltaT);
    virtual bool isTrajectoryOutsideTimeLimits();
    virtual void setTrajectoryToMaxTime();
    virtual void timeUpdated();
    virtual bool performReactionEvent(uint r);
    virtual void applyReflectingBarriers();
    virtual bool isTrajectoryOutsideStateLimits();
    virtual void reactionCountsUpdated(uint r);
    virtual void speciesCountsUpdated();
    virtual void trajectoryFinished();

public:
    // Output options.
    bool writeInitialTrajectoryState, writeFinalTrajectoryState;
    bool writeSpeciesTimeSeries, writeDegreeAdvancementTimeSeries, writeOrderParameterTimeSeries;
    double degreeAdvancementWriteInterval, orderParameterWriteInterval, speciesWriteInterval;

protected:
    RandomGenerator::Distributions neededDists;
    RandomGenerator * rng;
    ReactionModel* reactionModel;
    bool hasUpdateReactionCountsListeners;
    bool hasUpdateSpeciesCountsListeners;

    // Order parameter function.
    int32_t numberOrderParameters;
    lm::oparam::OrderParameterFunction** orderParameterFunctions;

    // Trajectory output.
    lm::message::WorkUnitOutput* output;

    // Trajectory status.
    lm::message::WorkUnitStatus::Status status;
    uint64_t trajectoryId;
    bool previouslyStarted;

    // Limits for the trajectory.
    lm::types::TrajectoryLimits timeLimits;
    lm::types::TrajectoryLimits stateLimits;
    lm::types::TrajectoryLimit limitReached;

    // Barriers for the trajectory.
    lm::types::TrajectoryBarriers reflectingBarriers;
    lm::types::TrajectoryBarriers trackingBarriers;

    // First passage time variables.
    int numberFptSpecies, numberFptOrderParameters;
    lm::me::FPTDeque* fptSpeciesValues;
    lm::me::FPTDeque* fptOrderParameterValues;

    // The current state.
    uint64_t totalSteps;
    int32_t* speciesCounts;
    int32_t* speciesCountsPrevious;
    bool trackDegreeAdvancements;
    uint64_t* degreeAdvancements;
    double* orderParameterValues;
    double* orderParameterValuesPrevious;

    double time;
    double timeStep;    // stores last time step calculated, used for building histogram

    // Variables for recording time series data.
    lm::me::TimeSeriesList<int32_t> speciesTimeSeries;
    lm::me::TimeSeriesList<uint64_t> degreeAdvancementTimeSeries;
    lm::me::TimeSeriesList<double> orderParameterTimeSeries;

    size_t numberTrackingBarriers;
    uint64_t* trackingBarrierPriorCrossings;
    vector<int32_t>* trackingBarrierTimesCounts;
    vector<double>* trackingBarrierTimesTimes;
    vector<uint64_t>* trackingBarrierTimesSteps;

    // Variables for recording trajectory histograms.
    lm::types::Tilings tilings;
    double histogramBeginTime, histogramEndTime;
    vector<lm::me::TrajectoryHistogram<double>*> tilingHistograms;
};

}
}

#endif
