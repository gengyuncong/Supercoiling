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

#ifndef LM_FFLUX_FFLUXTRAJECTORYLIST_H_
#define LM_FFLUX_FFLUXTRAJECTORYLIST_H_

//#include <google/protobuf/repeated_field.h>
//#include <map>
//#include <string>
//#include <vector>

//#include "lm/ffpilot/FFPilotTrajectory.h"
//#include "lm/ffpilot/FFPilotPhaseZeroTrajectory.h"
//#include "lm/ffpilot/input/FFPilotInput.h"
//#include "lm/ffpilot/input/FFPilotPhase.pb.h"
//#include "lm/ffpilot/io/FFPilotPhaseOutputWrap.h"
//#include "lm/input/Input.h"
//#include "lm/input/ReactionModel.pb.h"
//#include "lm/io/SpeciesTimeSeries.pb.h"
//#include "lm/io/TrajectoryState.pb.h"
//#include "lm/limit/LimitTrackingListWrap.h"
//#include "lm/message/Communicator.h"
//#include "lm/message/FinishedWorkUnit.pb.h"
//#include "lm/message/Message.pb.h"
//#include "lm/message/WorkUnitOutput.pb.h"
//#include "lm/message/WorkUnitStatus.pb.h"
//#include "lm/rng/XORShift.h"
//#include "lm/tiling/Tilings.h"
//#include "lm/trajectory/Trajectory.h"
#include "lm/trajectory/TrajectoryList.h"
//#include "lm/Types.h"

namespace lm {
namespace ffpilot {

class FFPilotTrajectoryList : public lm::trajectory::TrajectoryList
{
//public:
//    typedef lm::ffpilot::input::FFPilotInput FFPilotInput;
//    typedef lm::ffpilot::input::FFPilotPhase FFPilotPhaseMsg;
//    typedef lm::ffpilot::input::FFPilotPhaseLimit FFPilotPhaseLimitMsg;

//    // ffpilotPhase n==0 constructor
//    FFPilotTrajectoryList(uint64_t count, uint64_t simulationPhaseID, const FFPilotPhaseMsg& ffpilotPhase,
//                        uint simultaneousWorkUnits, const FFPilotInput& input,
//                        const lm::input::Basin& basin);

//    // ffpilotPhase n>0 constructor
//    FFPilotTrajectoryList(uint64_t count, uint64_t simulationPhaseID, const FFPilotPhaseMsg& ffpilotPhase,
//                        uint simultaneousWorkUnits, const FFPilotInput& input,
//                        const lm::protowrap::FFPilotPhaseOutputWrap& previousPhaseOutput);

//    // ffpilotPhase custom constructor
//    FFPilotTrajectoryList(uint64_t count, uint64_t simulationPhaseID, const FFPilotPhaseMsg& ffpilotPhase,
//                        uint simultaneousWorkUnits, const FFPilotInput& input);

//    virtual ~FFPilotTrajectoryList() {}

//    virtual void workUnitPartFinished(const lm::message::WorkUnitStatus& wusMsg, lm::trajectory::Trajectory* traj);

//    // accessors
//    virtual uint64_t getSimultaneousTrajectories();
//    virtual uint64_t getTrajectoriesToStart();

//    // mutators
//    virtual void initTrajectories(uint64_t trajectoriesToStart);
//    virtual void recycleFFPilotTrajectory(lm::trajectory::Trajectory* traj);

//protected:
//    // initializers
//    template <typename InputIterator> lm::trajectory::Trajectory* initFFPilotPhaseZeroTrajectory(const lm::input::Input& input, InputIterator speciesStart, InputIterator speciesEnd, double startTime, uint64_t phase, uint64_t id=DEFAULT_TRAJECTORY_ID)
//    {
//        return initTrajectory(new lm::ffpilot::FFPilotPhaseZeroTrajectory(input, speciesStart, speciesEnd, startTime, phase, resolveTrajectoryID(id)));
//    }
//    template <typename InputIterator> lm::trajectory::Trajectory* initFFPilotTrajectory(const lm::input::Input& input, InputIterator speciesStart, InputIterator speciesEnd, double startTime, uint64_t phase, uint64_t id=DEFAULT_TRAJECTORY_ID)
//    {
//        return initTrajectory(new lm::ffpilot::FFPilotTrajectory(input, speciesStart, speciesEnd, startTime, phase, resolveTrajectoryID(id)));
//    }
//    virtual void initTrajectoriesCyclic(uint64_t trajectoriesToStart);
//    virtual void initTrajectoriesUniformRandom(uint64_t trajectoriesToStart);
//    virtual lm::trajectory::Trajectory* recycleTrajectoryCyclic(uint64_t oldID);
//    virtual lm::trajectory::Trajectory* recycleTrajectoryUniformRandom(uint64_t oldID);

//protected:
//    const FFPilotInput& input;
//    const FFPilotPhaseMsg& ffpilotPhase;

//    // this is a pointer (and not a ref) because in some cases it has to be set to NULL
//    const lm::protowrap::FFPilotPhaseOutputWrap* previousPhaseOutputPtr;
//    lm::protowrap::FFPilotPhaseOutputWrap previousPhaseOutputCustomWrap;
//    lm::ffpilot::io::FFPilotPhaseOutput previousPhaseOutputCustom;

//    lm::limit::LimitTrackingListWrap limitTrackingListWrap;

//    size_t cyclicCounter;

//    uint simultaneousWorkUnits;
};

}
}
#endif
