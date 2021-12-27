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

#ifndef LM_FFLUX_INPUT_FFLUXINPUT_H_
#define LM_FFLUX_INPUT_FFLUXINPUT_H_

//#include <google/protobuf/repeated_field.h>
//#include <list>
//#include <map>
//#include <string>
//#include <vector>

//#include "lm/ffpilot/input/FFPilotOptions.pb.h"
//#include "lm/ffpilot/input/FFPilotPhase.pb.h"
//#include "lm/ffpilot/input/FFPilotPhaseLimit.pb.h"
//#include "lm/ffpilot/input/FFPilotStage.pb.h"
//#include "lm/ffpilot/io/FFPilotStageOutput.pb.h"
//#include "lm/input/Input.h"
//#include "lm/io/hdf5/SimulationFile.h"
//#include "lm/types/BoundaryConditions.pb.h"
//#include "lm/input/DiffusionModel.pb.h"
//#include "lm/types/OrderParameters.pb.h"
//#include "lm/input/OutputOptions.pb.h"
//#include "lm/input/ReactionModel.pb.h"
//#include "lm/input/TrajectoryLimits.pb.h"
//#include "lm/limit/TrajectoryLimits.h"
//#include "lm/oparam/OParams.h"
//#include "lm/protowrap/Repeated.h"
//#include "lm/tiling/Tilings.h"
//#include "lm/Types.h"

//namespace lm {
//namespace ffpilot {
//namespace input {

//class FFPilotInput : public lm::input::Input
//{
//public:
//    static bool registered;
//    static bool registerClass();
//    static void* allocateObject();

//public:
//    typedef lm::ffpilot::input::FFPilotStage FFPilotStageMsg;
//    typedef google::protobuf::RepeatedPtrField<FFPilotStageMsg> FFPilotStageRepeated;
//    typedef lm::ffpilot::input::FFPilotPhase FFPilotPhaseMsg;
//    typedef lm::ffpilot::input::FFPilotPhaseLimit FFPilotPhaseLimitMsg;

////    FFPilotInput(const lm::io::hdf5::Hdf5File& file);
//    FFPilotInput();
//    virtual ~FFPilotInput() {};
//    virtual void init();

//// reinitializers
//    virtual void reinitOptions(int64_t phaseID);
//    virtual void reinitOutputOptions(const std::string& recordNamePrefix, bool isPilotStage);
//    virtual void reinitTrajectoryLimits(const FFPilotPhaseMsg& ffpilotPhase, const lm::tiling::Tiling& tiling);
//    virtual void reinitTrajectoryLimitsPhaseZero(const FFPilotPhaseMsg& ffpilotPhase, const lm::tiling::Tiling& tiling);

//// accessors
//    const lm::ffpilot::input::FFPilotOptions& getFFPilotOptionsMsg() const {return inputMsg.ffpilot_options();}
//    const FFPilotStageRepeated& getFFPilotStages() const {return inputMsg.ffpilot_stages();}

//    uint64_t phaseZeroBurnInCount() const {return getFFPilotOptionsMsg().phase_zero_burn_in_count();}
//    double phaseZeroSamplingMultiplier() const {return getFFPilotOptionsMsg().phase_zero_sampling_multiplier();}
//    uint64_t pilotStageCount() const {return getFFPilotOptionsMsg().pilot_stage_count();}
//    double errorGoal() const {return getFFPilotOptionsMsg().error_goal();}
//    double errorGoalConfidence() const {return getFFPilotOptionsMsg().error_goal_confidence();}
//    bool minimizeCost() const {return getFFPilotOptionsMsg().minimize_cost();}
//    uint64_t productionStageCountMinimum() const {return getFFPilotOptionsMsg().production_stage_count_minimum();}

//    bool hasErrorGoal() const {return getFFPilotOptionsMsg().has_error_goal();}
//    bool hasFFPilotOptions() const {return inputMsg.has_ffpilot_options();}
//    bool hasFFPilotStages() const {return not getFFPilotStages().empty();}

//    // TODO: make this private
//    virtual void initSanityCheck();
//protected:

//    virtual void initFFPilotOptions();
////    virtual void initFFPilotStages();
////    virtual void initFFPilotStagesCustom();
////    virtual void initFFPilotStagesFromScratch();

//};

//}
//}
//}

#endif /* LM_FFLUX_INPUT_FFLUXINPUT_H_ */
