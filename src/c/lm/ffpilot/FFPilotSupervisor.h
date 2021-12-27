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

#ifndef LM_FFPILOT_FPILOTSUPERVISOR_H_
#define LM_FFPILOT_FPILOTSUPERVISOR_H_

#include <string>

////#include <google/protobuf/repeated_field.h>

////#include "hrtime.h"
#include "lm/Types.h"
#include "lm/ffpilot/FFPilotTrajectoryList.h"
#include "lm/ffpilot/FFPilotData.h"
#include "lm/input/OutputOptions.pb.h"
#include "lm/io/BarrierCrossingTimes.pb.h"
#include "lm/oparam/OrderParameterFunction.h"

////#include "lm/ffpilot/input/FFPilotPhase.pb.h"
////#include "lm/ffpilot/io/FFPilotPhaseOutput.pb.h"
////#include "lm/ffpilot/io/FFPilotPhaseOutputWrap.h"
////#include "lm/ffpilot/io/FFPilotStageOutput.pb.h"
////#include "lm/ffpilot/io/FFPilotStageOutputWrap.h"
////#include "lm/ffpilot/input/FFPilotInput.h"
////#include "lm/input/Options.pb.h"
////#include "lm/Iterator.h"
////#include "lm/message/Message.pb.h"
////#include "lm/protowrap/Repeated.h"
#include "lm/rng/RandomGenerator.h"
#include "lm/rng/XORShift.h"
#include "lm/simulation/StagedSimulationSupervisor.h"
#include "lm/types/TrajectoryLimits.pb.h"

namespace lm {
namespace ffpilot {

class FFPilotSupervisor : public lm::simulation::StagedSimulationSupervisor
{
public:
////    typedef lm::ffpilot::input::FFPilotStage FFPilotStageMsg;
////    typedef google::protobuf::RepeatedPtrField<FFPilotStageMsg> FFPilotStageRepeated;
////    typedef lm::ffpilot::input::FFPilotPhase FFPilotPhaseMsg;
////    typedef lm::protowrap::Repeated<FFPilotPhaseMsg> FFPilotPhasesWrap;
////    typedef lm::ffpilot::input::FFPilotPhaseLimit FFPilotPhaseLimitMsg;

////    typedef lm::protowrap::Repeated<lm::ffpilot::io::FFPilotStageOutput> FFPilotStageOutputsWrap;
////    typedef lm::protowrap::Repeated<lm::ffpilot::input::FFPilotPhaseLimit> FFPilotPhaseLimitsWrap;
////    typedef lm::protowrap::Repeated<lm::ffpilot::io::FFPilotPhaseOutput> FFPilotPhaseOutputsWrap;
////    typedef lm::protowrap::Repeated<lm::ffpilot::io::FFPilotPhaseOutputList> FFPilotPhaseOutputListsWrap;

    static bool registered;
    static bool registerClass();
    static void* allocateObject();

////    virtual int getRecvSleepMilliseconds();

public:
    FFPilotSupervisor();
    virtual ~FFPilotSupervisor();
    virtual std::string getClassName();
    virtual void printPerformanceStatistics();

protected:
    virtual void setInput(lm::input::Input* input);
    virtual int getNumberStages();
    virtual int getNumberPhases(int stage);

////    // setup methods that run once at the beginning of the simulation
    virtual void startSimulation();
////    virtual void initSimulationStageList();
////    virtual void initSimulationStageListCustom();
////    virtual FFPilotStageMsg* buildProductionStage(FFPilotStageMsg* productionStage, uint64_t replicateID, const lm::tiling::Tiling& tiling, int64_t basinIndex);
////    virtual void addTiling(FFPilotStageMsg* stage, const lm::tiling::Tiling& tiling, int64_t basinIndex);
////    virtual FFPilotStageMsg* buildPilotStage(FFPilotStageMsg* pilotStage, FFPilotStageMsg* productionStage);
////    virtual void addFFPilotPhases(FFPilotStageMsg* stage, FFPhaseEnums::TrajectoryGeneration trajGeneration, FFPhaseEnums::TrajectoryDuplication trajDuplication);
////    virtual void addOptions(FFPilotPhaseMsg* phase, const FFPilotStageMsg& stage);

////    // setup methods that run at the start of every ffpilot stage
    virtual void startStage(int stage);
    virtual void startPilotStage(int stage);
    virtual void startProductionStage(int stage);

////    virtual void addFFPilotStageOutput();
////    template <typename Value>
////    void* addFFPilotPhaseLimit(FFPilotPhaseMsg* ffpilotPhase, FFPhaseLimEnums::StopCondition stopCondition, Value value);
////    template <typename Value> void addFFPilotPhaseLimitsForPilotStage(FFPilotStageMsg* stage, FFPhaseLimEnums::StopCondition stopCondition, Value phaseZeroValue, Value value);
////    // TODO: spin this function off as part of an FFPilotPhase wrapper
////    static void buildFFPilotPhaseLimitEventsPerTrajectory(lm::ffpilot::input::FFPilotPhaseLimit* ffpilotPhaseLimit, const FFPilotPhaseMsg& ffpilotPhase, uint simultaneousWorkUnits);
////    static void buildFFPilotPhaseLimitTrajectoriesToRun(lm::ffpilot::input::FFPilotPhaseLimit* ffpilotPhaseLimit, const FFPilotPhaseMsg& ffpilotPhase, uint simultaneousWorkUnits);
////    void repeatFFPilotPhaseLimits(FFPilotStageMsg* stage, const lm::ffpilot::input::FFPilotPhaseLimit& limitToRepeat);
//////    template <typename ValueT> void repeatFFPilotPhaseLimits(lm::protowrap::Repeated<lm::ffpilot::input::FFPilotPhaseLimit>::iterator begin, const lm::ffpilot::input::FFPilotPhaseLimit& limitToRepeat);
////    virtual void addFFPilotPhaseLimitsFromInput(FFPilotStageMsg* productionStage);
////    virtual void addFFPilotPhaseLimitsFromStageOutput(FFPilotStageMsg* productionStage, const lm::protowrap::FFPilotStageOutputWrap& stageOutput);

////    // setup methods that run at the start of every ffpilot phase
    virtual void startPhase(int stage, int phase);
////    virtual void addFFPilotPhaseOutput();

    virtual void buildPhase(int stage, int phase);
    virtual void buildPhase0(int stage);
    virtual void buildPhase0FallbackTiling(int stageIndex);
    virtual void buildPhase0FallbackBasin(int stageIndex);
    virtual void buildPhaseN(int stage, int phase);
    virtual void createPhaseNTrajectory(int stageIndex, int phaseIndex);

    virtual bool assignWork(int stage, int phase);
    virtual void buildRunWorkUnitHeader(lm::message::RunWorkUnit* msg);
    virtual void buildRunWorkUnitParts(lm::message::RunWorkUnit* msg);
    virtual bool processEvent(lm::message::Message& msg);
    virtual void receivedFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg);
    virtual void processPhaseZeroBarrierCrossingTimesFallbackTiling(StageData& stage, PhaseData& phase, int64_t workUnitId, const lm::io::BarrierCrossingTimes& data0, const lm::io::BarrierCrossingTimes& data1);
    virtual void processPhaseZeroBarrierCrossingTimesFallbackBasin(StageData& stage, PhaseData& phase, int64_t workUnitId, const lm::io::BarrierCrossingTimes& data0);
    virtual void processPhaseNTrajectoryFinished(int stageIndex, int phaseIndex, const lm::io::TrajectoryState& state, bool success);


////    // methods that control what happens at the end of a ffpilotPhase
////    virtual bool _terminateSimulationPhase();
////    virtual void printFFPilotLimitProgress();
    virtual void finishPhase(int stage, int phase);
////    virtual void sendSimulationPhaseOutput();
////    virtual bool incrementSimulationPhase();

    // Methods that control what happens at the end of a stage.
    virtual void finishStage(int stage);
    virtual void finishPilotStage(int stage);
    virtual void savePilotStageOutput(int stageIndex, StageData& stageData, ndarray<uint64_t>& optimizedCounts);
    virtual void finishProductionStage(int stage);
    virtual void saveProductionStageOutput(int stageIndex, StageData& stageData, ndarray<double>& fpts);

    // Methods implementing the FFPilot optimizing equation.
    virtual ndarray<uint64_t> calculateOptimizedTrajectoryCounts(StageData& pilotStageData);
    virtual ndarray<uint64_t> calculateCountsMinimizeCost(double errorGoal, double errorGoalConfidence, const ndarray<double>& weights, const ndarray<double>& variances, const ndarray<double>& costs);
    virtual ndarray<uint64_t> calculateCountsMinimizeCounts(double errorGoal, double errorGoalConfidence, const ndarray<double>& weights, const ndarray<double>& variances);
    virtual ndarray<double> calculatePhaseFactors(const ndarray<double>& weights, const ndarray<double>& variances, const ndarray<double>& costs);

    // Ethods implementing production stage data analysis.
    virtual ndarray<double> calculateFirstPassageTimes(StageData& stageData);

////    virtual bool incrementSimulationStage();

////    // methods that run at the end of the entire simulation
    virtual void finishSimulation();

////    // methods that handle setting up RunWorkUnit messages
////    virtual const lm::input::Options& getOptions() {return currentPhase().options();}
////    virtual const lm::input::OutputOptions& getOutputOptions() {return currentPhase().output_options();}

////    // methods that handle FinishedWorkUnit messages
////    virtual void receivedFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg);
////    virtual void receivedFinishedWorkUnitPart(const lm::message::WorkUnitStatus& wusMsg);
////    virtual void receivedFinishedWorkUnitPartPhaseZero(const lm::message::WorkUnitStatus& wusMsg);

////    // accessors
////    virtual const FFPilotPhaseMsg& currentPhase() const {return *currentFFPilotPhaseIter;}
////    virtual int64_t currentFFPilotPhaseID() const {return currentPhase().phase_id();}
////    virtual std::string phaseInfo(bool path = false, const FFPilotPhaseMsg* phase = NULL, const FFPilotStageMsg* stage = NULL) const;
////    virtual const lm::ffpilot::input::FFPilotPhaseLimit& currentPhaseLimit() const {return currentPhase().ffpilot_phase_limit();}
////    virtual const lm::protowrap::FFPilotPhaseOutputWrap& currentPhaseOutput() const {return *currentFFPilotPhaseOutputWrapPtr;}
////    virtual int64_t finalFFPilotPhaseID() const {return currentStage().ffpilot_phases_size() - 1;}
////    virtual bool isCurrentPhaseLast() const {return isLast(currentFFPilotPhaseIter, currentStage().ffpilot_phases());} //{return currentStage().ffpilot_phases().end()==currentFFPilotPhaseIter;}
////    virtual const lm::protowrap::FFPilotPhaseOutputWrap& previousPhaseOutput() const {return *previousFFPilotPhaseOutputWrapPtr;}

////    virtual const FFPilotStageMsg& currentStage() const {return *currentFFPilotStageIter;}
////    virtual int64_t currentStageIndex() const {return _currentStageIndex;}
////    virtual std::string currentStageInfo(bool path=false, const FFPilotStageMsg* stage=NULL) const;
////    virtual const lm::protowrap::FFPilotStageOutputWrap& currentStageOutput() const {return currentFFPilotStageOutputWrap;}
////    virtual int getStageCount() const {return ffpilotStages.size();}
////    virtual bool isCurrentStageLast() const {return isLast(currentFFPilotStageIter, ffpilotStages);}  //{return currentFFPilotStageIter==ffpilotStageExecutionOrder.end();}

////    virtual const lm::tiling::Tiling& currentTiling() const {return currentTilingWrap;}

////    virtual std::string stageLogPilot(const lm::protowrap::FFPilotStageOutputWrap& stageOutput, double errorGoal, double errorGoalConfidence, vector<uint64_t>& trajectoryCounts) const;
////    virtual std::string stageLogProduction(const lm::protowrap::FFPilotStageOutputWrap& stageOutput) const;

////    // mutators
////    virtual FFPilotPhaseMsg* mutableCurrentPhase() {return &*currentFFPilotPhaseIter;}
////    virtual lm::ffpilot::input::FFPilotPhaseLimit* mutableCurrentPhaseLimit() {return mutableCurrentPhase()->mutable_ffpilot_phase_limit();}
////    virtual lm::protowrap::FFPilotPhaseOutputWrap* mutableCurrentPhaseOutput() {return currentFFPilotPhaseOutputWrapPtr;}

////    virtual FFPilotStageMsg* mutableCurrentStage() {return &*currentFFPilotStageIter;}
////    virtual lm::protowrap::FFPilotStageOutputWrap* mutableCurrentStageOutput() {return &currentFFPilotStageOutputWrap;}

////    virtual lm::tiling::Tiling* mutableCurrentTiling() {return &currentTilingWrap;}
////    virtual void setCurrentTiling(lm::input::Tiling* newCurrentTilingMsg) {currentTilingWrap.init(newCurrentTilingMsg, input->getOrderParameters());}

////    // setters/destructors to help with shadowing pointers in the base class
////    virtual void setInput(lm::input::Input* newInput);
////    virtual void destructInput() {lm::simulation::SimulationSupervisor::destructInput(); input = NULL;}
////    virtual void destructTrajectoryList() {lm::simulation::SimulationSupervisor::destructTrajectoryList(); trajectoryList = NULL;}

protected:
    lm::rng::RandomGenerator* random;
    int32_t orderParameterIndex;
    lm::oparam::OrderParameterFunction* orderParameterFunction;
    int32_t tilingIndex;
    ndarray<double>* tilingEdges;

    FFPilotData data;

    uint64_t nextTrajectoryId;
    FFPilotTrajectoryList* trajectoryList;
    lm::types::TrajectoryLimits trajectoryLimits;
    lm::types::TrajectoryBarriers trajectoryBarriers;
    lm::input::OutputOptions outputOptions;

////    FFPilotStageRepeated ffpilotStages;
////    FFPilotStageRepeated::iterator currentFFPilotStageIter;
////    FFPilotPhasesWrap::iterator currentFFPilotPhaseIter;

////    lm::tiling::Tiling currentTilingWrap;

////    // phase output messages
////    FFPilotPhaseOutputListsWrap::WrappedField ffpilotPhaseOutputListsMsg;
////    FFPilotPhaseOutputListsWrap ffpilotPhaseOutputListsWrap;
////    FFPilotPhaseOutputsWrap currentFFPilotPhaseOutputsWrap;

////    lm::message::Message ffpilotPhaseOutputContainingMsg;
////    lm::protowrap::FFPilotPhaseOutputWrap _ffpilotPhaseOutputWrap_0;
////    lm::protowrap::FFPilotPhaseOutputWrap _ffpilotPhaseOutputWrap_1;
////    lm::protowrap::FFPilotPhaseOutputWrap* previousFFPilotPhaseOutputWrapPtr;
////    lm::protowrap::FFPilotPhaseOutputWrap* currentFFPilotPhaseOutputWrapPtr;

////    // stage output messages
////    lm::message::Message ffpilotStageOutputRawContainingMsg;
////    lm::message::Message ffpilotStageOutputSummaryContainingMsg;
////    FFPilotStageOutputsWrap::WrappedField ffpilotStageOutputsMsg;
////    FFPilotStageOutputsWrap ffpilotStageOutputsWrap;
////    lm::protowrap::FFPilotStageOutputWrap currentFFPilotStageOutputWrap;

////    // flag that indicates that a phase has ended
////    bool ffpilotPhaseTerminated;

////    /*
////     * - flags that prevent output from being sent multiple times
////     *     - this is important when the supervisor has to cycle through the finishSimulationPhase() and finishSimulationStage multiple times, as at program end
////     */
////    bool simulationPhaseOutputSent;
////    bool simulationStageOutputSent;

////    // shadowing ptrs from the base class
////    lm::ffpilot::input::FFPilotInput* input;

};

}
}

#endif
