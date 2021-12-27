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

#ifdef OPT_AVX

#ifndef LM_AVX_GILLESPIEDSOLVERAVX_H_
#define LM_AVX_GILLESPIEDSOLVERAVX_H_

//#include <deque>
//#include <vector>

//#include <immintrin.h>

//#include "lm/ClassFactory.h"
//#include "lm/cme/GillespieDSolver.h"
//#include "lm/me/FPTDeque.h"
//#include "lm/message/WorkUnitOutput.pb.h"
//#include "lm/message/WorkUnitStatus.pb.h"
//#include "lm/rng/RandomGenerator.h"

//namespace lm {
//namespace avx {

//class GillespieDSolverAVX : public lm::cme::GillespieDSolver
//{
//public:
//    static bool registered;
//    static bool registerClass();
//    static void* allocateObject();

//public:
//    GillespieDSolverAVX();
//    virtual ~GillespieDSolverAVX();

//    virtual uint getSimultaneousTrajectories();
//    virtual void setLimits(const lm::types::TrajectoryLimits& limits);
//    virtual void setOrderParameters(const lm::types::OrderParameters& opsBuf);
//    virtual void setReactionModel(const lm::input::ReactionModel& rm);
//    virtual void setSparses();
//    virtual void reset();
//    virtual void getState(lm::io::TrajectoryState* state, uint trajectoryNumber=0);
//    virtual void setState(const lm::io::TrajectoryState& state, uint trajectoryNumber=0);
//    virtual uint64_t generateTrajectory(uint64_t maxSteps);
//    virtual lm::message::WorkUnitOutput* getOutput(uint trajectoryNumber=0);
//    virtual lm::message::WorkUnitStatus::Status getStatus(uint trajectoryNumber=0);

//protected:
//    template <typename Value, typename StorageValue>
//    inline double initWriteIntervalAVX(double interval, Value* valueArray, int valueSize, std::vector<StorageValue>* valueVector, std::vector<double>* timeVector, int i)
//    {
//        if (not previouslyStarted[i] and (writeInitialTrajectoryState or almostDivisible(((double*)&time)[i], interval)))
//        {
//            // if output of the initial state has been requested, do that
//            for (uint j=0; j<valueSize; j++) valueVector[i].push_back(lround(valueArray[j*DOUBLES_PER_AVX+i]));
//            timeVector[i].push_back(((double*)&time)[i]);
//        }
//        // set the next write time. If the next write time happens to be the current time, skip it (since it was already handled either just now or at the end of the previous work unit)
//        return roundNextMultiple(((double*)&time)[i], interval);
//    }

//    virtual void allocateRngBuffers();
//    virtual void deallocateRngBuffers();
//    virtual void updateAllPropensities();
//    void updateTimeAVX(avxd nextTimeStep, avxd nextTime);
//    void callUpdateTimeListenersAVX();
//    void updatePropensities(avxd time, uint* sourceReaction);
//    void performReactionEventAVX(uint* reactionsToPerform);
//    void callUpdateSpeciesCountsListenersAVX(uint* reactionsToPerform);
//    bool isTrajectoryOutsideLimitsAVX();
//    void copyOutputToBaseSolver(uint trajectoryNumber);
//    void copyTrajectoryStateToBaseSolver(uint trajectoryNumber);
//    void copyTrajectoryStateFromBaseSolver(uint trajectoryNumber);
//    void copyOutputFromBaseSolver(uint trajectoryNumber);

//protected:
//    // If the trajectory has been initialized.
//    bool initialized[DOUBLES_PER_AVX];

//    // Trajectory output.
//    lm::message::WorkUnitOutput* output[DOUBLES_PER_AVX];

//    // Trajectory status.
//    lm::message::WorkUnitStatus::Status status[DOUBLES_PER_AVX];

//    // Limits for the trajectory.
//    lm::types::TrajectoryLimit::LimitType limitTypeReached[DOUBLES_PER_AVX];
//    double* limitValues;

//    // First passage time variables.
//    int fptAllocatedValues;
//    lm::me::FPTDeque* fptValues;
//    double* fptMinValues;
//    double* fptMaxValues;

//    // First passage time order parameter variables.
//    uint numberFptOPValues;
//    double* fptOPMinValuesAchieved;
//    double* fptOPMaxValuesAchieved;
//    std::deque<double>* fptOPValues;
//    std::deque<double>* fptOPTimes;

//    // The current state.
//    uint64_t trajectoryId[DOUBLES_PER_AVX];
//    bool previouslyStarted[DOUBLES_PER_AVX];
//    double* speciesCounts;
//    double* propensities;
//    double* orderParameterValues;
//    double* orderParameterPreviousValues;
//    uint64_t* degreeAdvancements;

//    avxd timeLimit;
//    avxd time;
//    avxd timeStep;

//    species_sparse_type* speciesSparse[DOUBLES_PER_AVX];
//    bool speciesSparseOwn;
//};

//}
//}

#endif
#endif
