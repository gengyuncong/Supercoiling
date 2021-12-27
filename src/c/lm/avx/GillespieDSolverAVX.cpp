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

//#include <cmath>
//#include <cstdio>
//#include <cstdlib>
//#include <limits>
//#include <list>
//#include <map>
//#include <string>
//#include <vector>

//#include "lm/ClassFactory.h"
//#include "lm/Exceptions.h"
//#include "lm/Tune.h"
//#include "lm/Math.h"
//#include "lm/avx/AVXMath.h"
//#include "lm/Print.h"
//#include "lm/Types.h"
//#include "lm/avx/GillespieDSolverAVX.h"
//#include "lm/cme/CMESolver.h"
//#include "lm/cme/ReactionModel.h"
//#include "lm/io/DegreeAdvancementTimeSeries.pb.h"
//#include "lm/io/FirstPassageTimes.pb.h"
//#include "lm/io/SpeciesTimeSeries.pb.h"
//#include "lm/io/TrajectorySparse.pb.h"
//#include "lm/message/Message.pb.h"
//#include "lm/message/WorkUnitOutput.pb.h"
//#include "lm/message/WorkUnitStatus.pb.h"
//#include "lm/rng/RandomGenerator.h"
//#include "lm/rng/XORShift.h"
//#include "lm/thread/Thread.h"
//#include "lm/thread/Worker.h"
//#include "lptf/Profile.h"
//#include "lptf/ProfileCodes.h"
//#include "robertslab/pbuf/NDArraySerializer.h"


//using std::deque;
//using std::list;
//using std::map;
//using std::pair;
//using std::string;
//using std::vector;
//using lm::rng::RandomGenerator;

//namespace lm {
//namespace avx {

//bool GillespieDSolverAVX::registered=GillespieDSolverAVX::registerClass();

//bool GillespieDSolverAVX::registerClass()
//{
//    lm::ClassFactory::getInstance().registerClass("lm::me::MESolver","lm::avx::GillespieDSolverAVX",&GillespieDSolverAVX::allocateObject);
//    return true;
//}

//void* GillespieDSolverAVX::allocateObject()
//{
//    // Allocate aligned memory for the object since it contains avxd variables.
//    GillespieDSolverAVX* obj;
//    POSIX_EXCEPTION_CHECK(posix_memalign((void**)&obj, DOUBLES_PER_AVX*sizeof(double), sizeof(GillespieDSolverAVX)));

//    // Construct the object in place.
//    new (obj) GillespieDSolverAVX();

//    return obj;
//}

//GillespieDSolverAVX::GillespieDSolverAVX()
//:GillespieDSolver(),limitValues(NULL),fptAllocatedValues(0),fptValues(NULL),fptMinValues(NULL),fptMaxValues(NULL),
// numberFptOPValues(0),fptOPMinValuesAchieved(NULL),fptOPMaxValuesAchieved(NULL),fptOPValues(NULL),fptOPTimes(NULL),
// speciesCounts(NULL),propensities(NULL),orderParameterValues(NULL),orderParameterPreviousValues(NULL),
// degreeAdvancements(NULL),speciesSparseOwn(false)
//{
//    // Initialize any avx variables.
//    timeLimit = _mm256_set1_pd(std::numeric_limits<double>::infinity());
//    time = _mm256_set1_pd(0.0);
//    timeStep = _mm256_set1_pd(0.0);

//    // Initialize any array variables.
//    for (int i=0; i<DOUBLES_PER_AVX; i++)
//    {
//        initialized[i] = false;
//        output[i] = new lm::message::WorkUnitOutput();
//        status[i] = lm::message::WorkUnitStatus::NONE;
//        limitTypeReached[i] = lm::types::TrajectoryLimit::NONE;
//        trajectoryId[i] = std::numeric_limits<uint64_t>::max();
//        previouslyStarted[i] = false;
//        speciesSparse[i] = NULL;
//    }
//}

//GillespieDSolverAVX::~GillespieDSolverAVX()
//{
//    // Free any array memory.
//    for (int i=0; i<DOUBLES_PER_AVX; i++)
//    {
//        if (output[i] != NULL) delete output[i]; output[i] = NULL;
//    }

//    // Free any memory associated with the species first passage times.
//    fptAllocatedValues = 0;
//    if (fptValues != NULL) delete[] fptValues; fptValues = NULL;
//    if (fptMinValues != NULL) free(fptMinValues); fptMinValues = NULL;
//    if (fptMaxValues != NULL) free(fptMaxValues); fptMaxValues = NULL;

//    // Free any memory associated with the order parameter first passage times.
//    numberFptOPValues = 0;
//    if (fptOPMinValuesAchieved != NULL) free(fptOPMinValuesAchieved); fptOPMinValuesAchieved = NULL;
//    if (fptOPMaxValuesAchieved != NULL) free(fptOPMaxValuesAchieved); fptOPMaxValuesAchieved = NULL;
//    if (fptOPValues != NULL) delete[] fptOPValues; fptOPValues = NULL;
//    if (fptOPTimes != NULL) delete[] fptOPTimes; fptOPTimes = NULL;

//    // Free any other memory.
//    if (degreeAdvancements != NULL) free(degreeAdvancements); degreeAdvancements = NULL;
//    if (limitValues != NULL) free(limitValues); limitValues = NULL;
//    if (propensities != NULL) free(propensities); propensities = NULL;
//    if (orderParameterValues != NULL) free(orderParameterValues); orderParameterValues = NULL;
//    if (orderParameterPreviousValues != NULL) free(orderParameterPreviousValues); orderParameterPreviousValues = NULL;
//    if (speciesCounts != NULL) free(speciesCounts); speciesCounts = NULL;

//    // speciesSparse are allocated one-by-one, delete the same way
//    if (speciesSparseOwn) for (int i=0; i < DOUBLES_PER_AVX; i++ ) if (speciesSparse[i] != NULL) {delete speciesSparse[i]; speciesSparse[i] = NULL;}
//}

//void GillespieDSolverAVX::allocateRngBuffers()
//{
//    if (expRngValues == NULL || rngValues == NULL)
//    {
//        POSIX_EXCEPTION_CHECK(posix_memalign((void**)&rngValues, DOUBLES_PER_AVX*sizeof(double), TUNE_LOCAL_RNG_CACHE_SIZE*sizeof(double)));
//        POSIX_EXCEPTION_CHECK(posix_memalign((void**)&expRngValues, DOUBLES_PER_AVX*sizeof(double), TUNE_LOCAL_RNG_CACHE_SIZE*sizeof(double)));
//        nextRngValue = TUNE_LOCAL_RNG_CACHE_SIZE;
//    }
//}

//void GillespieDSolverAVX::deallocateRngBuffers()
//{
//    // Delete the rng caches.
//    if (expRngValues != NULL) free(expRngValues); expRngValues = NULL;
//    if (rngValues != NULL) free(rngValues); rngValues = NULL;
//    nextRngValue = 0;
//}

//uint GillespieDSolverAVX::getSimultaneousTrajectories()
//{
//    return DOUBLES_PER_AVX;
//}

//void GillespieDSolverAVX::setLimits(const lm::types::TrajectoryLimits& lm)
//{
//    GillespieDSolver::setLimits(lm);

//    // (Re)set the time limit.
//    timeLimit = _mm256_set1_pd(CMESolver::timeLimit);

//    // Free any previous limit values.
//    if (limitValues != NULL) free(limitValues); limitValues = NULL;

//    if (numberLimits > 0)
//    {
//        // Allocate space for the limit values.
//        POSIX_EXCEPTION_CHECK(posix_memalign((void**)&limitValues, DOUBLES_PER_AVX*sizeof(double), numberLimits*DOUBLES_PER_AVX*sizeof(double)));

//        // Copy the limit values into the avx buffer.
//        for (int i=0; i<(int)numberLimits; i++)
//        {
//            for (int j=0; j<DOUBLES_PER_AVX; j++)
//                if (limits[i].type == lm::types::TrajectoryLimit::SPECIES)
//                    limitValues[i*DOUBLES_PER_AVX+j] = double(limits[i].ivalue);
//                else if (limits[i].type == lm::types::TrajectoryLimit::DEGREE_ADVANCEMENT)
//                    limitValues[i*DOUBLES_PER_AVX+j] = double(limits[i].uvalue);
//                else
//                    limitValues[i*DOUBLES_PER_AVX+j] = limits[i].dvalue;
//        }
//    }
//}

//void GillespieDSolverAVX::setOrderParameters(const lm::types::OrderParameters& ops)
//{
//    GillespieDSolver::setOrderParameters(ops);

//    // Free any previous state.
//    if (orderParameterValues != NULL) free(orderParameterValues); orderParameterValues = NULL;
//    if (orderParameterPreviousValues != NULL) free(orderParameterPreviousValues); orderParameterPreviousValues = NULL;

//    // Allocate space for the order parameters.
//    POSIX_EXCEPTION_CHECK(posix_memalign((void**)&orderParameterValues, DOUBLES_PER_AVX*sizeof(double), numberOrderParameters*DOUBLES_PER_AVX*sizeof(double)));
//    POSIX_EXCEPTION_CHECK(posix_memalign((void**)&orderParameterPreviousValues, DOUBLES_PER_AVX*sizeof(double), numberOrderParameters*DOUBLES_PER_AVX*sizeof(double)));
//}

//void GillespieDSolverAVX::setReactionModel(const lm::input::ReactionModel& rm)
//{
//    GillespieDSolver::setReactionModel(rm);

//    // Free any previous state.
//    if (degreeAdvancements != NULL) free(degreeAdvancements); propensities = NULL;
//    if (propensities != NULL) free(propensities); propensities = NULL;
//    if (speciesCounts != NULL) free(speciesCounts); speciesCounts = NULL;

//    // Allocate degree advancements table. We don't yet know how many DAs we're going to be tracking, so make sure there's space for all of them
//    POSIX_EXCEPTION_CHECK(posix_memalign((void**)&degreeAdvancements, DOUBLES_PER_AVX*sizeof(double), reactionModel->numberReactions*DOUBLES_PER_AVX*sizeof(uint64_t)));

//    // Allocate reaction propensities table.
//    POSIX_EXCEPTION_CHECK(posix_memalign((void**)&propensities, DOUBLES_PER_AVX*sizeof(double), reactionModel->numberReactions*DOUBLES_PER_AVX*sizeof(double)));

//    // Allocate species counts table.
//    POSIX_EXCEPTION_CHECK(posix_memalign((void**)&speciesCounts, DOUBLES_PER_AVX*sizeof(double), reactionModel->numberSpecies*DOUBLES_PER_AVX*sizeof(double)));
//}

//void GillespieDSolverAVX::setSparses()
//{
//    // call the base class non-AVX method to ensure sparses are set correctly when we need to fall back
//    CMESolver::setSparses();

//    if (binSpecies)
//    {
//        if (workUnitCondenseOutput)
//        {
//            // if condensing output, have all trajectories bin into the single speciesSparse from the base class
//            for (int i=0; i < DOUBLES_PER_AVX; i++) if (speciesSparse[i] == NULL) speciesSparse[i] = CMESolver::speciesSparse;
//        }
//        else
//        {
//            // otherwise initialize any sparse histograms and helpers, if needed and they haven't already been taken care of
//            for (int i=0; i < DOUBLES_PER_AVX; i++) {if (speciesSparse[i] == NULL) speciesSparse[i] = new species_sparse_type(); speciesSparseOwn = true;}
//        }
//        hasUpdateTimeListeners = true;
//    }
//}

//void GillespieDSolverAVX::reset()
//{
//    GillespieDSolver::reset();

//    // Reset the species counts.
//    for (int i=0; i<(int)reactionModel->numberSpecies*DOUBLES_PER_AVX; i++)
//    {
//        speciesCounts[i] = 0.0;
//    }

//    // Reset the propensities.
//    for (int i=0; i<(int)reactionModel->numberReactions*DOUBLES_PER_AVX; i++)
//    {
//        propensities[i] = 0.0;
//    }

//    // Reset the degree advancements.
//    for (int i=0; i<(int)numberDegreeAdvancements*DOUBLES_PER_AVX; i++)
//    {
//        degreeAdvancements[i] = 0;
//    }

//    // Reset any array variables.
//    for (int i=0; i<DOUBLES_PER_AVX; i++)
//    {
//        initialized[i] = false;
//        if (output[i] != NULL) delete output[i];
//        output[i] = new lm::message::WorkUnitOutput();
//        status[i] = lm::message::WorkUnitStatus::NONE;
//        limitTypeReached[i] = lm::types::TrajectoryLimit::NONE;
//        trajectoryId[i] = std::numeric_limits<uint64_t>::max();
//        previouslyStarted[i] = false;
//    }

//    // Reset the time.
//    time = _mm256_set1_pd(0.0);
//    timeStep = _mm256_set1_pd(0.0);

//    // Reset the order parameters.
//    for (int i=0; i<numberOrderParameters*DOUBLES_PER_AVX; i++)
//    {
//        orderParameterValues[i] = 0.0;
//        orderParameterPreviousValues[i] = 0.0;
//    }

//    // Reset the fpt values.
//    fptAllocatedValues = 0;
//    if (fptValues != NULL) delete[] fptValues; fptValues = NULL;
//    if (fptMinValues != NULL) free(fptMinValues); fptMinValues = NULL;
//    if (fptMaxValues != NULL) free(fptMaxValues); fptMaxValues = NULL;

//    // Reset the fpt order parameter values
//    numberFptOPValues = 0;
//    if (fptOPMinValuesAchieved != NULL) free(fptOPMinValuesAchieved); fptOPMinValuesAchieved = NULL;
//    if (fptOPMaxValuesAchieved != NULL) free(fptOPMaxValuesAchieved); fptOPMaxValuesAchieved = NULL;
//    if (fptOPValues != NULL) delete[] fptOPValues; fptOPValues = NULL;
//    if (fptOPTimes != NULL) delete[] fptOPTimes; fptOPTimes = NULL;
//}

//void GillespieDSolverAVX::getState(lm::io::TrajectoryState* state, uint trajectoryNumber)
//{
//    if (trajectoryNumber >= getSimultaneousTrajectories()) throw lm::InvalidArgException("trajectoryNumber", "exceeded the maximum number of simultaneous trajectories",trajectoryNumber,getSimultaneousTrajectories());

//    // Setup the base class state with the specified trajextory.
//    copyTrajectoryStateToBaseSolver(trajectoryNumber);

//    // Run the base get state method last.
//    CMESolver::getState(state, trajectoryNumber);
//}

//void GillespieDSolverAVX::setState(const lm::io::TrajectoryState& state, uint trajectoryNumber)
//{
//    // Run the base set state first.
//    CMESolver::setState(state, trajectoryNumber);

//    // Allocate space for the fpt values, if necessary.
//    if (numberFptSpecies > 0)
//    {
//        // See if we have not yet allocated space.
//        if (fptAllocatedValues == 0)
//        {
//            fptValues = new lm::me::FPTDeque[numberFptSpecies*DOUBLES_PER_AVX];
//            fptAllocatedValues = numberFptSpecies;
//            POSIX_EXCEPTION_CHECK(posix_memalign((void**)&fptMinValues, DOUBLES_PER_AVX*sizeof(double), fptAllocatedValues*DOUBLES_PER_AVX*sizeof(double)));
//            POSIX_EXCEPTION_CHECK(posix_memalign((void**)&fptMaxValues, DOUBLES_PER_AVX*sizeof(double), fptAllocatedValues*DOUBLES_PER_AVX*sizeof(double)));
//        }
//        //Otherwise, make sure the sizes match.
//        else if (fptAllocatedValues != numberFptSpecies)
//        {
//            throw Exception("Mismatch between the number of fpt tracked values between trajectories",fptAllocatedValues,numberFptSpecies);
//        }
//    }

//    // Allocate space for the fpt order parameter values, if necessary.
//    if (numberFptTrackedOrderParameters > 0)
//    {
//        // See if we have not yet allocated space.
//        if (numberFptOPValues == 0)
//        {
//            numberFptOPValues = numberFptTrackedOrderParameters;
//            POSIX_EXCEPTION_CHECK(posix_memalign((void**)&fptOPMinValuesAchieved, DOUBLES_PER_AVX*sizeof(double), numberFptOPValues*DOUBLES_PER_AVX*sizeof(double)));
//            POSIX_EXCEPTION_CHECK(posix_memalign((void**)&fptOPMaxValuesAchieved, DOUBLES_PER_AVX*sizeof(double), numberFptOPValues*DOUBLES_PER_AVX*sizeof(double)));
//            fptOPValues = new deque<double>[numberFptOPValues*DOUBLES_PER_AVX];
//            fptOPTimes = new deque<double>[numberFptOPValues*DOUBLES_PER_AVX];
//        }
//            //Otherwise, make sure the sizes match.
//        else if (numberFptOPValues != numberFptTrackedOrderParameters)
//        {
//            throw Exception("Mismatch between the number of fpt tracked order parameter values between trajectories",numberFptOPValues,numberFptTrackedOrderParameters);
//        }
//    }

//    // Copy the state from the base class.
//    copyTrajectoryStateFromBaseSolver(trajectoryNumber);

//    // Mark that this trajectory was initialized.
//    initialized[trajectoryNumber] = true;
//}

//void GillespieDSolverAVX::copyTrajectoryStateToBaseSolver(uint trajectoryNumber)
//{
//    // Set the status.
//    CMESolver::status = status[trajectoryNumber];

//    // Set the limit reached.
//    CMESolver::limitTypeReached = limitTypeReached[trajectoryNumber];

//    // Set the trajectory id.
//    CMESolver::trajectoryId = trajectoryId[trajectoryNumber];

//    // Set the trajectory started flag.
//    CMESolver::previouslyStarted = previouslyStarted[trajectoryNumber];

//    // Set the species counts.
//    for (uint i=0; i<reactionModel->numberSpecies; i++)
//    {
//        CMESolver::speciesCounts[i] = lround(speciesCounts[i*DOUBLES_PER_AVX+trajectoryNumber]);
//    }

//    // Set the time.
//    CMESolver::time = ((double*)&time)[trajectoryNumber];
//    CMESolver::timeStep = ((double*)&timeStep)[trajectoryNumber];

//    // Set the order parameters.
//    for (int i=0; i<numberOrderParameters; i++)
//    {
//        CMESolver::orderParameterValues[i] = orderParameterValues[i*DOUBLES_PER_AVX+trajectoryNumber];
//        CMESolver::orderParameterPreviousValues[i] = orderParameterPreviousValues[i*DOUBLES_PER_AVX+trajectoryNumber];
//    }

//    // Set the first passage times.
//    for (int i=0; i<numberFptSpecies; i++)
//    {
//        CMESolver::fptValues[i] = fptValues[i*DOUBLES_PER_AVX+trajectoryNumber];
//    }

//    // Set the order parameter first passage times.
//    for (uint i=0; i<numberFptOPValues; i++)
//    {
//        fptTrackedOrderParameters[i].minValueAchieved = lround(fptOPMinValuesAchieved[i*DOUBLES_PER_AVX+trajectoryNumber]);
//        fptTrackedOrderParameters[i].maxValueAchieved = lround(fptOPMaxValuesAchieved[i*DOUBLES_PER_AVX+trajectoryNumber]);
//        fptTrackedOrderParameters[i].fptValues = fptOPValues[i*DOUBLES_PER_AVX+trajectoryNumber];
//        fptTrackedOrderParameters[i].fptTimes = fptOPTimes[i*DOUBLES_PER_AVX+trajectoryNumber];
//    }

//    // Set the degree advancements.
//    for (uint i=0; i<numberDegreeAdvancements; i++)
//    {
//        CMESolver::degreeAdvancements[i] = degreeAdvancements[i*DOUBLES_PER_AVX+trajectoryNumber];
//    }

//    // Set the histogram bin values.
//    //TODO: implement
//}

//void GillespieDSolverAVX::copyTrajectoryStateFromBaseSolver(uint trajectoryNumber)
//{
//    // Set the status.
//    status[trajectoryNumber] = CMESolver::status;

//    // Set the limit reached.
//    limitTypeReached[trajectoryNumber] = CMESolver::limitTypeReached;

//    // Set the trajectory id.
//    trajectoryId[trajectoryNumber] = CMESolver::trajectoryId;

//    // Set the trajectory started flag.
//    previouslyStarted[trajectoryNumber] = CMESolver::previouslyStarted;

//    // Set the species counts.
//    for (uint i=0; i<reactionModel->numberSpecies; i++)
//    {
//        speciesCounts[i*DOUBLES_PER_AVX+trajectoryNumber] = double(CMESolver::speciesCounts[i]);
//    }

//    // Set the time.
//    ((double*)&time)[trajectoryNumber] = CMESolver::time;
//    ((double*)&timeStep)[trajectoryNumber] = CMESolver::timeStep;

//    // Set the order parameters.
//    for (int i=0; i<numberOrderParameters; i++)
//    {
//        orderParameterValues[i*DOUBLES_PER_AVX+trajectoryNumber] = CMESolver::orderParameterValues[i];
//        orderParameterPreviousValues[i*DOUBLES_PER_AVX+trajectoryNumber] = CMESolver::orderParameterPreviousValues[i];
//    }

//    // Set the first passage times.
//    for (int i=0; i<numberFptSpecies; i++)
//    {
//        fptValues[i*DOUBLES_PER_AVX+trajectoryNumber] = CMESolver::fptValues[i];
//        fptMinValues[i*DOUBLES_PER_AVX+trajectoryNumber] = double(CMESolver::fptValues[i].minValue);
//        fptMaxValues[i*DOUBLES_PER_AVX+trajectoryNumber] = double(CMESolver::fptValues[i].maxValue);
//    }

//    // Set the degree advancements.
//    for (uint i=0; i<numberDegreeAdvancements; i++)
//    {
//        degreeAdvancements[i*DOUBLES_PER_AVX+trajectoryNumber] = CMESolver::degreeAdvancements[i];
//    }

//    // Set the order parameter first passage times.
//    for (uint i=0; i<numberFptOPValues; i++)
//    {
//        fptOPMinValuesAchieved[i*DOUBLES_PER_AVX+trajectoryNumber] = double(fptTrackedOrderParameters[i].minValueAchieved);
//        fptOPMaxValuesAchieved[i*DOUBLES_PER_AVX+trajectoryNumber] = double(fptTrackedOrderParameters[i].maxValueAchieved);
//        fptOPValues[i*DOUBLES_PER_AVX+trajectoryNumber] = fptTrackedOrderParameters[i].fptValues;
//        fptOPTimes[i*DOUBLES_PER_AVX+trajectoryNumber] = fptTrackedOrderParameters[i].fptTimes;
//    }

//    // Set the histogram bin values.
//    //TODO: implement
//}

//void GillespieDSolverAVX::copyOutputToBaseSolver(uint trajectoryNumber)
//{
//    CMESolver::output->CopyFrom(*output[trajectoryNumber]);
//}


//void GillespieDSolverAVX::copyOutputFromBaseSolver(uint trajectoryNumber)
//{
//    output[trajectoryNumber]->CopyFrom(*CMESolver::output);
//}

//lm::message::WorkUnitOutput* GillespieDSolverAVX::getOutput(uint trajectoryNumber)
//{
//    if (trajectoryNumber >= getSimultaneousTrajectories()) throw lm::InvalidArgException("trajectoryNumber", "exceeded the maximum number of simultaneous trajectories",trajectoryNumber,getSimultaneousTrajectories());

//    // Get the output pointer.
//    lm::message::WorkUnitOutput* ret = output[trajectoryNumber];

//    // Forget about the pointer, since the caller is now responsible for it.
//    output[trajectoryNumber] = NULL;

//    return ret;

//}

//lm::message::WorkUnitStatus::Status GillespieDSolverAVX::getStatus(uint trajectoryNumber)
//{
//    if (trajectoryNumber >= getSimultaneousTrajectories()) throw lm::InvalidArgException("trajectoryNumber", "exceeded the maximum number of simultaneous trajectories",trajectoryNumber,getSimultaneousTrajectories());
//    return status[trajectoryNumber];
//}

//uint64_t GillespieDSolverAVX::generateTrajectory(uint64_t maxSteps)
//{
//    // See how many trajectories were initialized.
//    uint numberInitialized=0;
//    for (int i=0; i<DOUBLES_PER_AVX; i++)
//        if (initialized[i])
//            numberInitialized++;

//    // If any trajectories were not initialized, run them all with the base Gillespie solver.
//    if (numberInitialized != DOUBLES_PER_AVX)
//    {
//        uint64_t steps=0;
//        for (int i=0; i<DOUBLES_PER_AVX; i++)
//        {
//            if (initialized[i])
//            {
//                Print::printf(Print::DEBUG, "GillespieDSolverAVX started without a full set of trajectories, running trajectory %llu with the GillespieDSolver.", trajectoryId[i]);
//                copyOutputToBaseSolver(i);
//                copyTrajectoryStateToBaseSolver(i);
//                GillespieDSolver::updateAllPropensities();
//                steps += GillespieDSolver::generateTrajectory(maxSteps);
//                copyTrajectoryStateFromBaseSolver(i);
//                copyOutputFromBaseSolver(i);
//            }
//        }

//        // Skip any necessary rng values to get us back into the proper alignment frame.
//        nextRngValue += nextRngValue%DOUBLES_PER_AVX;

//        return steps;
//    }

//    if (reactionModel == NULL) throw Exception("GillespieDSolverAVX did not have a reaction model.");
//    if (propensities == NULL) throw Exception("GillespieDSolverAVX state was not initialized.");

//    // Make sure we have propensity functions for every reaction.
//    for (uint i=0; i<reactionModel->numberReactions; i++)
//        if (reactionModel->propensityFunctions[i] == NULL)
//            throw Exception("A reaction did not have a valid propensity function",i);

//    // Create local copies of the data for efficiency.
//    const uint numberReactions = reactionModel->numberReactions;

//    // Set the propensities to their initial values.
//    updateAllPropensities();

//    // Initialize the total propensity.
//    avxd totalPropensity = _mm256_setzero_pd();
//    for (uint i=0; i<numberReactions; i++)
//    {
//        avxd propensity = _mm256_load_pd(&propensities[i*DOUBLES_PER_AVX]);
//        totalPropensity = _mm256_add_pd(totalPropensity,propensity);
//    }

//    avxd eps = _mm256_set1_pd(EPS);

//    // If we are writing degree advancement steps, create the data set and get the write interval.
//    avxd nextDegreeAdvancementWriteTime;
//    vector<uint64_t> degreeAdvancementTimeSeriesCounts[DOUBLES_PER_AVX];
//    vector<double> degreeAdvancementTimeSeriesTimes[DOUBLES_PER_AVX];
//    if (writeDegreeAdvancementTimeSeries)
//    {
//        for (int i=0; i<DOUBLES_PER_AVX; i++)
//        {
//            ((double*)&nextDegreeAdvancementWriteTime)[i] = initWriteIntervalAVX(degreeAdvancementWriteInterval, degreeAdvancements, numberDegreeAdvancements, degreeAdvancementTimeSeriesCounts, degreeAdvancementTimeSeriesTimes, i);
//        }
//    }

//    // If we are writing order parameter time steps, create the data set and get the write interval.
//    avxd nextOrderParameterWriteTime;
//    vector<double> orderParameterTimeSeriesCounts[DOUBLES_PER_AVX], orderParameterTimeSeriesTimes[DOUBLES_PER_AVX];
//    if (writeOrderParameterTimeSeries)
//    {
//        for (int i=0; i<DOUBLES_PER_AVX; i++)
//        {
//            ((double*)&nextOrderParameterWriteTime)[i] = initWriteIntervalAVX(orderParameterWriteInterval, orderParameterValues, numberOrderParameters, orderParameterTimeSeriesCounts, orderParameterTimeSeriesTimes, i);
//        }
//    }

//    // If we are writing species time steps, create the data set and get the write interval.
//    avxd nextSpeciesWriteTime;
//    vector<int32_t> speciesTimeSeriesCounts[DOUBLES_PER_AVX];
//    vector<double> speciesTimeSeriesTimes[DOUBLES_PER_AVX];
//    if (writeSpeciesTimeSeries)
//    {
//        for (int i=0; i<DOUBLES_PER_AVX; i++)
//        {
//            ((double*)&nextSpeciesWriteTime)[i] = initWriteIntervalAVX(speciesWriteInterval, speciesCounts, reactionModel->numberSpecies, speciesTimeSeriesCounts, speciesTimeSeriesTimes, i);
//        }
//    }

//    // Run the direct method.
//    Print::printf(Print::DEBUG, "Running Gillespie direct AVX simulation for %d steps with %d species, %d reactions, %d species limits\n", maxSteps, reactionModel->numberSpecies, reactionModel->numberReactions, numberLimits);
//    PROF_BEGIN(PROF_SIM_EXECUTE);
//    uint64_t steps=0;
//    int allFalse;
//    int trueMask;
//    avxd comp;
//    avxd nextTimeStep;
//    avxd nextTime;
//    while (true)
//    {
//        // See if we have finished the steps.
//        if (steps >= maxSteps)
//        {
//            for (int i=0; i<DOUBLES_PER_AVX; i++)
//                status[i] = lm::message::WorkUnitStatus::STEPS_FINISHED;
//            break;
//        }

//        // Increment the steps.
//        steps++;

//        // See if we need to update our rng caches.
//        if (nextRngValue >= TUNE_LOCAL_RNG_CACHE_SIZE)
//        {
//            rng->getRandomDoubles(rngValues,TUNE_LOCAL_RNG_CACHE_SIZE, true);
//            rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE, true);
//            nextRngValue=0;
//        }

//        // Get the random values for this iteration though the loop.
//        avxd randomValue = _mm256_load_pd(&rngValues[nextRngValue]);
//        avxd expRandomValue = _mm256_load_pd(&expRngValues[nextRngValue]);

//        // Go to the next rng pair.
//        nextRngValue+=DOUBLES_PER_AVX;

//        // Calculate the time to the next reaction.
//        nextTimeStep = _mm256_div_pd(expRandomValue, totalPropensity);
//        nextTime = _mm256_add_pd(time,nextTimeStep);

//         // If any new time is past the end time, we are done.
//        comp = _mm256_cmp_pd(nextTime, timeLimit, _CMP_GE_OQ);
//        allFalse = _mm256_testz_pd(comp,comp);
//        if (!allFalse)
//        {
//            // Get a bitmask of all values that were true.
//            trueMask = _mm256_movemask_pd(comp);

//            // Go through the mask.
//            for (int i=0; i<DOUBLES_PER_AVX; i++)
//            {
//                // If this element was true, set that the max time limit was reached.
//                if (trueMask&(1<<i))
//                {
//                    status[i] = lm::message::WorkUnitStatus::LIMIT_REACHED;
//                    limitTypeReached[i] = lm::types::TrajectoryLimit::TIME;
//                }
//                else
//                {
//                    // Otherwise set that we finished steps.
//                    status[i] = lm::message::WorkUnitStatus::STEPS_FINISHED;
//                }
//            }
//            break;
//        }

//        // Otherwise, set the time and timestep since no reactions were over the max time.
//        updateTimeAVX(nextTimeStep, nextTime);

//        // If we are writing degree advancements, write out any time steps before this event occurred.
//        if (writeDegreeAdvancementTimeSeries)
//        {
//            // Loop until we have finished writing out all elements.
//            while (true)
//            {
//                // See if any elements still need time steps written.
//                comp = _mm256_cmp_pd(nextDegreeAdvancementWriteTime, _mm256_add_pd(time, eps), _CMP_LE_OQ);
//                trueMask = _mm256_movemask_pd(comp);
//                if (!trueMask) break;

//                // Go through the mask.
//                for (int i=0; i<DOUBLES_PER_AVX; i++)
//                {
//                    // If this element was true, write its counts.
//                    if (trueMask&(1<<i))
//                    {
//                        // Record the degree advancement counts.
//                        for (uint j=0; j<numberDegreeAdvancements; j++) degreeAdvancementTimeSeriesCounts[i].push_back(degreeAdvancements[j*DOUBLES_PER_AVX+i]);
//                        degreeAdvancementTimeSeriesTimes[i].push_back(((double*)&nextDegreeAdvancementWriteTime)[i]);
//                        ((double*)&nextDegreeAdvancementWriteTime)[i] += degreeAdvancementWriteInterval;
//                    }
//                }
//            }
//        }

//        // If we are writing order parameters, write out any time steps before this event occurred.
//        if (writeOrderParameterTimeSeries)
//        {
//            // Loop until we have finished writing out all elements.
//            while (true)
//            {
//                // See if any elements still need time steps written.
//                comp = _mm256_cmp_pd(nextOrderParameterWriteTime, _mm256_add_pd(time, eps), _CMP_LE_OQ);
//                trueMask = _mm256_movemask_pd(comp);
//                if (!trueMask) break;

//                // Go through the mask.
//                for (int i=0; i<DOUBLES_PER_AVX; i++)
//                {
//                    // If this element was true, write its counts.
//                    if (trueMask&(1<<i))
//                    {
//                        // Record the order parameter values.
//                        for (uint j=0; j<numberOrderParameters; j++) orderParameterTimeSeriesCounts[i].push_back(orderParameterValues[j*DOUBLES_PER_AVX+i]);
//                        orderParameterTimeSeriesTimes[i].push_back(((double*)&nextOrderParameterWriteTime)[i]);
//                        ((double*)&nextOrderParameterWriteTime)[i] += orderParameterWriteInterval;
//                    }
//                }
//            }
//        }

//        // If we are writing species counts, write out any time steps before this event occurred.
//        if (writeSpeciesTimeSeries)
//        {
//            // Loop until we have finished writing out all elements.
//            while (true)
//            {
//                // See if any elements still need time steps written.
//                comp = _mm256_cmp_pd(nextSpeciesWriteTime, _mm256_add_pd(time, eps), _CMP_LE_OQ);
//                trueMask = _mm256_movemask_pd(comp);
//                if (!trueMask) break;

//                // Go through the mask.
//                for (int i=0; i<DOUBLES_PER_AVX; i++)
//                {
//                    // If this element was true, write its counts.
//                    if (trueMask&(1<<i))
//                    {
//                        // Record the species counts.
//                        for (uint j=0; j<reactionModel->numberSpecies; j++) speciesTimeSeriesCounts[i].push_back(lround(speciesCounts[j*DOUBLES_PER_AVX+i]));
//                        speciesTimeSeriesTimes[i].push_back(((double*)&nextSpeciesWriteTime)[i]);
//                        ((double*)&nextSpeciesWriteTime)[i] += speciesWriteInterval;
//                    }
//                }
//            }
//        }

//        // Calculate a random propensity to figure out the reaction.
//        avxd rngPropensity = _mm256_mul_pd(randomValue, totalPropensity);

////        {
////        double* res = (double*)&totalPropensity;
////        printf("%8.2f %8.2f %8.2f %8.2f\n", res[0], res[1], res[2], res[3]);
////        printf("RNG1: %8.2f %8.2f %8.2f %8.2f (%d)\n", rngValues[rngNext], rngValues[rngNext+1], rngValues[rngNext+2], rngValues[rngNext+3], rngNext);
////        res = (double*)&rngValue;
////        printf("RNG2: %8.2f %8.2f %8.2f %8.2f\n", res[0], res[1], res[2], res[3]);
////        }

//        // Figure out which reaction it was.
//        uint reactionsSelected = 0;
//        uint reactionsToPerform[DOUBLES_PER_AVX];
//        memset(reactionsToPerform, 0xFF, DOUBLES_PER_AVX*sizeof(uint));
//        for (uint r=0; r<numberReactions; r++)
//        {
//            // Load the propensities.
//            avxd propensity = _mm256_load_pd(&propensities[r*DOUBLES_PER_AVX]);

//            // Compare the random propensity with this entry.
//            comp = _mm256_cmp_pd(rngPropensity, propensity, _CMP_LT_OQ);
//            allFalse = _mm256_testz_pd(comp,comp);

////            printf("Reaction: %d\n",r);
////            res = (double*)&rngPropensity;
////            printf("%8.2f %8.2f %8.2f %8.2f\n", res[0], res[1], res[2], res[3]);
////            res = (double*)&propensity;
////            printf("%8.2f %8.2f %8.2f %8.2f\n", res[0], res[1], res[2], res[3]);
////            res = (double*)&comp;
////            printf("%8.2f %8.2f %8.2f %8.2f\n", res[0], res[1], res[2], res[3]);
////            printf("All false: %d\n",allFalse);

//            // If any were true, figure out which.
//            if (!allFalse)
//            {
//                // Get a bitmask of all values that were true.
//                trueMask = _mm256_movemask_pd(comp);

//                // Go through the mask.
//                for (int i=0; i<DOUBLES_PER_AVX; i++)
//                {
//                    // If this element was true, save the reaction and set the random propensity to inf.
//                    if (trueMask&(1<<i))
//                    {
//                        reactionsToPerform[i] = r;
//                        ((double*)&rngPropensity)[i] = std::numeric_limits<double>::infinity();
//                        reactionsSelected++;
//                    }
//                }

//                // If we found all of the reactions, we are done.
//                if (reactionsSelected == DOUBLES_PER_AVX)
//                    break;
//            }


//            // Subtract this entry from the random propensity and loop again.
//            rngPropensity = _mm256_sub_pd(rngPropensity, propensity);
//        }

//        // If any reactions were not set, they must be the last reaction.
//        if (reactionsSelected != DOUBLES_PER_AVX)
//        {
//            for (int i=0; i<DOUBLES_PER_AVX; i++)
//            {
//                if (reactionsToPerform[i] >= 0xFFFFFFFF)
//                    reactionsToPerform[i] = numberReactions-1;
//            }
//        }

//        //printf("Reaction to perform: %d %d %d %d\n", reactionsToPerform[0], reactionsToPerform[1], reactionsToPerform[2], reactionsToPerform[3]);

//        // Update the species counts.
//        performReactionEventAVX(reactionsToPerform);

//        // If we are outside of the limits, stop the trajectory.
//        if (numberLimits > 0 && isTrajectoryOutsideLimitsAVX()) break;

//        // Update the propensites given the reaction that occurred.
//        updatePropensities(time, reactionsToPerform);

//        // Recalculate the total propensity.
//        totalPropensity = _mm256_setzero_pd();
//        for (uint i=0; i<numberReactions; i++)
//        {
//            avxd propensity = _mm256_load_pd(&propensities[i*DOUBLES_PER_AVX]);
//            totalPropensity = _mm256_add_pd(totalPropensity,propensity);
//        }

//        // If any total propensities is zero, stop the trajectory.
//        avxd comp = _mm256_cmp_pd(totalPropensity, _mm256_setzero_pd(), _CMP_LE_OQ);
//        allFalse = _mm256_testz_pd(comp,comp);
//        if (!allFalse)
//        {
//            // Get a bitmask of all values that were true.
//            trueMask = _mm256_movemask_pd(comp);

//            // Go through the mask.
//            for (int i=0; i<DOUBLES_PER_AVX; i++)
//            {
//                // If this element was true, save the reaction and set the random propensity to inf.
//                if (trueMask&(1<<i))
//                {
//                    // If we have a time limit, say that we reached it.
//                    if (((double*)&timeLimit)[i] < std::numeric_limits<double>::infinity())
//                    {
//                        ((double*)&timeStep)[i] = ((double*)&timeLimit)[i]-((double*)&time)[i];
//                        ((double*)&time)[i] = ((double*)&timeLimit)[i];
//                        status[i] = lm::message::WorkUnitStatus::LIMIT_REACHED;
//                        limitTypeReached[i] = lm::types::TrajectoryLimit::TIME;
//                    }

//                    // Otherwise, zero propensity is an error.
//                    else
//                    {
//                        status[i] = lm::message::WorkUnitStatus::ERROR;
//                    }
//                }
//                else
//                {
//                    status[i] = lm::message::WorkUnitStatus::STEPS_FINISHED;
//                }
//            }
//            break;
//        }

//        //Print::printf(Print::VERBOSE_DEBUG, "Step %d: time=%e, count=%d, prop=%e, totprop=%e",steps,time,speciesCounts[0],propensities[0],totalPropensity);

////        double* p = (double*)&time;
////        printf("Time:             %8.2f %8.2f %8.2f %8.2f\n", p[0], p[1], p[2], p[3]);
////        p = (double*)speciesCounts;
////        printf("Species counts:   %8.2f %8.2f %8.2f %8.2f\n", p[0], p[1], p[2], p[3]);
////        //printf("Species counts:   %8.2f %8.2f %8.2f %8.2f\n", p[4], p[5], p[6], p[7]);
////        p = (double*)&totalPropensity;
////        printf("Total Propensity: %8.2f %8.2f %8.2f %8.2f\n", p[0], p[1], p[2], p[3]);
////        printf("-----------------\n");
//    }
//    PROF_END(PROF_SIM_EXECUTE);

////    double* p = (double*)&time;
////    printf("Final Time:             %8.2f %8.2f %8.2f %8.2f\n", p[0], p[1], p[2], p[3]);
////    p = (double*)speciesCounts;
////    printf("Final Species counts:   %8.2f %8.2f %8.2f %8.2f\n", p[0], p[1], p[2], p[3]);
////    //printf("Final Species counts:   %8.2f %8.2f %8.2f %8.2f\n", p[4], p[5], p[6], p[7]);
////    p = (double*)&totalPropensity;
////    printf("Final Total Propensity: %8.2f %8.2f %8.2f %8.2f\n", p[0], p[1], p[2], p[3]);

//    // Finalize each of the trajectories.
//    for (int i=0; i<DOUBLES_PER_AVX; i++)
//    {
//        // See if we finished all of the steps.
//        if (status[i] == lm::message::WorkUnitStatus::STEPS_FINISHED)
//        {
//            Print::printf(Print::DEBUG, "Generated trajectory %llu with %llu steps through time %e.", trajectoryId[i], steps, ((double*)&time)[i]);
//        }

//        // If we finished the total time, write out the remaining time steps.
//        else if (status[i] == lm::message::WorkUnitStatus::LIMIT_REACHED && limitTypeReached[i] == lm::types::TrajectoryLimit::TIME)
//        {
//            ((double*)&time)[i] = ((double*)&timeLimit)[i];
//            Print::printf(Print::DEBUG, "Generated trajectory %llu through time %e.", trajectoryId[i], ((double*)&time)[i]);

//            if (writeDegreeAdvancementTimeSeries)
//            {
//                while (((double*)&nextDegreeAdvancementWriteTime)[i] <= ((double*)&time)[i] + EPS)
//                {
//                    // Record the degree advancements.
//                    for (uint j=0; j<numberDegreeAdvancements; j++) degreeAdvancementTimeSeriesCounts[i].push_back(degreeAdvancements[j*DOUBLES_PER_AVX+i]);
//                    degreeAdvancementTimeSeriesTimes[i].push_back(((double*)&nextDegreeAdvancementWriteTime)[i]);
//                    ((double*)&nextDegreeAdvancementWriteTime)[i] += degreeAdvancementWriteInterval;
//                }
//            }

//            if (writeOrderParameterTimeSeries)
//            {
//                while (((double*)&nextOrderParameterWriteTime)[i] <= ((double*)&time)[i] + EPS)
//                {
//                    // Record the degree advancements.
//                    for (uint j=0; j<numberOrderParameters; j++) orderParameterTimeSeriesCounts[i].push_back(orderParameterValues[j*DOUBLES_PER_AVX+i]);
//                    orderParameterTimeSeriesTimes[i].push_back(((double*)&nextOrderParameterWriteTime)[i]);
//                    ((double*)&nextOrderParameterWriteTime)[i] += orderParameterWriteInterval;
//                }
//            }

//            if (writeSpeciesTimeSeries)
//            {
//                while (((double*)&nextSpeciesWriteTime)[i] <= ((double*)&time)[i] + EPS)
//                {
//                    // Record the species counts.
//                    for (uint j=0; j<reactionModel->numberSpecies; j++) speciesTimeSeriesCounts[i].push_back(lround(speciesCounts[j*DOUBLES_PER_AVX+i]));
//                    speciesTimeSeriesTimes[i].push_back(((double*)&nextSpeciesWriteTime)[i]);
//                    ((double*)&nextSpeciesWriteTime)[i] += speciesWriteInterval;
//                }
//            }
//        }

//        // If we hit a limit, write out the final trajectory state if requested
//        if (status[i] == lm::message::WorkUnitStatus::LIMIT_REACHED and writeFinalTrajectoryState)
//        {
//            // write out the last degree advancement count, or ensure that it already has been
//            if (writeDegreeAdvancementTimeSeries and (degreeAdvancementTimeSeriesTimes[i].empty() or degreeAdvancementTimeSeriesTimes[i].back() + EPS < ((double*)&time)[i]))
//            {
//                for (uint j=0; j<numberDegreeAdvancements; j++) degreeAdvancementTimeSeriesCounts[i].push_back(degreeAdvancements[j*DOUBLES_PER_AVX+i]);
//                degreeAdvancementTimeSeriesTimes[i].push_back(((double*)&time)[i]);
//            }

//            // write out the last order parameter value, or ensure that it already has been
//            if (writeOrderParameterTimeSeries and (orderParameterTimeSeriesTimes[i].empty() or orderParameterTimeSeriesTimes[i].back() + EPS < ((double*)&time)[i]))
//            {
//                for (uint j=0; j<numberOrderParameters; j++) orderParameterTimeSeriesCounts[i].push_back(orderParameterValues[j*DOUBLES_PER_AVX+i]);
//                orderParameterTimeSeriesTimes[i].push_back(((double*)&time)[i]);
//            }

//            // write out the last species count, or ensure that it already has been
//            if (writeSpeciesTimeSeries and (speciesTimeSeriesTimes[i].empty() or speciesTimeSeriesTimes[i].back() + EPS < ((double*)&time)[i]))
//            {
//                for (uint j=0; j<reactionModel->numberSpecies; j++) speciesTimeSeriesCounts[i].push_back(lround(speciesCounts[j*DOUBLES_PER_AVX+i]));
//                speciesTimeSeriesTimes[i].push_back(((double*)&time)[i]);
//            }
//        }

//        // If we have any degree advancement time series data, add them to the output message.
//        daTimeSeriesWrap.set_arrays_in_output_msg(output[i], degreeAdvancementTimeSeriesCounts[i], degreeAdvancementTimeSeriesTimes[i], trajectoryId[i], numberDegreeAdvancements);

//        // If we have any order parameter time series data, add them to the output message.
//        opTimeSeriesWrap.set_arrays_in_output_msg(output[i], orderParameterTimeSeriesCounts[i], orderParameterTimeSeriesTimes[i], trajectoryId[i], numberOrderParameters);

//        // If we have any order parameter time series data, add them to the output message.
//        //speciesTimeSeriesWrap.set_arrays_in_output_msg(output[i], speciesTimeSeriesCounts[i], speciesTimeSeriesTimes[i], trajectoryId[i], reactionModel->numberSpecies);

//        // If we have any species time series data, add them to the output message.
//        if (speciesTimeSeriesCounts[i].size() > 0 || speciesTimeSeriesTimes[i].size() > 0)
//        {
//            // Mark that the message does contain some data.
//            output[i]->set_has_output(true);

//            // Make sure the arrays are of a consistent size.
//            if (speciesTimeSeriesCounts[i].size() == speciesTimeSeriesTimes[i].size()*reactionModel->numberSpecies)
//            {
//                lm::io::SpeciesTimeSeries* speciesTimeSeriesDataSet = output[i]->mutable_species_time_series();
//                speciesTimeSeriesDataSet->set_trajectory_id(trajectoryId[i]);
//                robertslab::pbuf::NDArraySerializer::serializeInto<double>(speciesTimeSeriesDataSet->mutable_times(), speciesTimeSeriesTimes[i].data(), utuple(speciesTimeSeriesTimes[i].size()));
//                robertslab::pbuf::NDArraySerializer::serializeInto<int32_t>(speciesTimeSeriesDataSet->mutable_counts(), speciesTimeSeriesCounts[i].data(), utuple(speciesTimeSeriesTimes[i].size(),reactionModel->numberSpecies));
//            }
//            else
//            {
//                Print::printf(Print::ERROR, "Species time series counts and time mismatch %d,%d,%d", speciesTimeSeriesCounts[i].size(), reactionModel->numberSpecies, speciesTimeSeriesTimes[i].size());
//            }
//        }

//        // If we have any species time series data, add them to the output message.
//        if (speciesSparse[i] != NULL and not (speciesSparse[i]->empty() or workUnitCondenseOutput))
//        {
//            // Mark that the message does contain some data.
//            output[i]->set_has_output(true);

//            lm::io::TrajectorySparse* trajSparseMsg = output[i]->mutable_work_unit_output_generic()->add_trajectory_sparses();

//            trajSparseMsg->set_trajectory_id(trajectoryId[i]);
//            speciesSparse[i]->serializeTo(trajSparseMsg->mutable_cme_sparse());
//            speciesSparse[i]->clear();
//        }

//        // If the simulation reached a limit and we are tracking first passage times, add them to the output message.
//        if (status[i] == lm::message::WorkUnitStatus::LIMIT_REACHED && numberFptSpecies > 0)
//        {
//            // Mark that the message does contain some data.
//            output[i]->set_has_output(true);

//            for (int j=0; j<numberFptSpecies; j++)
//            {
//                fptValues[j*DOUBLES_PER_AVX+i].serializeInto(output[i]->add_first_passage_times());
//                Print::printf(Print::DEBUG, "avx gillespie writing fpt for %d %d,%d\n", i, (int)output[i]->first_passage_times(j).trajectory_id(), output[i]->first_passage_times(j).species()); fflush(stdout);
//            }
//        }

//        if (status[i] == lm::message::WorkUnitStatus::LIMIT_REACHED && numberFptTrackedOrderParameters > 0)
//        {
//            // Mark that the message does contain some data.
//            output[i]->set_has_output(true);

//            for (int j=0; j<numberFptOPValues; j++)
//            {
//                int dataIndex = j*DOUBLES_PER_AVX+i;
//                fptTrackedOrderParameters[j].serializeTo(output[i]->add_order_parameter_first_passage_times(), trajectoryId[i], fptOPValues[dataIndex], fptOPTimes[dataIndex]);
//            }
//        }
//    }

//    return steps*DOUBLES_PER_AVX;
//}

//void GillespieDSolverAVX::updateTimeAVX(avxd nextTimeStep, avxd nextTime)
//{
//    timeStep = nextTimeStep;
//    time = nextTime;

//    if (hasUpdateTimeListeners) callUpdateTimeListenersAVX();
//}

//void GillespieDSolverAVX::callUpdateTimeListenersAVX()
//{
//    if (binSpecies)
//    {
//        for (int i=0; i<DOUBLES_PER_AVX; i++)
//        {
//            // contruct an observation vector from the species counts
//            for (uint j=0; j<reactionModel->numberSpecies; j++) speciesObs->push_back(lround(speciesCounts[j*DOUBLES_PER_AVX+i]));

//            // bin the species counts, weighted by the timestep
//            speciesSparse[i]->addObs(*speciesObs, ((double*)&timeStep)[i]);

//            // clear the obs vector for reuse
//            speciesObs->clear();
//        }
//    }
//}

//void GillespieDSolverAVX::updateAllPropensities()
//{
//    // Update the propensities.
//    for (uint i=0; i<reactionModel->numberReactions; i++)
//    {
//        avxd propensity = reactionModel->propensityFunctions[i]->calculateAvx(time, speciesCounts, reactionModel->numberSpecies);
//        _mm256_store_pd(&propensities[i*DOUBLES_PER_AVX], propensity);
//    }
//}

//void GillespieDSolverAVX::updatePropensities(avxd time, uint* sourceReaction)
//{
//    updateAllPropensities();
//    // TODO: implement
////    // Update the propensities of the dependent reactions.
////    for (uint i=0; i<reactionModel->numberDependentReactions[sourceReaction]; i++)
////    {
////        uint r = reactionModel->dependentReactions[sourceReaction][i];
////        propensities[r] = reactionModel->propensityFunctions[i]->calculateAvx(time, speciesCounts, reactionModel->numberSpecies);
////    }
//}

//void GillespieDSolverAVX::performReactionEventAVX(uint* reactionsToPerform)
//{
//    // Update the counts according to the dependency tables.
//    for (uint j=0; j<DOUBLES_PER_AVX; j++)
//    {
//        uint r = reactionsToPerform[j];
//        for (int i=0; i<(int)reactionModel->numberDependentSpecies[r]; i++)
//        {
//            speciesCounts[(reactionModel->dependentSpecies[r][i])*DOUBLES_PER_AVX+j] += double(reactionModel->dependentSpeciesChange[r][i]);
//        }
//    }
//    if (hasUpdateSpeciesCountsListeners) callUpdateSpeciesCountsListenersAVX(reactionsToPerform);
//}

//void GillespieDSolverAVX::callUpdateSpeciesCountsListenersAVX(uint* reactionsToPerform)
//{
//    // Update the degree advancement, if we are tracking it.
//    if (numberDegreeAdvancements > 0)
//    {
//        for (uint i=0; i<DOUBLES_PER_AVX; i++)
//        {
//            uint r = reactionsToPerform[i];
//            degreeAdvancements[r*DOUBLES_PER_AVX+i]++;
//        }
//    }

//    // Update the first passage time tables.
//    for (int i=0; i<numberFptSpecies; i++)
//    {
//        uint speciesIndex = fptValues[i*DOUBLES_PER_AVX+0].species;
//        avxd counts = _mm256_load_pd(&speciesCounts[speciesIndex*DOUBLES_PER_AVX]);

//        // Check if any trajectory went below the previous min or above the previous max.
//        avxd comp = _mm256_cmp_pd(counts, _mm256_load_pd(&fptMinValues[i*DOUBLES_PER_AVX]), _CMP_LT_OQ);
//        int minAllFalse = _mm256_testz_pd(comp,comp);
//        comp = _mm256_cmp_pd(counts, _mm256_load_pd(&fptMaxValues[i*DOUBLES_PER_AVX]), _CMP_GT_OQ);
//        int maxAllFalse = _mm256_testz_pd(comp,comp);
//        if (!minAllFalse || !maxAllFalse)
//        {
//            // At least on trajectory was outside the range, so update all of them.
//            for (int j=0; j<DOUBLES_PER_AVX; j++)
//            {
//                // Update the fpt.
//                fptValues[i*DOUBLES_PER_AVX+j].insert(lround(((double*)&counts)[j]), ((double*)&time)[j]);

//                // Reload the min and and max.
//                fptMinValues[i*DOUBLES_PER_AVX+j] = double(fptValues[i*DOUBLES_PER_AVX+j].minValue);
//                fptMaxValues[i*DOUBLES_PER_AVX+j] = double(fptValues[i*DOUBLES_PER_AVX+j].maxValue);
//            }
//        }
//    }

//    // Update any order parameters.
//    if (numberOrderParameters > 0)
//    {
//        // Copy the old order parameters.
//        memcpy(orderParameterPreviousValues, orderParameterValues, numberOrderParameters*sizeof(double)*DOUBLES_PER_AVX);

//        // Update any order parameters.
//        for (int i=0; i<numberOrderParameters; i++)
//        {
//            avxd value = orderParameterFunctions[i]->calculateAvx(time, speciesCounts, reactionModel->numberSpecies);
//            _mm256_store_pd(&orderParameterValues[i*DOUBLES_PER_AVX], value);
//        }
//    }

//    // Update the order parameter first passage time tables.
//    for (int i=0; i<numberFptOPValues; i++)
//    {
//        int allFalse;
//        int trueMask;
//        uint opValIndex = fptTrackedOrderParameters[i].oparamID*DOUBLES_PER_AVX;
//        uint fptopIndex = i*DOUBLES_PER_AVX;

//        // rounding version
//        //avxd counts = _mm256_round_pd(_mm256_load_pd(&orderParameterValues[opValIndex]), _MM_FROUND_TO_ZERO);

//        // Check if we went below the previous min
//        AVX_COMP_MASK_ALLFALSE(_CMP_LT_OQ, &orderParameterValues[opValIndex], &fptOPMinValuesAchieved[fptopIndex], &trueMask, &allFalse)
//        if (!allFalse)
//        {
//            // Go through the mask.
//            for (int j=0; j<DOUBLES_PER_AVX; j++)
//            {
//                // If this element was true, update the fpt tables.
//                if (trueMask&(1<<j))
//                {
//                    fptOPMinValuesAchieved[fptopIndex + j] = orderParameterValues[opValIndex + j];
//                    fptOPValues[fptopIndex + j].push_front(fptOPMinValuesAchieved[fptopIndex + j]);
//                    fptOPTimes[fptopIndex + j].push_front(((double*)&time)[j]);
//                }
//            }
//        }

//        // Check if we went above the previous max
//        AVX_COMP_MASK_ALLFALSE(_CMP_GT_OQ, &orderParameterValues[opValIndex], &fptOPMaxValuesAchieved[fptopIndex], &trueMask, &allFalse)
//        if (!allFalse)
//        {
//            // Go through the mask.
//            for (int j=0; j<DOUBLES_PER_AVX; j++)
//            {
//                // If this element was true, update the fpt tables.
//                if (trueMask&(1<<j))
//                {
//                    fptOPMaxValuesAchieved[fptopIndex + j] = orderParameterValues[opValIndex + j];
//                    fptOPValues[fptopIndex + j].push_back(fptOPMaxValuesAchieved[fptopIndex + j]);
//                    fptOPTimes[fptopIndex + j].push_back(((double*)&time)[j]);
//                }
//            }
//        }
//    }

//    // Update any tilingHists.
//}

//bool GillespieDSolverAVX::isTrajectoryOutsideLimitsAVX()
//{
//    // Go through the limits.
//    for (uint i=0; i<numberLimits; i++)
//    {
//        int outsideLimitMask=0;
//        lm::limit::TrajectoryLimit& l = limits[i];
//        avxd limitValue = _mm256_load_pd(&limitValues[i*DOUBLES_PER_AVX]);
//        avxd comp1;
//        avxd comp2;

//        switch (l.type)
//        {
//        case TrajLimEnums::NONE: throw Exception("GillespieDSolverAVX tried to check a limit that did not have an associated LimitType"); break;
//        case TrajLimEnums::TIME: throw Exception("GillespieDSolverAVX reached a time limit that was mixed in with the other limits"); break;

//        case TrajLimEnums::SPECIES:
//            switch (l.stoppingCondition)
//            {
//            case TrajLimEnums::MIN:
//                if (l.includeEndpoint)
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&speciesCounts[l.valueID * DOUBLES_PER_AVX]), limitValue, _CMP_LE_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1);
//                }
//                else
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&speciesCounts[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_LT_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1);
//                }
//                break;
//            case TrajLimEnums::MAX:
//                if (l.includeEndpoint)
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&speciesCounts[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_GE_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1);
//                }
//                else
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&speciesCounts[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_GT_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1);
//                }
//                break;
//            case TrajLimEnums::INCREASING: throw Exception("unimplemented"); break;
//            case TrajLimEnums::DECREASING: throw Exception("unimplemented"); break;
//            }
//            break;

//        case TrajLimEnums::ORDER_PARAMETER:
//            switch (l.stoppingCondition)
//            {
//            case TrajLimEnums::MIN:
//                if (l.includeEndpoint)
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterValues[l.valueID * DOUBLES_PER_AVX]), limitValue, _CMP_LE_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1);
//                }
//                else
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterValues[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_LT_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1);
//                }
//                break;
//            case TrajLimEnums::MAX:
//                if (l.includeEndpoint)
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterValues[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_GE_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1);
//                }
//                else
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterValues[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_GT_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1);
//                }
//                break;
//            case TrajLimEnums::DECREASING:
//                if (l.includeEndpoint)
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterPreviousValues[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_GT_OQ);
//                    comp2 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterValues[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_LE_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1)&_mm256_movemask_pd(comp2);
//                }
//                else
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterPreviousValues[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_GE_OQ);
//                    comp2 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterValues[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_LT_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1)&_mm256_movemask_pd(comp2);
//                }
//                break;
//            case TrajLimEnums::INCREASING:
//                if (l.includeEndpoint)
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterPreviousValues[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_LT_OQ);
//                    comp2 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterValues[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_GE_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1)&_mm256_movemask_pd(comp2);
//                }
//                else
//                {
//                    comp1 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterPreviousValues[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_LE_OQ);
//                    comp2 = _mm256_cmp_pd(_mm256_load_pd(&orderParameterValues[l.valueID*DOUBLES_PER_AVX]), limitValue, _CMP_GT_OQ);
//                    outsideLimitMask = _mm256_movemask_pd(comp1)&_mm256_movemask_pd(comp2);
//                }
//                break;
//            }
//            break;

//        case TrajLimEnums::DEGREE_ADVANCEMENT:
//            switch (l.stoppingCondition)
//            {
//            case TrajLimEnums::MIN: throw Exception("unimplemented"); break;
//            case TrajLimEnums::MAX: throw Exception("unimplemented"); break;
//            case TrajLimEnums::INCREASING: throw Exception("unimplemented"); break;
//            case TrajLimEnums::DECREASING: throw Exception("unimplemented"); break;
//            }
//            break;
            
//        default:
//            break;
//        }

//        // Check if the limit was triggered.
//        if (outsideLimitMask)
//        {
//            // Go through the mask.
//            for (int j=0; j<DOUBLES_PER_AVX; j++)
//            {
//                // If this element was true, set that a limit was reached and make a note of which one.
//                if (outsideLimitMask&(1<<j))
//                {
//                    status[j] = lm::message::WorkUnitStatus::LIMIT_REACHED;
//                    limitTypeReached[j] = l.type;
//                }
//                else
//                {
//                    // Otherwise set that we finished steps.
//                    status[j] = lm::message::WorkUnitStatus::STEPS_FINISHED;
//                }
//            }
//            return true;
//        }
//    }
//    return false;
//}


//}
//}

#endif
