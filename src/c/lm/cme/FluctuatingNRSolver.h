/*
 * Copyright 2019 Johns Hopkins University
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
 * Author(s): Elijah Roberts
 */

#ifndef LM_CME_FLUCTUATINGNRSOLVER_H_
#define LM_CME_FLUCTUATINGNRSOLVER_H_

#include "lm/cme/NextReactionSolver.h"
#include "lm/me/TimeSeriesList.h"
#include "lm/rng/RandomGenerator.h"
#include "robertslab/Types.h"

using lm::reaction::ReactionQueue;
using lm::rng::RandomGenerator;

namespace lm {
namespace cme {

class StochasticProcess;

class FluctuatingNRSolver : public NextReactionSolver
{
public:
    static bool registered;
    static bool registerClass();
    static void* allocateObject();

public:
    FluctuatingNRSolver();
    virtual ~FluctuatingNRSolver();

    virtual void setOutputOptions(const lm::input::OutputOptions& outputOptions);
    virtual void setReactionModel(const lm::input::ReactionModel& rm);
    virtual void reset();
    virtual void getState(lm::io::TrajectoryState* state, uint trajectoryNumber=0);
    virtual void setState(const lm::io::TrajectoryState& state, uint trajectoryNumber=0);
    virtual void trajectoryStarted();
    virtual void timeUpdated();
    virtual void trajectoryFinished();
    virtual uint64_t generateTrajectory(uint64_t maxSteps);

protected:
    // These methods are not virtual, because they are called in the constructor and destructor.
    void allocateRngBuffers();
    void deallocateRngBuffers();

    void updateReactionPropensities();
    void rescaleReactionEventsInQueue();

protected:
    double* normRngValues;
    size_t nextNormRngValue;

    // Variables for storing the depencdncies.
    size_t numberNoisyReactions;
    uint* noisyReactions;
    vector<utuple> noisyReactionDependencies;
    vector<double> noisyReactionOriginalConstants;

    double nextProcessUpdateTime;
    double processUpdateInterval;
    size_t numberProcesses;
    StochasticProcess** processes;

    bool writeProcessTimeSeries;
    double processWriteInterval;
    lm::me::TimeSeriesList<double> processTimeSeries;
};

class StochasticProcess
{
public:
    StochasticProcess();
    virtual ~StochasticProcess();
    virtual void reset(double normRng)=0;
    virtual void update(double dt, double normRng)=0;
    double value;
};

class OrnsteinUhlenbeckProcess : public StochasticProcess
{
public:
    OrnsteinUhlenbeckProcess(double variance, double tau);
    virtual ~OrnsteinUhlenbeckProcess();
    virtual void reset(double normRng);
    virtual void update(double dt, double normRng);
private:
    double variance;
    double tau;
    double _value;
};

class LogOrnsteinUhlenbeckProcess : public StochasticProcess
{
public:
    LogOrnsteinUhlenbeckProcess(double variance, double tau);
    virtual ~LogOrnsteinUhlenbeckProcess();
    virtual void reset(double normRng);
    virtual void update(double dt, double normRng);
private:
    double variance;
    double tau;
    double _value;
};

}
}

#endif
