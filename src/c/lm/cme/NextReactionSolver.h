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

#ifndef LM_CME_NEXTREACTIONSOLVER_H_
#define LM_CME_NEXTREACTIONSOLVER_H_

#include "lm/cme/CMESolver.h"
#include "lm/reaction/ReactionQueue.h"

using lm::reaction::ReactionQueue;

namespace lm {
namespace cme {

class NextReactionSolver : public CMESolver
{
public:
    static bool registered;
    static bool registerClass();
    static void* allocateObject();

public:
    NextReactionSolver();
    NextReactionSolver(RandomGenerator::Distributions neededDists);
    virtual ~NextReactionSolver();

    virtual void reset();
    virtual void getState(lm::io::TrajectoryState* state, uint trajectoryNumber=0);
    virtual void setState(const lm::io::TrajectoryState& state, uint trajectoryNumber=0);
    virtual uint64_t generateTrajectory(uint64_t maxSteps);

protected:
    // These methods are not virtual, because they are called in the constructor and destructor.
    void allocateRngBuffers();
    void deallocateRngBuffers();

    void updateAllReactionEventsInQueue();
    void updateReactionEventsInQueue(uint sourceReaction);

protected:
    double* expRngValues;
    size_t nextExpRngValue;
    ReactionQueue * reactionQueue;
};

}
}

#endif
