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

#ifndef LM_CME_GILLESPIEDSOLVER_H_
#define LM_CME_GILLESPIEDSOLVER_H_

#include "lm/cme/CMESolver.h"

namespace lm {
namespace cme {

class GillespieDSolver : public CMESolver
{
public:
    static bool registered;
    static bool registerClass();
    static void* allocateObject();

public:
    GillespieDSolver();
    virtual ~GillespieDSolver();

    virtual void reset();
    virtual void getState(lm::io::TrajectoryState* state, uint trajectoryNumber=0);
    virtual void setState(const lm::io::TrajectoryState& state, uint trajectoryNumber=0);
    virtual uint64_t generateTrajectory(uint64_t maxSteps);

protected:
    virtual void allocateRngBuffers();
    virtual void deallocateRngBuffers();
    virtual void updateAllPropensities();
    inline void updatePropensities(uint r);

protected:
    double* rngValues;
    double* expRngValues;
    size_t nextRngValue;
    double * propensities;
};

}
}

#endif /* LM_CME_GILLESPIEDSOLVER_H_ */
