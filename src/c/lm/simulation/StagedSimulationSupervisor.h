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

#ifndef LM_SIMULATION_PHASEDSIMULATIONSUPERVISOR_H
#define LM_SIMULATION_PHASEDSIMULATIONSUPERVISOR_H

#include "lm/simulation/SimulationSupervisor.h"

namespace lm {
namespace simulation {

class StagedSimulationSupervisor : public lm::simulation::SimulationSupervisor
{
public:
    StagedSimulationSupervisor();
    virtual ~StagedSimulationSupervisor();
    virtual void printPerformanceStatistics();

protected:

    virtual int getNumberStages()=0;
    virtual int getNumberPhases(int stage)=0;
    virtual int getCurrentStage();
    virtual int getCurrentPhase();
    virtual void startSimulation();
    virtual bool assignWork();
    virtual void finishSimulation();
    virtual void startStage(int stage)=0;
    virtual void startPhase(int stage, int phase)=0;
    virtual bool assignWork(int stage, int phase)=0;
    virtual void finishPhase(int stage, int phase)=0;
    virtual void finishStage(int stage)=0;

protected:
    int currentStage;
    int currentPhase;
};

}
}

#endif
