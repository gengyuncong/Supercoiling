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

#include "lm/Print.h"

#include "lm/simulation/SimulationSupervisor.h"
#include "lm/simulation/StagedSimulationSupervisor.h"

namespace lm {
namespace simulation {

StagedSimulationSupervisor::StagedSimulationSupervisor()
:currentStage(0),currentPhase(0)
{
}

StagedSimulationSupervisor::~StagedSimulationSupervisor()
{
}

int StagedSimulationSupervisor::getCurrentStage()
{
    return currentStage;
}

int StagedSimulationSupervisor::getCurrentPhase()
{
    return currentPhase;
}

void StagedSimulationSupervisor::startSimulation()
{
    // Call the base class method.
    SimulationSupervisor::startSimulation();

    // Start the first stage.
    currentStage = 0;
    startStage(currentStage);

    // Start the first phase.
    currentPhase = 0;
    startPhase(currentStage, currentPhase);
}

bool StagedSimulationSupervisor::assignWork()
{
    // Try to assign work for the phase.
    if (!assignWork(currentStage, currentPhase))
    {
        // Signal that the phase is finished.
        finishPhase(currentStage, currentPhase);

        // We are done with the phase, increment the phase counter.
        currentPhase++;

        // See if we have another phase.
        if (currentPhase < getNumberPhases(currentStage))
        {
            // Start the next phase.
            startPhase(currentStage, currentPhase);

            // Assign some work.
            return assignWork();
        }
        else
        {
            // Otherwise, we are finished with the stage.
            finishStage(currentStage);
            currentStage++;

            // See if we have another stage.
            if (currentStage < getNumberStages())
            {
                // Start the next stage.
                startStage(currentStage);

                // Start the first phase.
                currentPhase = 0;
                startPhase(currentStage, currentPhase);
                return assignWork();

            }
            else
            {
                // Otherwise, we are done with the simulation.
                return false;
            }
        }
    }
    return true;
}

void StagedSimulationSupervisor::finishSimulation()
{
    // Call the base class method.
    SimulationSupervisor::finishSimulation();
}

void StagedSimulationSupervisor::printPerformanceStatistics()
{
    SimulationSupervisor::printPerformanceStatistics();
    Print::printf(Print::INFO, "Simulation status: stage %d of %d, phase %d of %d", currentStage, getNumberStages(), currentPhase, getNumberPhases(currentStage));
}


}
}
