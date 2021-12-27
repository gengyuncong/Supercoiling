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


#include "lm/Math.h"
#include "lm/ffpilot/FFPilotMath.h"
#include "lm/ffpilot/StageData.h"

namespace lm {
namespace ffpilot {

StageData::StageData()
:numberPhases(0),phaseData(NULL)
{
}

StageData::~StageData()
{
    if (phaseData != NULL) delete[] phaseData; phaseData = NULL;
}

void StageData::initialize(size_t numberPhases, bool pilot, bool forward, size_t linkedProductionStage)
{
    this->numberPhases = numberPhases;
    this->pilot = pilot;
    this->forward = forward;
    this->linkedProductionStage = linkedProductionStage;

    if (!this->pilot && linkedProductionStage != std::numeric_limits<size_t>::max()) THROW_EXCEPTION(lm::RuntimeException, "an ffpilot production stage cannot have a linked production phase");

    if (phaseData != NULL) THROW_EXCEPTION(lm::RuntimeException, "attempted to reinitialize an ffpilot stage");
    phaseData = new PhaseData[numberPhases];

    // Make sure the statistics arrays are the right size.
    phaseCosts.reshape(utuple(static_cast<uint>(numberPhases)));
    phaseWeights.reshape(utuple(static_cast<uint>(numberPhases)));
    phaseWeightVariances.reshape(utuple(static_cast<uint>(numberPhases)));
}

size_t StageData::getNumberPhases()
{
    return numberPhases;
}

PhaseData& StageData::getPhaseData(size_t phaseIndex)
{
    if (phaseData == NULL) THROW_EXCEPTION(lm::RuntimeException, "phase data not initialized");
    if (phaseIndex >= numberPhases) THROW_EXCEPTION(lm::RuntimeException, "phase index %d greater than the number of phases %d", phaseIndex, numberPhases);

    return phaseData[phaseIndex];
}

PhaseData& StageData::getPhaseData(int phaseIndex)
{
    return getPhaseData(static_cast<size_t>(phaseIndex));
}

void StageData::calculateStageStatistics()
{
    // Call the statistics for each phase.
    for (size_t i=0; i<numberPhases; i++)
    {
        // Calculate the phase cost.
        ndarray<uint64_t> costs = phaseData[i].getCosts();
        phaseCosts[static_cast<uint>(i)] = costs.mean();

        // Calculate the phase weights.
        ndarray<double> values = phaseData[i].getValues();
        phaseWeights[static_cast<uint>(i)] = values.mean();
        phaseWeightVariances[static_cast<uint>(i)] = values.variance();

        // If this is the pilot stage and not phase zero, change the phase weights to their lower bound estimate.
        if (pilot && i > 0)
        {
            // Get the estimate.
            phaseWeights[static_cast<uint>(i)] = bernouliConfidenceIntervalLowerBound(phaseWeights[static_cast<uint>(i)], static_cast<double>(phaseData[i].trajectoryData.size()), 0.99);

            // Make sure the estimate is bounded at some minimum value.
            phaseWeights[static_cast<uint>(i)] = max(phaseWeights[static_cast<uint>(i)], 1e-4);
        }
    }
}

}
}
