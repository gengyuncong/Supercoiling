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
 * Author(s): Elijah Roberts, Max Klein
 */

#include <cmath>
#include <map>
#include <set>

#include "lm/cme/CMETrajectory.h"
#include "lm/types/OrderParameters.pb.h"
#include "lm/input/CMERestart.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/io/OrderParameterFirstPassageTimes.pb.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/trajectory/Trajectory.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using std::map;
using std::set;
using robertslab::pbuf::NDArraySerializer;

namespace lm {
namespace cme {

CMETrajectory::CMETrajectory(uint64_t id)
:Trajectory(id)
{
    // Initialize the total step count.
    state.mutable_cme_state()->set_total_steps(0);
}

CMETrajectory::~CMETrajectory()
{
}

void CMETrajectory::initializeSpeciesCountsState(const lm::input::ReactionModel& rm)
{
    // Clear the species counts.
    state.mutable_cme_state()->mutable_species_counts()->Clear();

    // Initialize the species counts.
    if (rm.initial_species_count().size() != static_cast<int>(rm.number_species())) THROW_EXCEPTION(lm::RuntimeException, "wrong number of initial species counts: %d %d", rm.initial_species_count().size(), rm.number_species());
    lm::io::SpeciesCounts* sc = state.mutable_cme_state()->mutable_species_counts();
    sc->set_trajectory_id(id);
    sc->set_number_entries(1);
    sc->set_number_species(static_cast<int>(rm.number_species()));
    for (int j=0; j<rm.initial_species_count().size(); j++)
    {
        sc->add_species_count(static_cast<int>(rm.initial_species_count(j)));
    }
    sc->add_time(0.0);
}

void CMETrajectory::initializeSpeciesCountsState(const ndarray<int32_t>& initialCounts, double time)
{
    // Clear the species counts.
    state.mutable_cme_state()->mutable_species_counts()->Clear();

    // Initialize the species counts.
    lm::io::SpeciesCounts* sc = state.mutable_cme_state()->mutable_species_counts();
    sc->set_trajectory_id(id);
    sc->set_number_entries(1);
    sc->set_number_species(int(initialCounts.shape[0]));
    for (uint j=0; j<initialCounts.shape[0]; j++)
    {
        sc->add_species_count(initialCounts[j]);
    }
    sc->add_time(time);
}

void CMETrajectory::initializeDegreeAdvancementsState(const lm::input::ReactionModel& rm)
{
    // Initialize the degree advancements to zero.
    ndarray<uint64_t> initalDegreeAdvancementCounts(utuple(rm.number_reactions()));
    NDArraySerializer::serializeInto(state.mutable_cme_state()->mutable_degree_advancements(), initalDegreeAdvancementCounts);
}

void CMETrajectory::initializeOrderParametersState(const lm::types::OrderParameters& op)
{
    lm::io::OrderParametersValues* opv = state.mutable_cme_state()->mutable_order_parameter_values();
    opv->set_trajectory_id(id);
    opv->set_number_entries(0);
    opv->set_number_order_parameters(op.order_parameter().size());
}

void CMETrajectory::initializeSpeciesFirstPassageTimeState(uint32_t speciesToTrack)
{
    lm::io::FirstPassageTimes* fpt = state.mutable_cme_state()->add_first_passage_times();
    fpt->set_trajectory_id(id);
    fpt->set_species(speciesToTrack);
    ndarray<int32_t> counts(utuple(0));
    ndarray<double> times(utuple(0));
    NDArraySerializer::serializeInto(fpt->mutable_counts(), counts);
    NDArraySerializer::serializeInto(fpt->mutable_first_passage_times(), times);
}

void CMETrajectory::initializeOrderParameterFirstPassageTimeState(uint32_t orderParameterToTrack)
{
    lm::io::OrderParameterFirstPassageTimes* opFPT = state.mutable_cme_state()->add_order_parameter_first_passage_times();
    opFPT->set_trajectory_id(id);
    opFPT->set_order_parameter_id(orderParameterToTrack);
    ndarray<int32_t> values(utuple(0));
    ndarray<double> times(utuple(0));
    NDArraySerializer::serializeInto(opFPT->mutable_order_parameter_value(), values);
    NDArraySerializer::serializeInto(opFPT->mutable_first_passage_time(), times);
}

}
}
