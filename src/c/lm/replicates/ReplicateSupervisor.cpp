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

#include <map>
#include <string>
#include <vector>

#include "hrtime.h"
#include "lm/ClassFactory.h"
#include "lm/Print.h"
#include "lm/cme/CMETrajectory.h"
#include "lm/io/OutputWriter.h"
#include "lm/main/Globals.h"
#include "lm/rdme/RDMETrajectory.h"
#include "lm/replicates/ReplicateSupervisor.h"
#include "lm/replicates/ReplicateTrajectoryList.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using std::map;
using std::string;
using std::vector;
using robertslab::pbuf::NDArraySerializer;

namespace lm {
namespace replicates {

bool ReplicateSupervisor::registered=ReplicateSupervisor::registerClass();

bool ReplicateSupervisor::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::simulation::SimulationSupervisor","lm::replicates::ReplicateSupervisor",&ReplicateSupervisor::allocateObject);
    return true;
}

void* ReplicateSupervisor::allocateObject()
{
    return new ReplicateSupervisor();
}

ReplicateSupervisor::ReplicateSupervisor()
:numberReplicates(0)
{
}

std::string ReplicateSupervisor::getClassName()
{
    return "lm::replicates::ReplicateSupervisor";
}

ReplicateSupervisor::~ReplicateSupervisor()
{
}

void ReplicateSupervisor::startSimulation()
{
    // Call the base class method
    TrajectoryListSupervisor::startSimulation();

    // Turn on writing of the first and last time points.
    input->mutable_output_options()->set_write_initial_trajectory_state(true);
    input->mutable_output_options()->set_write_final_trajectory_state(true);

    // Set any simulation options.
    input->mutable_simulation_options()->set_parts_per_work_unit(replicateBatchSize);

    Print::printf(Print::INFO, "Replicate supervisor started.");
}

void ReplicateSupervisor::buildTrajectoryList()
{
    // Create the new trajectory list.
    trajectoryList = new ReplicateTrajectoryList();

    // If we have restart data, get the arrays.
    ndarray<uint64_t>* restart_trajectories = NULL;
    ndarray<double>* restart_times = NULL;
    ndarray<int32_t>* restart_counts = NULL;
    if (input->has_cme_restart() && input->cme_restart().has_restart_trajectory_ids() && input->cme_restart().has_restart_times() && input->cme_restart().has_restart_species_counts())
    {
        // Unpack the arrays.
        restart_trajectories = NDArraySerializer::deserializeAllocate<uint64_t>(input->cme_restart().restart_trajectory_ids());
        restart_times = NDArraySerializer::deserializeAllocate<double>(input->cme_restart().restart_times());
        restart_counts = NDArraySerializer::deserializeAllocate<int32_t>(input->cme_restart().restart_species_counts());
        Print::printf(Print::INFO, "Replicate supervisor loaded restart data for %d replicates.", restart_trajectories->shape[0]);
    }

    // Add the trajectories to the list.
    for (vector<uint64_t>::iterator it = ::replicates.begin(); it != ::replicates.end(); it++)
    {
        uint64_t replicateId = *it;

        if (input->has_reaction_model())
        {
            // Create a trajectory.
            lm::cme::CMETrajectory* trajectory;
            lm::rdme::RDMETrajectory* rdmeTrajectory = NULL;

            // Create the right type of trajectory object.
            if (input->has_diffusion_model())
            {
                rdmeTrajectory = new lm::rdme::RDMETrajectory(replicateId);
                trajectory = rdmeTrajectory;
            }
            else
            {
                trajectory = new lm::cme::CMETrajectory(replicateId);
            }

            // Initialize the trajectory from the model.
            trajectory->initializeSpeciesCountsState(input->reaction_model());
            if (input->has_diffusion_model() && rdmeTrajectory != NULL)
                rdmeTrajectory->initializeLatticeState(input->diffusion_model());

            // If we have restart data, overwrite any necessary trajectory initialization.
            if (restart_trajectories != NULL)
            {
                for (uint i=0; i<restart_trajectories->shape[0]; i++)
                {
                    if (restart_trajectories->get(i) == replicateId && i < restart_times->shape[0] && i < restart_counts->shape[0])
                    {
                        ndarray<int32_t> counts(utuple(restart_counts->shape[1]));
                        for (uint j=0; j<restart_counts->shape[1]; j++)
                            counts[j] = restart_counts->get(utuple(i,j));
                        trajectory->initializeSpeciesCountsState(counts, restart_times->get(i));

                        // Initialize the restart lattice, if necessary.
                        if (input->has_diffusion_model() && rdmeTrajectory != NULL && input->has_rdme_restart())
                        {
                            // Initialize the particles.
                            ndarray<uint8_t>* restart_lattice_particles = NDArraySerializer::deserializeAllocate<uint8_t>(input->rdme_restart().restart_lattice(static_cast<int>(i)).particles());
                            rdmeTrajectory->initializeLatticeParticlesState(*restart_lattice_particles);
                            delete restart_lattice_particles;

                            // Initialize the site, if they are present.
                            if (input->rdme_restart().restart_lattice(static_cast<int>(i)).has_sites())
                            {
                                ndarray<uint8_t>* restart_lattice_sites = NDArraySerializer::deserializeAllocate<uint8_t>(input->rdme_restart().restart_lattice(static_cast<int>(i)).sites());
                                rdmeTrajectory->initializeLatticeSitesState(*restart_lattice_sites);
                                delete restart_lattice_sites;
                            }
                        }
                        break;
                    }
                }
            }

            // Initialize any degree advancements state.
            if (input->has_output_options() && input->output_options().has_degree_advancement_write_interval()) trajectory->initializeDegreeAdvancementsState(input->reaction_model());

            // Initialize any order parameters state.
            if (input->has_order_parameters()) trajectory->initializeOrderParametersState(input->order_parameters());

            // Initialize any fpt tracking state.
            if (input->has_output_options())
            {
                for (int i=0; i<input->output_options().fpt_species_to_track().size(); i++)
                    trajectory->initializeSpeciesFirstPassageTimeState(input->output_options().fpt_species_to_track(i));
                for (int i=0; i<input->output_options().fpt_order_parameter_to_track().size(); i++)
                    trajectory->initializeOrderParameterFirstPassageTimeState(input->output_options().fpt_order_parameter_to_track(i));
            }

            // Add the trajectory to the list.
            trajectoryList->addTrajectory(replicateId, trajectory);
        }
    }

    numberReplicates += trajectoryList->size();

    // Free any data.
    if (restart_counts != NULL) delete restart_counts;
    if (restart_times != NULL) delete restart_times;
    if (restart_trajectories != NULL) delete restart_trajectories;
}

void ReplicateSupervisor::buildTrajectoryOptions()
{
    trajectoryLimits.Clear();
    trajectoryBarriers.Clear();

    // See if we have a time limit.
    if (input->has_simulation_options() && input->simulation_options().has_time_limit())
    {
        // Add the time limit to the trajectory limits.
        lm::types::TrajectoryLimit* limit = trajectoryLimits.add_limits();
        limit->set_id(trajectoryLimits.limits_size()-1);
        limit->set_type(lm::types::TrajectoryLimit::TIME);
        limit->set_stopping_condition(lm::types::TrajectoryLimit::MAX_INCLUSIVE);
        limit->set_stopping_value_double(input->simulation_options().time_limit());
    }

    // See if we have any species upper limits.
    if (input->has_simulation_options() && input->simulation_options().species_upper_limit().size() > 0)
    {
        for (int i=0; i<input->simulation_options().species_upper_limit().size(); i++)
        {
            lm::input::SpeciesLimit speciesLimit = input->simulation_options().species_upper_limit(i);

            // Add the species limit to the trajectory limits.
            lm::types::TrajectoryLimit* trajLimit = trajectoryLimits.add_limits();
            trajLimit->set_id(trajectoryLimits.limits_size()-1);
            trajLimit->set_type(lm::types::TrajectoryLimit::SPECIES);
            trajLimit->set_type_arg(speciesLimit.species_index());
            trajLimit->set_stopping_condition(lm::types::TrajectoryLimit::MAX_INCLUSIVE);
            trajLimit->set_stopping_value_int(speciesLimit.limit_value());
        }
    }

    // See if we have any species lower limits.
    if (input->has_simulation_options() && input->simulation_options().species_lower_limit().size() > 0)
    {
        for (int i=0; i<input->simulation_options().species_lower_limit().size(); i++)
        {
            lm::input::SpeciesLimit speciesLimit = input->simulation_options().species_lower_limit(i);

            // Add the species limit to the trajectory limits.
            lm::types::TrajectoryLimit* trajLimit = trajectoryLimits.add_limits();
            trajLimit->set_id(trajectoryLimits.limits_size()-1);
            trajLimit->set_type(lm::types::TrajectoryLimit::SPECIES);
            trajLimit->set_type_arg(speciesLimit.species_index());
            trajLimit->set_stopping_condition(lm::types::TrajectoryLimit::MIN_INCLUSIVE);
            trajLimit->set_stopping_value_int(speciesLimit.limit_value());
        }
    }

    // See if we have any species reflecting barriers.
    if (input->has_simulation_options() && input->simulation_options().species_reflecting_barrier().size() > 0)
    {
        for (int i=0; i<input->simulation_options().species_reflecting_barrier().size(); i++)
        {
            lm::input::SpeciesLimit speciesLimit = input->simulation_options().species_reflecting_barrier(i);

            // Add the species barrier to the trajectory barriers.
            lm::types::TrajectoryBarrier* trajBarrier = trajectoryBarriers.add_barrier();
            trajBarrier->set_type(lm::types::TrajectoryBarrier::SPECIES);
            trajBarrier->set_type_arg(speciesLimit.species_index());
            trajBarrier->set_behavior(lm::types::TrajectoryBarrier::REFLECTING);
            trajBarrier->set_behavior_value_int(speciesLimit.limit_value());
        }
    }

    // See if we have any order parameter reflecting barriers.
    if (input->has_simulation_options() && input->simulation_options().order_parameter_reflecting_barrier().size() > 0)
    {
        for (int i=0; i<input->simulation_options().order_parameter_reflecting_barrier().size(); i++)
        {
            lm::input::OrderParameterLimit opLimit = input->simulation_options().order_parameter_reflecting_barrier(i);

            // Add the order parameter barrier to the trajectory barriers.
            lm::types::TrajectoryBarrier* trajBarrier = trajectoryBarriers.add_barrier();
            trajBarrier->set_type(lm::types::TrajectoryBarrier::ORDER_PARAMETER);
            trajBarrier->set_type_arg(opLimit.order_parameter_index());
            trajBarrier->set_behavior(lm::types::TrajectoryBarrier::REFLECTING);
            trajBarrier->set_behavior_value_double(opLimit.limit_value());
        }
    }

    // See if we have any species tracking barriers.
    if (input->has_simulation_options() && input->simulation_options().species_tracking_barrier().size() > 0)
    {
        for (int i=0; i<input->simulation_options().species_tracking_barrier().size(); i++)
        {
            lm::input::SpeciesLimit speciesLimit = input->simulation_options().species_tracking_barrier(i);

            // Add the species barrier to the trajectory barriers.
            lm::types::TrajectoryBarrier* trajBarrier = trajectoryBarriers.add_barrier();
            trajBarrier->set_type(lm::types::TrajectoryBarrier::SPECIES);
            trajBarrier->set_type_arg(speciesLimit.species_index());
            trajBarrier->set_behavior(lm::types::TrajectoryBarrier::TRACKING);
            trajBarrier->set_behavior_value_int(speciesLimit.limit_value());
        }
    }

    // See if we have any order parameter tracking barriers.
    if (input->has_simulation_options() && input->simulation_options().order_parameter_tracking_barrier().size() > 0)
    {
        for (int i=0; i<input->simulation_options().order_parameter_tracking_barrier().size(); i++)
        {
            lm::input::OrderParameterLimit opLimit = input->simulation_options().order_parameter_tracking_barrier(i);

            // Add the order parameter barrier to the trajectory barriers.
            lm::types::TrajectoryBarrier* trajBarrier = trajectoryBarriers.add_barrier();
            trajBarrier->set_type(lm::types::TrajectoryBarrier::ORDER_PARAMETER);
            trajBarrier->set_type_arg(opLimit.order_parameter_index());
            trajBarrier->set_behavior(lm::types::TrajectoryBarrier::TRACKING);
            trajBarrier->set_behavior_value_double(opLimit.limit_value());
        }
    }

    // See if we have any barrier crossing limits.
    if (input->has_simulation_options() && input->simulation_options().tracking_barrier_crossing_limit().size() > 0)
    {
        for (int i=0; i<input->simulation_options().tracking_barrier_crossing_limit().size(); i++)
        {
            lm::input::BarrierLimit speciesLimit = input->simulation_options().tracking_barrier_crossing_limit(i);

            // Add the crossing limit to the trajectory limits.
            lm::types::TrajectoryLimit* trajLimit = trajectoryLimits.add_limits();
            trajLimit->set_id(trajectoryLimits.limits_size()-1);
            trajLimit->set_type(lm::types::TrajectoryLimit::BARRIER_CROSSING);
            trajLimit->set_type_arg(speciesLimit.barrier_index()+static_cast<uint32_t>(input->simulation_options().species_reflecting_barrier().size()+input->simulation_options().order_parameter_reflecting_barrier().size()));
            trajLimit->set_stopping_condition(lm::types::TrajectoryLimit::MAX_INCLUSIVE);
            trajLimit->set_stopping_value_int(speciesLimit.limit_value());

            if (trajLimit->type_arg() >= static_cast<uint32_t>(trajectoryBarriers.barrier().size())) THROW_EXCEPTION(lm::RuntimeException, "Invalid barrier crossing limit, refers to nonexistent barrier.");
            if (trajectoryBarriers.barrier(static_cast<int32_t>(trajLimit->type_arg())).behavior() != lm::types::TrajectoryBarrier::TRACKING) THROW_EXCEPTION(lm::RuntimeException, "Invalid barrier crossing limit, refers to nontracking barrier.");
        }
    }



}

void ReplicateSupervisor::finishSimulation()
{
    // Call the base class method
    TrajectoryListSupervisor::finishSimulation();

    Print::printf(Print::INFO, "Replicate supervisor finished %lld replicates in %0.3e seconds.", numberReplicates, convertHrToSeconds(getHrTime()-simulationStartTime));
}

void ReplicateSupervisor::printPerformanceStatistics()
{
    TrajectoryListSupervisor::printPerformanceStatistics();
}

}
}
