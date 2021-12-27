/*
 * Copyright 2016-2018 Johns Hopkins University
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

#include <cmath>
#include <limits>
#include <map>
#include <string>

#include "hrtime.h"
#include "lm/ClassFactory.h"
#include "lm/Math.h"
#include "lm/Print.h"
#include "lm/cme/CMETrajectory.h"
#include "lm/types/TrajectoryLimits.pb.h"
#include "lm/io/OutputWriter.h"
#include "lm/io/TrajectoryState.pb.h"
#include "lm/main/Globals.h"
#include "lm/message/Message.pb.h"
#include "lm/message/FinishedWorkUnit.pb.h"
#include "lm/message/RunWorkUnit.pb.h"
#include "lm/message/StartedWorkUnit.pb.h"
#include "lm/microenv/METrajectoryList.h"
#include "lm/microenv/MicroenvironmentSupervisor.h"
#include "lm/microenv/PDETrajectoryList.h"
#include "lm/pde/DiffusionPDETrajectory.h"
#include "lm/resource/ResourceMap.h"
#include "lm/simulation/SimulationSupervisor.h"
#include "lm/slot/SlotList.h"
#include "robertslab/Types.h"
#include "robertslab/pbuf/NDArraySerializer.h"

#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"

using std::map;
using std::string;
using lm::message::Communicator;
using lm::message::Endpoint;
using lm::resource::ResourceMap;

namespace lm {
namespace microenv {

bool MicroenvironmentSupervisor::registered=MicroenvironmentSupervisor::registerClass();

bool MicroenvironmentSupervisor::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::simulation::SimulationSupervisor","lm::microenv::MicroenvironmentSupervisor",&MicroenvironmentSupervisor::allocateObject);
    return true;
}

void* MicroenvironmentSupervisor::allocateObject()
{
    return new MicroenvironmentSupervisor();
}

MicroenvironmentSupervisor::MicroenvironmentSupervisor()
:numberReplicates((uint)::replicates.size()),currentReplicateIndex(0),numberTimesteps(0),currentTimestep(0),tau(0.0),maxTime(0.0),
meTrajectoryList(NULL),pdeTrajectoryList(NULL),
gridSpacing(0.0),gridElementVolume(0.0),numberCells(0),cellCoordinates(NULL),cellGridPoints(NULL),cellVolumes(NULL),cellPreviousCounts(NULL),cellCurrentCounts(NULL),cellFlux(NULL),
stats_pdeWorkUnitsSteps(0),stats_pdeWorkUnitsTime(0.0),stats_timesteps(0),stats_timestepStartTime(0),stats_timestepTotalTime(0),stats_timestepPDETime(0),stats_timestepMETime(0),stats_timestepReconcileTime(0)
{
}

MicroenvironmentSupervisor::~MicroenvironmentSupervisor()
{
    if (meTrajectoryList != NULL) delete meTrajectoryList; meTrajectoryList = NULL;
    if (pdeTrajectoryList != NULL) delete pdeTrajectoryList; pdeTrajectoryList = NULL;
    if (cellCoordinates != NULL) delete cellCoordinates; cellCoordinates = NULL;
    if (cellGridPoints != NULL) delete cellGridPoints; cellGridPoints = NULL;
    if (cellVolumes != NULL) delete cellVolumes; cellVolumes = NULL;
    if (cellPreviousCounts != NULL) delete cellPreviousCounts; cellPreviousCounts = NULL;
    if (cellCurrentCounts != NULL) delete cellCurrentCounts; cellCurrentCounts = NULL;
    if (cellFlux != NULL) delete cellFlux; cellFlux = NULL;
}

std::string MicroenvironmentSupervisor::getClassName()
{
    return "lm::microenv::MicroenvironmentSupervisor";
}

void MicroenvironmentSupervisor::startSimulation()
{
    // Get the tau.
    if (!input->has_microenv_model()) THROW_EXCEPTION(lm::RuntimeException, "MicroenvironmentSupervisor requires microenv_model in the input");
    tau = input->microenv_model().synchronization_timestep();

    // Get the grid properties.
    utuple gridShape(uint(input->microenv_model().grid_shape().size()), (const uint32_t*)(input->microenv_model().grid_shape().data()));
    gridSpacing = input->microenv_model().grid_spacing();
    gridElementVolume = pow(gridSpacing, 3.0)*1000;

    // Get the cell coordinates and volume.
    numberCells = input->microenv_model().number_cells();
    if (numberCells > 0)
    {
        cellCoordinates = robertslab::pbuf::NDArraySerializer::deserializeAllocate<double>(input->microenv_model().cell_coordinates());
        cellGridPoints = new ndarray<uint32_t>(utuple(numberCells,3));
        cellVolumes = robertslab::pbuf::NDArraySerializer::deserializeAllocate<double>(input->microenv_model().cell_volume());
        cellPreviousCounts = new ndarray<int32_t>(utuple(numberCells,1));
        cellCurrentCounts = new ndarray<int32_t>(utuple(numberCells,1));
        cellFlux = new ndarray<int32_t>(utuple(numberCells,1));

        // Calculate the cell grid points.
        for (uint32_t i=0; i<numberCells; i++)
        {
            // Get the nearest grid point to the cell.
            uint32_t x = uint32_t(round((*cellCoordinates)[utuple(i,0U)]/gridSpacing));
            uint32_t y = uint32_t(round((*cellCoordinates)[utuple(i,1U)]/gridSpacing));
            uint32_t z = uint32_t(round((*cellCoordinates)[utuple(i,2U)]/gridSpacing));

            // Verify that the cell falls in the diffusion grid.
            if (x >= gridShape[0] || y >= gridShape[1] || z >= gridShape[2])
            {
                Print::printf(Print::FATAL, "Cell %d was located at %e,%e,%e (%d,%d,%d), which is off the diffusion grid (%d,%d,%d).", i, (*cellCoordinates)[utuple(i,0U)], (*cellCoordinates)[utuple(i,1U)], (*cellCoordinates)[utuple(i,2U)], x, y, z, gridShape[0], gridShape[1], gridShape[2]);
                throw RuntimeException("MicroenvironmentSupervisor encountered a critical error in the configuration.");
            }

            // Save the cell's grid point.
            (*cellGridPoints)[utuple(i,0U)] = x;
            (*cellGridPoints)[utuple(i,1U)] = y;
            (*cellGridPoints)[utuple(i,2U)] = z;
        }
    }

    // Figure out how many timesteps we need to perform.
    if (!input->has_simulation_options() || !input->simulation_options().has_time_limit()) throw RuntimeException("MicroenvironmentSupervisor requires max_time in the input");
    maxTime = input->simulation_options().time_limit();
    numberTimesteps = uint(ceil((maxTime/tau)-EPS));

    // Turn on writing of the first time point.
    input->mutable_output_options()->set_write_initial_trajectory_state(true);

    PROF_BEGIN(PROF_MENV_RUN_SIM);
    simulationStartTime=getHrTime();

    Print::printf(Print::INFO, "Microenvironment supervisor started: %d replicates", numberReplicates);

    // Initialize the first replicate.
    startReplicate();

    // Call the base class method
    SimulationSupervisor::startSimulation();
}

void MicroenvironmentSupervisor::startReplicate()
{
    PROF_BEGIN(PROF_MENV_START_REPLICATE);

    // Set the output prefix to the replicate number.
    input->mutable_output_options()->set_record_name_prefix("/Replicates/"+std::to_string(::replicates[currentReplicateIndex]));

    // Build the trajectory lists.
    buildTrajectoryLists();

    // Rest the timestep counter.
    currentTimestep = 0;

    // If we have cells, initialize them from the diffusion grid.
    if (numberCells > 0)
    {
        // Initialize the diffusing species counts from the diffusion grid.
        (*cellFlux) = 0;
        pdeTrajectoryList->reconcileDiffusionGrid(gridElementVolume, cellGridPoints, cellVolumes, cellPreviousCounts, cellFlux, 0);

        // Copy the current counts of the diffusing species.
        meTrajectoryList->copySpeciesCountFrom(*cellPreviousCounts, 0, 0);
    }

    PROF_END(PROF_MENV_START_REPLICATE);
}

void MicroenvironmentSupervisor::buildTrajectoryLists()
{
    // Free the old trajectory lists, if they exist.
    if (meTrajectoryList != NULL) delete meTrajectoryList; meTrajectoryList = NULL;
    if (pdeTrajectoryList != NULL) delete pdeTrajectoryList; pdeTrajectoryList = NULL;

    // Allocate the new lists.
    pdeTrajectoryList = new PDETrajectoryList();
    meTrajectoryList = new METrajectoryList();

    // Create the diffusion pde trajectory and add it to the PDE list.
    lm::pde::DiffusionPDETrajectory* pdeTrajectory = new lm::pde::DiffusionPDETrajectory(numberCells);
    pdeTrajectory->initializeConcentrations(input->microenv_model().initial_concentrations(0));
    pdeTrajectoryList->addTrajectory(pdeTrajectory->getID(), pdeTrajectory);

    // Go through each cell in the microenvironment and add it to the ME list.
    if (numberCells > 0)
    {
        ndarray<uint32_t>* initialCounts = robertslab::pbuf::NDArraySerializer::deserializeAllocate<uint32_t>(input->microenv_model().cell_initial_species_counts());
        uint numberSpecies = initialCounts->shape[1];
        for (uint i=0; i<numberCells; i++)
        {
            // Create a new trajectory for the cell.
            uint64_t id = i;
            lm::cme::CMETrajectory* trajectory = new lm::cme::CMETrajectory(id);

            // Set the initial species counts for the cell.
            ndarray<int32_t> cellInitialCounts(numberSpecies);
            for (uint j=0; j<numberSpecies; j++)
                cellInitialCounts[j] = int32_t(initialCounts->get(utuple(i,j)));
            trajectory->initializeSpeciesCountsState(cellInitialCounts, 0.0);

            meTrajectoryList->addTrajectory(id, trajectory);
        }
        delete initialCounts;
    }
}

void MicroenvironmentSupervisor::continueReplicate()
{
    PROF_BEGIN(PROF_MENV_CONT_REPLICATE);

    hrtime t0 = getHrTime();

    // If we have cells, reconcile them with the diffusion grid.
    if (numberCells > 0)
    {
        // Copy the current counts of the diffusing species.
        meTrajectoryList->copySpeciesCountInto(cellCurrentCounts, 0, 0);

        // Calculate the flux into or out of the diffusion grid over the last timestep.
        cellFlux->equalsDifference(*cellCurrentCounts, *cellPreviousCounts);

        // Go through each cell and reconcile it with the diffusion grid.
        pdeTrajectoryList->reconcileDiffusionGrid(gridElementVolume, cellGridPoints, cellVolumes, cellCurrentCounts, cellFlux, 0);

        // Set the new counts of the difusing species.
        meTrajectoryList->copySpeciesCountFrom(*cellCurrentCounts, 0, 0);

        // Swap the current and previous counts.
        ndarray<int32_t>* tmp = cellPreviousCounts;
        cellPreviousCounts = cellCurrentCounts;
        cellCurrentCounts = tmp;
    }

    // Update the trajectory lists to run for another timestep.
    meTrajectoryList->restartFinishedTrajectories();
    pdeTrajectoryList->restartFinishedTrajectories();

    // Record how long it took to reconcile.
    stats_timestepReconcileTime += getHrTime()-t0;

    PROF_END(PROF_MENV_CONT_REPLICATE);
}

bool MicroenvironmentSupervisor::assignWork()
{
    // Assign any trajectories still waiting to be run for this time step.
    while (slotList->hasFreeSlots() && (meTrajectoryList->anyWaiting() || pdeTrajectoryList->anyWaiting()))
    {
        // Create the run work unit message.
        lm::message::Message msg;
        lm::message::RunWorkUnit* rwuMsg = msg.mutable_run_work_unit();

        // Get the free slot.
        const lm::slot::Slot slot = slotList->getFreeSlot();

        // Build the work unit parts.
        if (meTrajectoryList->anyWaiting())
            buildRunWorkUnit(rwuMsg, true);
        else
            buildRunWorkUnit(rwuMsg, false);

        // Run the work unit.
        slotList->runWorkUnit(communicator, &msg);
    }

    // See if we are done with the timestep.
    if (meTrajectoryList->allFinished() && pdeTrajectoryList->allFinished())
    {
        // Increment the timestep.
        currentTimestep++;

        // See if we are done with the replicate.
        if (currentTimestep >= numberTimesteps)
        {
            // Increment the replicate.
            currentReplicateIndex++;

            // See if we are done with all the replicates.
            if (currentReplicateIndex >= numberReplicates)
            {
                // We are done with all the work.
                return false;
            }
            else
            {
                // Start working on the next replicate.
                startReplicate();

                // Assign work for the new replicate.
                return assignWork();
            }
        }
        else
        {
            // Start the next timestep in the replicate.
            continueReplicate();

            // Assign work for the new timestep.
            return assignWork();
        }
    }
    else
    {
        // Otherwise, we are still working on the timestep.
        return true;
    }
}

void MicroenvironmentSupervisor::buildRunWorkUnit(lm::message::RunWorkUnit* msg, bool me)
{
    PROF_BEGIN(PROF_MENV_BUILD_WORK_UNIT);

    // Call the base class method.
    SimulationSupervisor::buildRunWorkUnitHeader(msg);

    // Set the output options.
    msg->mutable_output_options()->CopyFrom(input->output_options());

    // Set the time limit.
    lm::types::TrajectoryLimit* limit = msg->mutable_trajectory_limits()->add_limits();
    limit->set_id(msg->mutable_trajectory_limits()->limits_size()-1);
    limit->set_type(lm::types::TrajectoryLimit::TIME);
    limit->set_stopping_condition(lm::types::TrajectoryLimit::MAX_INCLUSIVE);
    limit->set_stopping_value_double((currentTimestep+1)*tau);

    if (me)
    {
        // Set the solver to be master equation.
        msg->set_solver_type(lm::types::SolverType::ME);

        // Set the maximum number of steps for the work unit.
        msg->set_max_steps(input->simulation_options().steps_per_work_unit_part());

        // Add the model.
        if (input->has_reaction_model()) msg->mutable_reaction_model()->CopyFrom(input->reaction_model());
        if (input->has_diffusion_model()) msg->mutable_diffusion_model()->CopyFrom(input->diffusion_model());
        if (input->has_order_parameters()) msg->mutable_order_parameters()->CopyFrom(input->order_parameters());

        // Add the parts.
        meTrajectoryList->buildWorkUnitParts(msg->work_unit_id(), msg, 250);
    }
    else
    {
        // Set the solver to be pde.
        msg->set_solver_type(lm::types::SolverType::DIFFUSION_PDE);

        // Set the maximum number of steps for the work unit.
        msg->set_max_steps(1000);

        // Add the model.
        if (input->has_microenv_model()) msg->mutable_microenv_model()->CopyFrom(input->microenv_model());

        // Add the parts.
        pdeTrajectoryList->buildWorkUnitParts(msg->work_unit_id(), msg, 1);
    }

    PROF_END(PROF_MENV_BUILD_WORK_UNIT);
}

void MicroenvironmentSupervisor::receivedFinishedWorkUnit(const lm::message::FinishedWorkUnit& msg)
{
    PROF_BEGIN(PROF_TRAJECTORY_LIST_WORK_UNIT_FINISHED);

    // Update the trajectory lists.
    if (msg.solver_type() == lm::types::SolverType::ME)
        meTrajectoryList->processWorkUnitFinished(msg.work_unit_id(), msg);
    else if (msg.solver_type() == lm::types::SolverType::DIFFUSION_PDE)
        pdeTrajectoryList->processWorkUnitFinished(msg.work_unit_id(), msg);
    else
        throw RuntimeException("MicroenvironmentSupervisor received a finished work unit message that could not be mapped to a trajectory list.");

    PROF_END(PROF_TRAJECTORY_LIST_WORK_UNIT_FINISHED);

    // Call the base class method.
    SimulationSupervisor::receivedFinishedWorkUnit(msg);
}

void MicroenvironmentSupervisor::finishSimulation()
{
    Print::printf(Print::INFO, "MicroenvironmentSupervisor supervisor finished %u timesteps for %u replicates in %0.2f seconds.", numberTimesteps, numberReplicates, convertHrToSeconds(getHrTime()-simulationStartTime));
    PROF_END(PROF_MENV_RUN_SIM);

    // Call the base class method.
    SimulationSupervisor::finishSimulation();
}

void MicroenvironmentSupervisor::printPerformanceStatistics()
{
    SimulationSupervisor::printPerformanceStatistics();
    Print::printf(Print::INFO, "MicroenvironmentSupervisor working on replicate %d/%d and timestep %d/%d.", currentReplicateIndex+1, numberReplicates, currentTimestep, numberTimesteps);
//    Print::printf(Print::INFO, "  Performed %lld timesteps in %0.3e seconds (%0.3e timesteps/second), average PDE: %0.4e s, ME: %0.4e s, reconcile: %0.4e s.",stats_timesteps,convertHrToSeconds(stats_timestepTotalTime),double(stats_timesteps)/convertHrToSeconds(stats_timestepTotalTime), convertHrToSeconds(stats_timestepPDETime)/double(stats_timesteps), convertHrToSeconds(stats_timestepMETime)/double(stats_timesteps), convertHrToSeconds(stats_timestepReconcileTime)/double(stats_timesteps));
//    Print::printf(Print::INFO, "  ME solvers performed %lld steps in %0.3e seconds (%0.3e steps/second).", stats_workUnitsSteps, stats_workUnitTime, double(stats_workUnitsSteps)/stats_workUnitTime);
//    Print::printf(Print::INFO, "  PDE solvers performed %lld steps in %0.3e seconds (%0.3e steps/second).", stats_pdeWorkUnitsSteps, stats_pdeWorkUnitsTime, double(stats_pdeWorkUnitsSteps)/stats_pdeWorkUnitsTime);
}

}
}
