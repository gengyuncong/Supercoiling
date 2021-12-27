/*
 * University of Illinois Open Source License
 * Copyright 2008-2012 Luthey-Schulten Group,
 * Copyright 2012-2019 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
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
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
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

#include "lm/ClassFactory.h"
#include "lm/Exceptions.h"
#include "lm/Tune.h"
#include "lm/Print.h"
#include "lm/cme/CMESolver.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/types/Lattice.pb.h"
#include "lm/io/LatticeTimeSeries.pb.h"
#include "lm/io/SpeciesCounts.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/message/WorkUnitOutput.pb.h"
#include "lm/rdme/Lattice.h"
#include "lm/rdme/ByteLattice.h"
#include "lm/rdme/NextSubvolumeSolver.h"
#include "lm/reaction/ReactionQueue.h"
#include "lm/rng/RandomGenerator.h"
#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using lm::rng::RandomGenerator;

namespace lm {
namespace rdme {

bool NextSubvolumeSolver::registered=NextSubvolumeSolver::registerClass();

bool NextSubvolumeSolver::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::me::MESolver","lm::rdme::NextSubvolumeSolver",&NextSubvolumeSolver::allocateObject);
    return true;
}

void* NextSubvolumeSolver::allocateObject()
{
    return new NextSubvolumeSolver();
}


NextSubvolumeSolver::NextSubvolumeSolver()
:RDMESolver((RandomGenerator::Distributions)(RandomGenerator::EXPONENTIAL|RandomGenerator::UNIFORM)),
expRngValues(NULL),uniRngValues(NULL),nextExpRngValue(0),nextUniRngValue(0),numberSubvolumes(0),latticeSpacingSquared(0.0),reactionQueue(NULL),currentSubvolumeSpeciesCounts(NULL)
{
    // Allocate the RNG buffers.
    allocateRngBuffers();
}

NextSubvolumeSolver::NextSubvolumeSolver(RandomGenerator::Distributions neededDists)
:RDMESolver((RandomGenerator::Distributions)(RandomGenerator::EXPONENTIAL|RandomGenerator::UNIFORM|neededDists)),
expRngValues(NULL),uniRngValues(NULL),nextExpRngValue(0),nextUniRngValue(0),numberSubvolumes(0),latticeSpacingSquared(0.0),reactionQueue(NULL),currentSubvolumeSpeciesCounts(NULL)
{
    // Allocate the RNG buffers.
    allocateRngBuffers();
}

NextSubvolumeSolver::~NextSubvolumeSolver()
{
    if (currentSubvolumeSpeciesCounts != NULL) delete[] currentSubvolumeSpeciesCounts; currentSubvolumeSpeciesCounts = NULL;
    if (reactionQueue != NULL) delete reactionQueue; reactionQueue = NULL;
    deallocateRngBuffers();
}

void NextSubvolumeSolver::allocateRngBuffers()
{
    if (expRngValues == NULL || uniRngValues == NULL)
    {
        expRngValues = new double[TUNE_LOCAL_RNG_CACHE_SIZE];
        uniRngValues = new double[TUNE_LOCAL_RNG_CACHE_SIZE];
        nextExpRngValue = TUNE_LOCAL_RNG_CACHE_SIZE;
        nextUniRngValue = TUNE_LOCAL_RNG_CACHE_SIZE;
    }
}

void NextSubvolumeSolver::deallocateRngBuffers()
{
    if (expRngValues != NULL) delete[] expRngValues; expRngValues = NULL;
    if (uniRngValues != NULL) delete[] uniRngValues; uniRngValues = NULL;
    nextExpRngValue = 0;
    nextUniRngValue = 0;
}

void NextSubvolumeSolver::setDiffusionModel(const lm::input::DiffusionModel& dm)
{
    RDMESolver::setDiffusionModel(dm);

    // Fill in some parameter variables.
    numberSubvolumes = lattice->getNumberSites();
    latticeSpacingSquared = diffusionModel->latticeSpacing*diffusionModel->latticeSpacing;

    // Fill in the boundary conditions.
    if (diffusionModel->boundaryConditions.axis_specific_boundaries())
    {
        bc[0]=diffusionModel->boundaryConditions.x_minus();
        bc[1]=diffusionModel->boundaryConditions.x_plus();
        bc[2]=diffusionModel->boundaryConditions.y_minus();
        bc[3]=diffusionModel->boundaryConditions.y_plus();
        bc[4]=diffusionModel->boundaryConditions.z_minus();
        bc[5]=diffusionModel->boundaryConditions.z_plus();
    }
    else
    {
        for (int j=0; j<NUM_NEIGHBORS; j++)
            bc[j]=diffusionModel->boundaryConditions.global();
    }

    // Allocate the subvolume species counts.
    if (currentSubvolumeSpeciesCounts != NULL) delete[] currentSubvolumeSpeciesCounts;
    currentSubvolumeSpeciesCounts = new int[reactionModel->numberSpecies];
}

void NextSubvolumeSolver::reset()
{
    RDMESolver::reset();

    // Reset the subvolume species counts.
    for (uint i=0; i<reactionModel->numberSpecies; i++)
        currentSubvolumeSpeciesCounts[i] = 0;

    // Create the reaction queue.
    PROF_BEGIN(PROF_NSM_INIT_QUEUE);
    if (reactionQueue != NULL) delete reactionQueue;
    reactionQueue = new ReactionQueue(numberSubvolumes);
    PROF_END(PROF_NSM_INIT_QUEUE);
}

void NextSubvolumeSolver::getState(lm::io::TrajectoryState* state, uint trajectoryNumber)
{
    RDMESolver::getState(state, trajectoryNumber);
}

void NextSubvolumeSolver::setState(const lm::io::TrajectoryState& state, uint trajectoryNumber)
{
    RDMESolver::setState(state, trajectoryNumber);

    // Initialize the reaction queue.
    updateAllSubvolumePropensities();
}

uint64_t NextSubvolumeSolver::generateTrajectory(uint64_t maxSteps)
{
    if (reactionModel == NULL) throw Exception("NextSubvolumeSolver did not have a reaction model.");
    if (diffusionModel == NULL) throw Exception("NextSubvolumeSolver did not have a diffusion model.");
    if (reactionQueue == NULL) throw Exception("NextSubvolumeSolver state was not initialized.");

    // Make sure we have propensity functions for every reaction.
    for (uint i=0; i<reactionModel->numberReactions; i++)
        if (reactionModel->propensityFunctions[i] == NULL)
            throw Exception("A reaction did not have a valid propensity function",i);

    // Make sure that the initial species counts agree with the actual number in the lattice.
    checkSpeciesCountsAgainstLattice();

    // Call the trajectory started method.
    trajectoryStarted();

    // Run the next subvolume method.
    Print::printf(Print::DEBUG, "Running next subvolume simulation for %d steps with %d species, %d reactions, %d subvolumes, %d site types, %d limits.", maxSteps, reactionModel->numberSpecies, reactionModel->numberReactions, numberSubvolumes, diffusionModel->numberSiteTypes, timeLimits.limits_size()+stateLimits.limits_size());
    PROF_BEGIN(PROF_SIM_EXECUTE);
    uint64_t steps=0;
    while (true)
    {
        // See if we have finished the steps.
        if (steps >= maxSteps)
        {
            status = lm::message::WorkUnitStatus::STEPS_FINISHED;
            break;
        }

        // Increment the steps.
        steps++;

        // Get the next subvolume with a reaction and the reaction time.
        lattice_size_t subvolume = reactionQueue->getNextReaction();
        double newTime = reactionQueue->getReactionEvent(subvolume).time;

        // If the the time until the next reaction is infinite, stop the simulation.
        if (newTime == std::numeric_limits<double>::infinity())
        {
            setTrajectoryToMaxTime();
            break;
        }

        // Update the time. If we are outside the time limits, stop the trajectory.
        if (performTimeIncrement(newTime-time))
            break;

        // Update the system with the reaction. If we are outside of the state limits, stop the trajectory.
        PROF_BEGIN(PROF_NSM_PERFORM_SUBVOLUME_EVENT);
        bool affectedNeighbor = false;
        lattice_size_t neighborSubvolume;
        if (performSubvolumeEvent(subvolume, affectedNeighbor, neighborSubvolume))
           break;
        PROF_END(PROF_NSM_PERFORM_SUBVOLUME_EVENT);

       // Update the propensity in the affected subvolumes.
       PROF_BEGIN(PROF_NSM_UPDATE_SUBVOLUME_PROPENSITY);
       updateSubvolumePropensity(subvolume);
       if (affectedNeighbor) updateSubvolumePropensity(neighborSubvolume);
       PROF_END(PROF_NSM_UPDATE_SUBVOLUME_PROPENSITY);
    }
    PROF_END(PROF_SIM_EXECUTE);
    Print::printf(Print::DEBUG, "Generated trajectory through time %e with %d steps and final status %d.", time, steps, status);

    // Make sure that the final species counts agree with the actual number in the lattice.
    checkSpeciesCountsAgainstLattice();

    // Call the trajectory finished method.
    trajectoryFinished();

    return steps;
}

void NextSubvolumeSolver::checkSpeciesCountsAgainstLattice()
{
    std::map<particle_t,uint> particleCounts = lattice->getParticleCounts();
    for (uint i=0; i<reactionModel->numberSpecies; i++)
    {
        if (speciesCounts[i] != int((particleCounts.count(i)>0)?particleCounts[i]:0))
            THROW_EXCEPTION(lm::RuntimeException, "Consistency error between species counts and lattice data: %d %d %d", i, speciesCounts[i], ((particleCounts.count(i)>0)?particleCounts[i]:0));
    }
}

void NextSubvolumeSolver::updateAllSubvolumePropensities()
{
    // Update all of the subvolumes.
    PROF_BEGIN(PROF_NSM_BUILD_QUEUE);
    for (lattice_size_t s=0; s<numberSubvolumes; s++)
    {
        double propensity = calculateSubvolumePropensity(s);
        double newTime = std::numeric_limits<double>::infinity();
        if (propensity > 0.0)
        {
            if (nextExpRngValue >= TUNE_LOCAL_RNG_CACHE_SIZE)
            {
                rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
                nextExpRngValue=0;
            }
            newTime = time+expRngValues[nextExpRngValue++]/propensity;
        }
        reactionQueue->updateReactionEvent(s, newTime, propensity);
    }
    PROF_END(PROF_NSM_BUILD_QUEUE);
}

void NextSubvolumeSolver::updateSubvolumePropensity(lattice_size_t subvolume)
{
    PROF_BEGIN(PROF_NSM_CALCULATE_SUBVOLUME_PROPENSITY);
    double propensity = calculateSubvolumePropensity(subvolume);
    PROF_END(PROF_NSM_CALCULATE_SUBVOLUME_PROPENSITY);
    double newTime = std::numeric_limits<double>::infinity();
    if (propensity > 0.0)
    {
        if (nextExpRngValue >= TUNE_LOCAL_RNG_CACHE_SIZE)
        {
            rng->getExpRandomDoubles(expRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
            nextExpRngValue=0;
        }
        newTime = time+expRngValues[nextExpRngValue++]/propensity;
    }
    PROF_BEGIN(PROF_NSM_UPDATE_QUEUE);
    reactionQueue->updateReactionEvent(subvolume, newTime, propensity);
    PROF_END(PROF_NSM_UPDATE_QUEUE);
}

double NextSubvolumeSolver::calculateSubvolumePropensity(lattice_size_t subvolume)
{
    // Update the species counts member for this subvolume.
    loadSubvolumeSpeciesCountsFromLattice(subvolume);

    // Get the site type for this subvolume.
    site_t sourceSite = lattice->getSiteType(subvolume);

    // Calculate all of the reaction propensities.
    double subvolumePropensity = 0.0;
    for (uint i=0; i<reactionModel->numberReactions; i++)
    {
        // Make sure the reaction can occur in this subvolume.
        if (diffusionModel->RL[i*diffusionModel->numberSiteTypes+sourceSite])
        {
            subvolumePropensity += reactionModel->propensityFunctions[i]->calculate(time, currentSubvolumeSpeciesCounts, reactionModel->numberSpecies);
        }
    }

    // Add in the diffusion propensity.
    subvolumePropensity+=calculateSubvolumeDiffusionPropensity(subvolume, sourceSite);

    // Add in any propensity for particle influx from the boundaries.
    subvolumePropensity+=calculateSubvolumeInfluxPropensity(subvolume);

    return subvolumePropensity;
}

double NextSubvolumeSolver::calculateSubvolumeDiffusionPropensity(lattice_size_t subvolume, site_t sourceSite)
{
    double subvolumePropensity=0.0;

    // Get the neighboring sites.
    lattice_size_t neighboringSubvolumes[NUM_NEIGHBORS];
    lattice->getNeighboringSites(subvolume, neighboringSubvolumes, false);

    // Calculate all of the diffusion propensities.
    for (uint i=0; i<reactionModel->numberSpecies; i++)
    {
        if (currentSubvolumeSpeciesCounts[i] > 0)
        {
            for (int j=0; j<NUM_NEIGHBORS; j++)
            {
                // See if the neighbor is a boundary.
                lattice_size_t neighborIndex=neighboringSubvolumes[j];
                if (neighborIndex == OUTSIDE_LATTICE_BOUNDARIES)
                {
                    if (bc[j] == lm::types::BoundaryConditions::REFLECTING)
                    {
                    }
                    else if (bc[j] == lm::types::BoundaryConditions::ABSORBING || bc[j] == lm::types::BoundaryConditions::FIXED_CONCENTRATION || bc[j] == lm::types::BoundaryConditions::FIXED_GRADIENT)
                    {
                        subvolumePropensity += (double(currentSubvolumeSpeciesCounts[i])) * (diffusionModel->DF[sourceSite*diffusionModel->numberSiteTypes*reactionModel->numberSpecies + sourceSite*reactionModel->numberSpecies + i]/latticeSpacingSquared);
                    }
                    else if (bc[j] == lm::types::BoundaryConditions::PERIODIC)
                    {
                        lattice_size_t neighboringSubvolumesPeriodic[NUM_NEIGHBORS];
                        lattice->getNeighboringSites(subvolume, neighboringSubvolumesPeriodic, true);
                        lattice_size_t neighborIndexPeriodic=neighboringSubvolumesPeriodic[j];
                        subvolumePropensity += (double(currentSubvolumeSpeciesCounts[i])) * (diffusionModel->DF[sourceSite*diffusionModel->numberSiteTypes*reactionModel->numberSpecies + lattice->getSiteType(neighborIndexPeriodic)*reactionModel->numberSpecies + i]/latticeSpacingSquared);
                    }
                }
                else
                {
                    subvolumePropensity += (double(currentSubvolumeSpeciesCounts[i])) * (diffusionModel->DF[sourceSite*diffusionModel->numberSiteTypes*reactionModel->numberSpecies + lattice->getSiteType(neighborIndex)*reactionModel->numberSpecies + i]/latticeSpacingSquared);
                }
            }
        }
    }
    return subvolumePropensity;
}

double NextSubvolumeSolver::calculateSubvolumeInfluxPropensity(lattice_size_t subvolume)
{
    if (diffusionModel->hasBoundaryInflux && lattice->isBoundarySite(subvolume))
    {
        return diffusionModel->boundaryInflux[subvolume];
    }

    return 0.0;
}

void  NextSubvolumeSolver::loadSubvolumeSpeciesCountsFromLattice(lattice_size_t subvolume)
{
    // Reset the species counts.
    for (uint i=0; i<reactionModel->numberSpecies; i++)
        currentSubvolumeSpeciesCounts[i] = 0;

    // Count the species that are in this subvolume.
    site_size_t numberParticles = lattice->getOccupancy(subvolume);
    for (site_size_t i=0; i<numberParticles; i++)
        currentSubvolumeSpeciesCounts[lattice->getParticle(subvolume, i)]++;
}

void NextSubvolumeSolver::saveSubvolumeSpeciesCountsToLattice(lattice_size_t subvolume)
{
    lattice->removeParticles(subvolume);
    for (uint i=0; i<reactionModel->numberSpecies; i++)
            addParticles(subvolume, i, uint(currentSubvolumeSpeciesCounts[i]));
}

bool NextSubvolumeSolver::performSubvolumeEvent(lattice_size_t subvolume, bool& affectedNeighbor, lattice_size_t& neighborSubvolume)
{
    // Only set the affected neighbor if this was a diffusion event.
    affectedNeighbor = false;

    // Update the species counts member for this subvolume.
    loadSubvolumeSpeciesCountsFromLattice(subvolume);

    // Stretch the random value across the total propensity range.
    if (nextUniRngValue >= TUNE_LOCAL_RNG_CACHE_SIZE)
    {
        rng->getRandomDoubles(uniRngValues,TUNE_LOCAL_RNG_CACHE_SIZE);
        nextUniRngValue=0;
    }
    double rngValue = uniRngValues[nextUniRngValue++]*reactionQueue->getReactionEvent(subvolume).propensity;

    // Get the site type for this subvolume.
    site_t sourceSite = lattice->getSiteType(subvolume);

    // See if it was a reaction that occurred.
    for (uint r=0; r<reactionModel->numberReactions; r++)
    {
        // Make sure the reaction can occur in this subvolume.
        if (diffusionModel->RL[r*diffusionModel->numberSiteTypes+sourceSite])
        {
            double reactionPropensity = reactionModel->propensityFunctions[r]->calculate(time, currentSubvolumeSpeciesCounts, reactionModel->numberSpecies);
            if (reactionPropensity > 0.0)
            {
                if (rngValue <= reactionPropensity)
                {
                    bool outsideLimit = performReactionEvent(r);
                    performReactionEventInCurrentSubvolume(r);
                    saveSubvolumeSpeciesCountsToLattice(subvolume);
                    return outsideLimit;
                }
                else
                {
                    rngValue -= reactionPropensity;
                }
            }
        }
    }

    // See if it was a diffusion event that occurred.
    tuple<bool> outcome = performSubvolumeDiffusionEvent(subvolume, sourceSite, affectedNeighbor, neighborSubvolume, rngValue);
    if (outcome[0]) return outcome[1];

    // See if it was an influx event that occurred.
    outcome = performSubvolumeInfluxEvent(subvolume, rngValue);
    if (outcome[0]) return outcome[1];

    THROW_EXCEPTION(lm::RuntimeException, "Unable to determine correct reaction or diffusion event in subvolume.");
}

void NextSubvolumeSolver::performReactionEventInCurrentSubvolume(uint r)
{
    for (uint i=0; i<reactionModel->numberDependentSpecies[r]; i++)
    {
        currentSubvolumeSpeciesCounts[reactionModel->dependentSpecies[r][i]] += reactionModel->dependentSpeciesChange[r][i];
    }
}

tuple<bool> NextSubvolumeSolver::performSubvolumeDiffusionEvent(lattice_size_t subvolume, site_t sourceSite, bool& affectedNeighbor, lattice_size_t& neighborSubvolume, double& rngValue)
{
    // Get the neighboring sites.
    lattice_size_t neighboringSubvolumes[NUM_NEIGHBORS];
    lattice->getNeighboringSites(subvolume, neighboringSubvolumes, false);

    // See if it was a diffusion event that occurred.
    for (uint i=0; i<reactionModel->numberSpecies; i++)
    {
        if (currentSubvolumeSpeciesCounts[i] > 0)
        {
            for (int j=0; j<NUM_NEIGHBORS; j++)
            {
                // See if the neighbor is a boundary.
                lattice_size_t neighborIndex=neighboringSubvolumes[j];
                if (neighborIndex == OUTSIDE_LATTICE_BOUNDARIES)
                {
                    if (bc[j] == lm::types::BoundaryConditions::REFLECTING)
                    {
                    }
                    else if (bc[j] == lm::types::BoundaryConditions::ABSORBING ||bc[j] == lm::types::BoundaryConditions::FIXED_CONCENTRATION || bc[j] == lm::types::BoundaryConditions::FIXED_GRADIENT)
                    {
                        double diffusionPropensity = ((double)currentSubvolumeSpeciesCounts[i]) * (diffusionModel->DF[sourceSite*diffusionModel->numberSiteTypes*reactionModel->numberSpecies + sourceSite*reactionModel->numberSpecies + i]/latticeSpacingSquared);

                        // See if this is the diffusion event that occurred.
                        if (rngValue <= diffusionPropensity)
                        {
                            speciesCounts[i]--;
                            if (hasUpdateSpeciesCountsListeners) speciesCountsUpdated();
                            currentSubvolumeSpeciesCounts[i]--;
                            saveSubvolumeSpeciesCountsToLattice(subvolume);
                            affectedNeighbor = false;
                            return tuple<bool>(true,isTrajectoryOutsideStateLimits());
                        }
                        else
                        {
                            rngValue -= diffusionPropensity;
                        }
                    }
                    else if (bc[j] == lm::types::BoundaryConditions::PERIODIC)
                    {
                        lattice_size_t neighboringSubvolumesPeriodic[NUM_NEIGHBORS];
                        lattice->getNeighboringSites(subvolume, neighboringSubvolumesPeriodic, true);
                        lattice_size_t neighborIndexPeriodic=neighboringSubvolumesPeriodic[j];
                        double diffusionPropensity = (double(currentSubvolumeSpeciesCounts[i])) * (diffusionModel->DF[sourceSite*diffusionModel->numberSiteTypes*reactionModel->numberSpecies + lattice->getSiteType(neighborIndexPeriodic)*reactionModel->numberSpecies + i]/latticeSpacingSquared);

                        // See if this is the diffusion event that occurred.
                        if (rngValue <= diffusionPropensity)
                        {
                            currentSubvolumeSpeciesCounts[i]--;
                            saveSubvolumeSpeciesCountsToLattice(subvolume);
                            affectedNeighbor = true;
                            neighborSubvolume = neighborIndexPeriodic;
                            addParticles(neighborSubvolume, i, 1);
                            return tuple<bool>(true,false);
                        }
                        else
                        {
                            rngValue -= diffusionPropensity;
                        }
                    }
                }
                else
                {
                    double diffusionPropensity = (double(currentSubvolumeSpeciesCounts[i])) * (diffusionModel->DF[sourceSite*diffusionModel->numberSiteTypes*reactionModel->numberSpecies + lattice->getSiteType(neighborIndex)*reactionModel->numberSpecies + i]/latticeSpacingSquared);

                    // See if this is the diffusion event that occurred.
                    if (rngValue <= diffusionPropensity)
                    {
                        currentSubvolumeSpeciesCounts[i]--;
                        saveSubvolumeSpeciesCountsToLattice(subvolume);
                        affectedNeighbor = true;
                        neighborSubvolume = neighboringSubvolumes[j];
                        addParticles(neighborSubvolume, i, 1);
                        return tuple<bool>(true,false);
                    }
                    else
                    {
                        rngValue -= diffusionPropensity;
                    }
                }
            }
        }
    }
    return tuple<bool>(false,false);
}

tuple<bool> NextSubvolumeSolver::performSubvolumeInfluxEvent(lattice_size_t subvolume, double& rngValue)
{
    if (diffusionModel->hasBoundaryInflux && lattice->isBoundarySite(subvolume))
    {
        double influxPropensity = diffusionModel->boundaryInflux[subvolume];
        if (rngValue <= influxPropensity)
        {
            speciesCounts[diffusionModel->boundaryConditions.boundary_species()]++;
            if (hasUpdateSpeciesCountsListeners) speciesCountsUpdated();
            currentSubvolumeSpeciesCounts[diffusionModel->boundaryConditions.boundary_species()]++;
            saveSubvolumeSpeciesCountsToLattice(subvolume);
            return tuple<bool>(true,isTrajectoryOutsideStateLimits());
        }
        else
        {
            rngValue -= influxPropensity;
        }
    }
    return tuple<bool>(false,false);
}


void NextSubvolumeSolver::addParticles(lattice_size_t subvolume, particle_t particle, uint count)
{
    for (uint i=0; i<count; i++)
    {
        try
        {
            lattice->addParticle(subvolume, particle);
        }
        catch (InvalidParticleException& e)
        {
            // We need to perform some overflow processing.
            lattice_size_t neighboringSubvolumes[NUM_NEIGHBORS];
            lattice->getNeighboringSites(subvolume, neighboringSubvolumes, diffusionModel->boundaryConditions.global() == lm::types::BoundaryConditions::PERIODIC);
            bool handled = false;
            for (uint i=0; i<NUM_NEIGHBORS && !handled; i++)
            {
                lattice_size_t neighborIndex=neighboringSubvolumes[i];
                if (neighborIndex != OUTSIDE_LATTICE_BOUNDARIES)
                {
                    if (lattice->getOccupancy(neighborIndex) < lattice->getMaxOccupancy() && lattice->getSiteType(neighborIndex) == lattice->getSiteType(subvolume))
                    {
                        lattice->addParticle(neighborIndex, particle);
                        handled = true;
                        Print::printf(Print::WARNING, "Handled overflow of particle type %d (%d total) from subvolume %d (type %d,occupancy %d) by moving to subvolume %d (type %d,occupancy %d).", particle, count, subvolume, lattice->getSiteType(subvolume), lattice->getOccupancy(subvolume), neighboringSubvolumes[i], lattice->getSiteType(neighboringSubvolumes[i]), lattice->getOccupancy(neighboringSubvolumes[i]));
                    }
                }
            }
            if (!handled) THROW_EXCEPTION(lm::RuntimeException, "Unable to handle overflow at site %d", subvolume);
        }
    }
}

}
}
