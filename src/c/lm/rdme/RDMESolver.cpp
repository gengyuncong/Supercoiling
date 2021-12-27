/*
 * University of Illinois Open Source License
 * Copyright 2008-2011 Luthey-Schulten Group,
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
 * Author(s): Elijah Roberts
 */

#include "lm/Exceptions.h"
#include "lm/Tune.h"
#include "lm/Print.h"
#include "lm/cme/CMESolver.h"
#include "lm/types/BoundaryConditions.pb.h"
#include "lm/input/DiffusionModel.pb.h"
#include "lm/types/Lattice.pb.h"
#include "lm/me/PropensityFunction.h"
#include "lm/rdme/Lattice.h"
#include "lm/rdme/ByteLattice.h"
#include "lm/rdme/DiffusionModel.h"
#include "lm/rdme/RDMESolver.h"
#include "lm/rng/RandomGenerator.h"
#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using lm::input::DiffusionModel;
using lm::rdme::Lattice;
using lm::rng::RandomGenerator;
using robertslab::pbuf::NDArraySerializer;

namespace lm {
namespace rdme {

RDMESolver::RDMESolver(RandomGenerator::Distributions neededDists)
:CMESolver(neededDists),diffusionModel(NULL),lattice(NULL)
{
}

RDMESolver::~RDMESolver()
{
    // Free any model memory.
    if (diffusionModel != NULL) delete diffusionModel; diffusionModel = NULL;

    // Free any memory associated with the state.
    if (lattice != NULL) delete lattice; lattice = NULL;
}

bool RDMESolver::needsDiffusionModel()
{
    return true;
}

void RDMESolver::setDiffusionModel(const lm::input::DiffusionModel& dm)
{
    CMESolver::setDiffusionModel(dm);

    // Validate the model.
    if (dm.number_species() != (int)reactionModel->numberSpecies) throw InvalidArgException("dm.number_species", "number of species in the diffusion model does not agree with the number in the reaction model");
    if (dm.number_reactions() != (int)reactionModel->numberReactions) throw InvalidArgException("dm.number_reactions", "number of reactions in the diffusion model does not agree with the number in the reaction model");
    if (dm.diffusion_matrix_size() != (dm.number_site_types()*dm.number_site_types()*dm.number_species())) throw InvalidArgException("dm", "diffusion matrix size does not agree with the number of species and site types");
    if (dm.reaction_location_matrix_size() != (dm.number_reactions() *dm.number_site_types())) throw InvalidArgException("dm", "reaction location matrix size does not agree with the number of reactions and site types");

    // Create the new model.
    if (diffusionModel != NULL) delete diffusionModel;
    diffusionModel = new DiffusionModel(dm);

    // Create the new lattice.
    if (lattice != NULL) delete lattice;
    lattice = new ByteLattice(diffusionModel->latticeXSize, diffusionModel->latticeYSize, diffusionModel->latticeZSize, diffusionModel->latticeSpacing, diffusionModel->particlesPerSite);

    // Update the propensity functions with the subvolume size.
    uint numberSubvolumes = lattice->getNumberSites();
    for (uint i=0; i<reactionModel->numberReactions; i++)
    {
        reactionModel->propensityFunctions[i]->changeVolume(1/double(numberSubvolumes));
//        if (reactionModel->propensityFunctions[i]->getOrder() >= 3)
//        {
//            Print::printf(Print::WARNING, "Propensity function of type %d was order %d, please ensure that this is a valid function for an RDME simulation.", reactionModel->propensityFunctions[i]->getType(),reactionModel->propensityFunctions[i]->getOrder());
//        }
    }
}

void RDMESolver::reset()
{
    if (diffusionModel == NULL || lattice == NULL) throw Exception("RDMESolver reset called before diffusion model was set.");

    CMESolver::reset();

    // Free any previous state.
    lattice->removeAllParticles();
    nextLatticeWriteTime = std::numeric_limits<double>::infinity();
    latticeTimeSeriesTimes.clear();
    latticeTimeSeriesLattices.clear();
}

void RDMESolver::setOutputOptions(const lm::input::OutputOptions& outputOptions)
{
    CMESolver::setOutputOptions(outputOptions);

    writeLatticeTimeSeries = outputOptions.has_lattice_write_interval();
    latticeWriteInterval = writeLatticeTimeSeries?outputOptions.lattice_write_interval():std::numeric_limits<double>::infinity();
}

void RDMESolver::getState(lm::io::TrajectoryState* state, uint trajectoryNumber)
{
    if (diffusionModel == NULL || lattice == NULL) throw Exception("RDMESolver get state called before diffusion model was set.");

    CMESolver::getState(state, trajectoryNumber);

    // Get the lattice sites.
    ndarray<uint8_t> sites(utuple(lattice->getSize().x,lattice->getSize().y,lattice->getSize().z), 0, ndarray_ArrayOrder::IMPL_ORDER);
    lattice->copySitesTo(&sites);
    state->mutable_rdme_state()->mutable_lattice()->set_allocated_sites(NDArraySerializer::serializeAllocate(sites));

    // Get the lattice particles.
    ndarray<uint8_t> particles(utuple(lattice->getSize().x,lattice->getSize().y,lattice->getSize().z,lattice->getMaxOccupancy()), 0, ndarray_ArrayOrder::IMPL_ORDER);
    lattice->copyParticlesTo(&particles);
    state->mutable_rdme_state()->mutable_lattice()->set_allocated_particles(NDArraySerializer::serializeAllocate(particles));
}

void RDMESolver::setState(const lm::io::TrajectoryState& state, uint trajectoryNumber)
{
    // Valdiate the state.
    if (diffusionModel == NULL || lattice == NULL) throw Exception("RDMESolver set state called before diffusion model was set.");
    if (!state.has_rdme_state()) throw Exception("State object does not contain rdme state to initialize the RDMESolver.");
    if (!state.rdme_state().lattice().has_sites()) throw Exception("State object does not contain the lattice sites to initialize the RDMESolver.");
    if (state.rdme_state().lattice().particles().shape(0) != lattice->getSize().x) throw Exception("State object and lattice have differing lattice x size",state.rdme_state().lattice().particles().shape(0),lattice->getSize().x);
    if (state.rdme_state().lattice().particles().shape(1) != lattice->getSize().y) throw Exception("State object and lattice have differing lattice y size",state.rdme_state().lattice().particles().shape(1),lattice->getSize().y);
    if (state.rdme_state().lattice().particles().shape(2) != lattice->getSize().z) throw Exception("State object and lattice have differing lattice z size",state.rdme_state().lattice().particles().shape(2),lattice->getSize().z);
    if (state.rdme_state().lattice().particles().shape(3) != lattice->getMaxOccupancy()) throw Exception("State object and lattice have differing number of particles per site",state.rdme_state().lattice().particles().shape(3),lattice->getMaxOccupancy());
    if (state.rdme_state().lattice().sites().shape(0) != lattice->getSize().x) throw Exception("State object and lattice have differing lattice sites x size",state.rdme_state().lattice().sites().shape(0),lattice->getSize().x);
    if (state.rdme_state().lattice().sites().shape(1) != lattice->getSize().y) throw Exception("State object and lattice have differing lattice sites y size",state.rdme_state().lattice().sites().shape(1),lattice->getSize().y);
    if (state.rdme_state().lattice().sites().shape(2) != lattice->getSize().z) throw Exception("State object and lattice have differing lattice sites z size",state.rdme_state().lattice().sites().shape(2),lattice->getSize().z);

    CMESolver::setState(state, trajectoryNumber);

    // Set the lattice sites.
    ndarray<uint8_t>* sites = NDArraySerializer::deserializeAllocate<uint8_t>(state.rdme_state().lattice().sites());
    lattice->copySitesFrom(sites);
    delete sites;

    // Set the lattice particles.
    ndarray<uint8_t>* particles = NDArraySerializer::deserializeAllocate<uint8_t>(state.rdme_state().lattice().particles());
    lattice->copyParticlesFrom(particles);
    delete particles;
}

void RDMESolver::trajectoryStarted()
{
    CMESolver::trajectoryStarted();

    //
    // Initialize any time series data.
    //

    // Get the interval for writing lattices.
    if (writeLatticeTimeSeries)
    {
        // Fill in the first set of values, if necessary.
        if (!previouslyStarted && writeInitialTrajectoryState)
        {
            PROF_BEGIN(PROF_NSM_SERIALIZE_LATTICE);

            // Get the lattice particles.
            ndarray<uint8_t> particles(utuple(lattice->getSize().x,lattice->getSize().y,lattice->getSize().z,lattice->getMaxOccupancy()));
            lattice->copyParticlesTo(&particles);
            latticeTimeSeriesLattices.push_back(robertslab::pbuf::NDArraySerializer::serializeAllocate(particles));
            PROF_END(PROF_NSM_SERIALIZE_LATTICE);
            latticeTimeSeriesTimes.push_back(time);
        }

        // Set the next write time.
        nextLatticeWriteTime = roundNextMultiple(time, latticeWriteInterval);
    }
}

void RDMESolver::timeUpdated()
{
    CMESolver::timeUpdated();

    //
    // Record any time series data.
    //

    // If we are writing lattice time steps, write out any lattice time steps before this event occurred.
    if (writeLatticeTimeSeries)
    {
        // Write time steps until the next write time is past the current time.
        while (nextLatticeWriteTime <= (time+EPS))
        {
            PROF_BEGIN(PROF_NSM_SERIALIZE_LATTICE);
            ndarray<uint8_t> particles(utuple(lattice->getSize().x,lattice->getSize().y,lattice->getSize().z,lattice->getMaxOccupancy()));
            lattice->copyParticlesTo(&particles);
            latticeTimeSeriesLattices.push_back(robertslab::pbuf::NDArraySerializer::serializeAllocate(particles));
            PROF_END(PROF_NSM_SERIALIZE_LATTICE);
            latticeTimeSeriesTimes.push_back(nextLatticeWriteTime);
            nextLatticeWriteTime += latticeWriteInterval;
        }
    }
}
void RDMESolver::trajectoryFinished()
{
    CMESolver::trajectoryFinished();

    //
    // Record any final time series data.
    //

    // If we hit any limit, write out the final time if requested.
    if (status == lm::message::WorkUnitStatus::LIMIT_REACHED && writeFinalTrajectoryState)
    {
        // write out the last lattice, or ensure that it already has been
        if (writeLatticeTimeSeries && (latticeTimeSeriesTimes.empty() || latticeTimeSeriesTimes.back() + EPS < time))
        {
            PROF_BEGIN(PROF_NSM_SERIALIZE_LATTICE);
            ndarray<uint8_t> particles(utuple(lattice->getSize().x,lattice->getSize().y,lattice->getSize().z,lattice->getMaxOccupancy()));
            lattice->copyParticlesTo(&particles);
            latticeTimeSeriesLattices.push_back(robertslab::pbuf::NDArraySerializer::serializeAllocate(particles));
            PROF_END(PROF_NSM_SERIALIZE_LATTICE);
            latticeTimeSeriesTimes.push_back(time);
        }
    }

    //
    // Add any saved records to the output.
    //

    // If we have any lattice time series data, add them to the output message.
    if (latticeTimeSeriesTimes.size() > 0 || latticeTimeSeriesLattices.size() > 0)
    {
        // Mark that the message does contain some data.
        output->set_has_output(true);

        // Make sure the arrays are of a consistent size.
        if (latticeTimeSeriesLattices.size() == latticeTimeSeriesTimes.size())
        {
            lm::io::LatticeTimeSeries* latticeTimeSeriesDataSet = output->mutable_lattice_time_series();
            latticeTimeSeriesDataSet->set_trajectory_id(trajectoryId);

            // Serialize the times.
            robertslab::pbuf::NDArraySerializer::serializeInto<double>(latticeTimeSeriesDataSet->mutable_times(), latticeTimeSeriesTimes.data(), utuple(static_cast<uint>(latticeTimeSeriesTimes.size())));

            // Serialize the lattice counts.
            for (uint i=0; i<latticeTimeSeriesLattices.size(); i++)
                latticeTimeSeriesDataSet->add_lattices()->set_allocated_particles(latticeTimeSeriesLattices[i]);
        }
        else
        {
            Print::printf(Print::ERROR, "Lattice time series counts and time mismatch %d,%d", latticeTimeSeriesLattices.size(), latticeTimeSeriesTimes.size());
        }
    }
}

}
}
