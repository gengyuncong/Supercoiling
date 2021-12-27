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

#ifndef LM_RDME_NEXTSUBVOLUMESOLVER_H_
#define LM_RDME_NEXTSUBVOLUMESOLVER_H_

#include "lm/io/TrajectoryState.pb.h"
#include "lm/rdme/Lattice.h"
#include "lm/rdme/RDMESolver.h"
#include "lm/reaction/ReactionQueue.h"
#include "lm/types/BoundaryConditions.pb.h"

using lm::reaction::ReactionQueue;

namespace lm {
namespace rdme {

class NextSubvolumeSolver : public RDMESolver
{
public:
    static bool registered;
    static bool registerClass();
    static void* allocateObject();

private:
    static const int NUM_NEIGHBORS=6;


public:
    NextSubvolumeSolver();
    NextSubvolumeSolver(RandomGenerator::Distributions neededDists);
    virtual ~NextSubvolumeSolver();
    virtual void setDiffusionModel(const lm::input::DiffusionModel& dm);
    virtual void reset();
    virtual void getState(lm::io::TrajectoryState* state, uint trajectoryNumber=0);
    virtual void setState(const lm::io::TrajectoryState& state, uint trajectoryNumber=0);
    virtual uint64_t generateTrajectory(uint64_t maxSteps);

private:
    void allocateRngBuffers();
    void deallocateRngBuffers();
    void checkSpeciesCountsAgainstLattice();

    void loadSubvolumeSpeciesCountsFromLattice(lattice_size_t subvolume);
    void saveSubvolumeSpeciesCountsToLattice(lattice_size_t subvolume);

    void updateAllSubvolumePropensities();
    void updateSubvolumePropensity(lattice_size_t subvolume);
    double calculateSubvolumePropensity(lattice_size_t s);
    double calculateSubvolumeDiffusionPropensity(lattice_size_t subvolume, site_t sourceSite);
    double calculateSubvolumeInfluxPropensity(lattice_size_t subvolume);

    bool performSubvolumeEvent(lattice_size_t s, bool& affectedNeighbor, lattice_size_t& neighborSubvolume);
    void performReactionEventInCurrentSubvolume(uint r);
    tuple<bool> performSubvolumeDiffusionEvent(lattice_size_t subvolume, site_t sourceSite, bool& affectedNeighbor, lattice_size_t& neighborSubvolume, double& rngValue);
    tuple<bool> performSubvolumeInfluxEvent(lattice_size_t subvolume, double& rngValue);
    void addParticles(lattice_size_t subvolume, particle_t particle, uint count);

private:
    double* expRngValues;
    double* uniRngValues;
    size_t nextExpRngValue;
    size_t nextUniRngValue;
    lattice_size_t numberSubvolumes;
    double latticeSpacingSquared;
    lm::types::BoundaryConditions::BoundaryConditionsType bc[NUM_NEIGHBORS];
    ReactionQueue* reactionQueue;
    int* currentSubvolumeSpeciesCounts;
};

}
}

#endif
