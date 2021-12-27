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
 * Author(s): Elijah Roberts
 */

#ifndef LM_RDME_RDMESOLVER_H_
#define LM_RDME_RDMESOLVER_H_

#include "lm/cme/CMESolver.h"
#include "lm/types/BoundaryConditions.pb.h"
#include "lm/input/DiffusionModel.pb.h"
#include "lm/rdme/DiffusionModel.h"
#include "lm/rdme/Lattice.h"
#include "robertslab/pbuf/NDArray.pb.h"

using lm::cme::CMESolver;
using lm::rdme::Lattice;

namespace lm {
namespace rdme {

class RDMESolver : public CMESolver
{
public:
    RDMESolver(RandomGenerator::Distributions neededDists);
    virtual ~RDMESolver();
    virtual bool needsDiffusionModel();
    virtual void setDiffusionModel(const lm::input::DiffusionModel& dm);
    virtual void reset();
    virtual void setOutputOptions(const lm::input::OutputOptions& outputOptions);
    virtual void getState(lm::io::TrajectoryState* state, uint trajectoryNumber=0);
    virtual void setState(const lm::io::TrajectoryState& state, uint trajectoryNumber=0);

protected:
    virtual void trajectoryStarted();
    virtual void timeUpdated();
    virtual void trajectoryFinished();

protected:

    // The diffusion model.
    DiffusionModel* diffusionModel;

    // The current state.
    Lattice* lattice;

    // Output options.
    bool writeLatticeTimeSeries;
    double latticeWriteInterval;
    double nextLatticeWriteTime;
    vector<double> latticeTimeSeriesTimes;
    vector<robertslab::pbuf::NDArray*> latticeTimeSeriesLattices;
};

}
}

#endif
