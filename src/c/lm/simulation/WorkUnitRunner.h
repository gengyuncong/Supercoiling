/*
 * University of Illinois Open Source License
 * Copyright 2011 Luthey-Schulten Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
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
 * Urbana-Champaign, nor the names of its contributors may be used to endorse or
 * promote products derived from this Software without specific prior written
 * permission.
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

#ifndef LM_MAIN_WORKUNITRUNNER_H_
#define LM_MAIN_WORKUNITRUNNER_H_

#include <string>
#include <vector>
#include "lm/main/Solver.h"
#include "lm/me/MESolver.h"
#include "lm/message/Communicator.h"
#include "lm/message/Endpoint.pb.h"
#include "lm/message/RunWorkUnit.pb.h"
#include "lm/message/StartWorkUnitRunner.pb.h"
#include "lm/pde/DiffusionPDESolver.h"
#include "lm/thread/Thread.h"
#include "lm/thread/Worker.h"
#include "lm/Types.h"

namespace lm {
namespace main {

class WorkUnitRunner : public lm::thread::Worker
{
public:
    WorkUnitRunner(const lm::message::StartWorkUnitRunner& properties);
    virtual ~WorkUnitRunner();
    virtual void wake();
    virtual int run();
    virtual void runWorkUnits(const lm::message::RunWorkUnit& msg);

protected:
//    virtual bool flushOutputCondensed();

protected:
    int id;
    uint64_t simulationPhaseID;

//    // data that gets cached for later use with condensed output
//    lm::message::Endpoint condensedOutputAddress;
//    std::vector<uint64_t> firstTrajectoryIds;
//    std::vector<uint64_t> lastTrajectoryIds;

    lm::message::Communicator* communicator;
    lm::message::StartWorkUnitRunner properties;
    lm::me::MESolver* meSolver;
    lm::pde::DiffusionPDESolver* pdeSolver;
};

}
}

#endif
