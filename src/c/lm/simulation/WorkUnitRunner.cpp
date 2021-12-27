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

#include <string>
#include <map>
#include <pthread.h>
#include <vector>

#include "hrtime.h"
#include "lm/ClassFactory.h"
#include "lm/Exceptions.h"
#include "lm/cme/CMESolver.h"
#include "lm/Print.h"
#if defined(OPT_CUDA)
#include "lm/Cuda.h"
#endif
#include "lm/me/MESolver.h"
#include "lm/message/FinishedWorkUnit.pb.h"
#include "lm/message/Message.pb.h"
#include "lm/message/RunWorkUnit.pb.h"
#include "lm/message/StartedWorkUnit.pb.h"
#include "lm/message/StartWorkUnitRunner.pb.h"
#include "lm/message/WorkUnit.pb.h"
#include "lm/pde/DiffusionPDESolver.h"
#include "lm/simulation/WorkUnitRunner.h"
#include "lm/thread/Thread.h"
#include "lm/thread/Worker.h"
#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"

using std::map;
using std::string;
using std::vector;
using lm::message::Communicator;
using lm::message::Endpoint;

namespace lm {
namespace main {

WorkUnitRunner::WorkUnitRunner(const lm::message::StartWorkUnitRunner& msg)
:id(-1),simulationPhaseID(0),communicator(NULL),properties(msg),meSolver(NULL),pdeSolver(NULL)
{
    id = msg.work_unit_runner_id();

    // Create the communicator.
    communicator = lm::message::Communicator::createObjectOfDefaultSubclass(false);
}

WorkUnitRunner::~WorkUnitRunner()
{
    if (meSolver != NULL) delete meSolver; meSolver = NULL;
    if (pdeSolver != NULL) delete pdeSolver; pdeSolver = NULL;
    if (communicator != NULL) delete communicator; communicator = NULL;
}

void WorkUnitRunner::wake()
{
    lm::message::Message msg;
    msg.mutable_ping_target()->set_id(0);
    communicator->sendMessage(communicator->getSourceAddress(), &msg);
}

int WorkUnitRunner::run()
{
    try
    {
        Print::printf(Print::INFO, "Work Unit runner %s started with %d cpu cores (affinity=%d) and %d gpus using solver %s.", Communicator::printableAddress(communicator->getSourceAddress()).c_str(), properties.cpu_size(), properties.use_cpu_affinity(), properties.gpu_size(),properties.me_solver().c_str());

        // Set the processor affinity.
        if (properties.use_cpu_affinity() && properties.cpu_size() > 0)
        {
            setAffinity(properties.cpu(0));
        }

        // Set the GPU affinity.
        #if defined(OPT_CUDA)
        if (properties.gpu_size() > 0)
        {
            Print::printf(Print::INFO, "Work Unit runner %s using gpu device %d.", Communicator::printableAddress(communicator->getSourceAddress()).c_str(), properties.gpu(0));
            lm::CUDA::setCurrentDevice(properties.gpu(0));
        }
        #endif

        // Create the solvers.
        meSolver = static_cast<lm::me::MESolver*>(lm::ClassFactory::getInstance().allocateObjectOfClass("lm::me::MESolver",properties.me_solver()));
        pdeSolver = static_cast<lm::pde::DiffusionPDESolver*>(lm::ClassFactory::getInstance().allocateObjectOfClass("lm::pde::DiffusionPDESolver", properties.diffusion_pde_solver()));

        // Set the solver resources.
        vector<int> cpus;
        vector<int> gpus;
        for (int i=0; i<properties.cpu_size(); i++) cpus.push_back(properties.cpu(i));
        for (int i=0; i<properties.gpu_size(); i++) gpus.push_back(properties.gpu(i));
        meSolver->setComputeResources(cpus, gpus);
        pdeSolver->setComputeResources(cpus, gpus);

        // Tell the supervisor the runner has started.
        lm::message::Message msgp;
        lm::message::StartedWorkUnitRunner* msg = msgp.mutable_started_work_unit_runner();
        msg->set_work_unit_runner_id(id);
        msg->mutable_address()->CopyFrom(communicator->getSourceAddress());
        msg->set_simultaneous_work_units(meSolver->getSimultaneousTrajectories());
        communicator->sendMessage(communicator->getHypervisorAddress(), &msgp);

        // Loop reading messages.
        lm::message::Message message;
        while (true)
        {
            // Read the next message.
            communicator->receiveMessage(&message);

            // Do something with the message.
            if (message.has_run_work_unit())
            {
                runWorkUnits(message.run_work_unit());
            }
            else if (message.has_ping_target())
            {
                if (!running)
                {
//                    // Final flush of the condensed output.
//                    flushOutputCondensed();

                    // If we are done running, stop the loop.
                    break;
                }
            }
            else
            {
                Print::printf(Print::ERROR, "Work unit runner received an unknown message: {\n%s}",message.DebugString().c_str());
            }

            // Clear the message object so it can be used again.
            message.Clear();
        }

        // Delete the solver.
        if (meSolver != NULL) delete meSolver; meSolver = NULL;
        if (pdeSolver != NULL) delete pdeSolver; pdeSolver = NULL;

        Print::printf(Print::INFO, "Work unit runner %s finished.", Communicator::printableAddress(communicator->getSourceAddress()).c_str());
        return 0;
    }
    catch (lm::Exception& e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (std::exception& e)
    {
        Print::printf(Print::FATAL, "Exception during execution: %s (%s:%d)", e.what(), __FILE__, __LINE__);
    }
    catch (...)
    {
        Print::printf(Print::FATAL, "Unknown Exception during execution (%s:%d)", __FILE__, __LINE__);
    }
    exit(-1);
}

void WorkUnitRunner::runWorkUnits(const lm::message::RunWorkUnit& rwuMsg)
{
    PROF_BEGIN(PROF_WORK_UNIT_RUN);

    // Tell the supervisor the work unit is started.
    lm::message::Message handshakeMsg;
    lm::message::StartedWorkUnit* swuMsg = handshakeMsg.mutable_started_work_unit();
    swuMsg->set_work_unit_id(rwuMsg.work_unit_id());
    communicator->sendMessage(rwuMsg.supervisor_address(), &handshakeMsg);

    // Create the output message.
    bool hasOutput = false;
    lm::message::Message outputParent;
    lm::message::ProcessWorkUnitOutput* pwuMsg = outputParent.mutable_process_work_unit_output();
    pwuMsg->set_work_unit_id(rwuMsg.work_unit_id());
    google::protobuf::RepeatedPtrField<lm::message::WorkUnitOutput>* outputParts = pwuMsg->mutable_part_output();

    // Create the finished work unit message.
    lm::message::Message finalStateMsg;
    lm::message::FinishedWorkUnit* fwuMsg = finalStateMsg.mutable_finished_work_unit();
    fwuMsg->set_work_unit_id(rwuMsg.work_unit_id());
    fwuMsg->set_solver_type(rwuMsg.solver_type());

    // Choose the solver to use.
    lm::main::Solver* solver = NULL;
    if (rwuMsg.solver_type() == lm::types::SolverType::ME)
    {
        // Set the model for the solver.
        if (meSolver->needsReactionModel())
        {
            if (rwuMsg.has_reaction_model())
                meSolver->setReactionModel(rwuMsg.reaction_model());
            else
                throw Exception("Work Unit runner terminating, solver requires a reaction model but none was specified", properties.me_solver().c_str());
        }
        if (meSolver->needsDiffusionModel())
        {
            if (rwuMsg.has_diffusion_model())
                meSolver->setDiffusionModel(rwuMsg.diffusion_model());
            else
                throw Exception("Work Unit runner terminating, solver requires a diffusion model but none was specified", properties.me_solver().c_str());
        }

        // Set the order parameters for the solver
        if (rwuMsg.has_order_parameters())
        {
            meSolver->setOrderParameters(rwuMsg.order_parameters());
        }

        solver = meSolver;
    }
    else if (rwuMsg.solver_type() == lm::types::SolverType::DIFFUSION_PDE)
    {
        if (rwuMsg.has_microenv_model())
            pdeSolver->setMicroenvironmentModel(rwuMsg.microenv_model());
        else
            throw Exception("Work Unit runner terminating, solver requires a microenvironment model but none was specified", properties.diffusion_pde_solver().c_str());

        solver = pdeSolver;
    }
    else
    {
        throw RuntimeException("RunWorkUnit message did not contain a valid solver type.");
    }

    // Set the limits.
    if (rwuMsg.has_trajectory_limits())
        solver->setLimits(rwuMsg.trajectory_limits());

    // Set the barriers.
    if (rwuMsg.has_trajectory_barriers())
        solver->setBarriers(rwuMsg.trajectory_barriers());

    // Set the output options.
    if (rwuMsg.has_output_options())
        solver->setOutputOptions(rwuMsg.output_options());

    uint64_t totalSteps=0;
    hrtime totalTime=0;
    for (int i=0; i<rwuMsg.part_size(); i+=solver->getSimultaneousTrajectories())
    {
        // Reset the solver.
        solver->reset();

        // Configure the solver state for each simultaneous trajectory.
        for (int j=0; j<(int)solver->getSimultaneousTrajectories() && (i+j)<rwuMsg.part_size(); j++)
        {
            solver->setState(rwuMsg.part(i+j).initial_state(), j);
        }

        PROF_BEGIN(PROF_WORK_UNIT_RUN_PART);

        // Run the work unit.
        hrtime t1=getHrTime();
        totalSteps += solver->generateTrajectory(rwuMsg.max_steps());
        totalTime += getHrTime()-t1;

        PROF_END(PROF_WORK_UNIT_RUN_PART);

        PROF_BEGIN(PROF_WORK_UNIT_SAVE_PART);

        // Save the status and the state.
        for (int j=0; j<(int)solver->getSimultaneousTrajectories() && (i+j)<rwuMsg.part_size(); j++)
        {
            lm::message::WorkUnitStatus* status = fwuMsg->add_part_status();
            status->set_status(solver->getStatus(j));
            if (status->status() == lm::message::WorkUnitStatus::LIMIT_REACHED) status->mutable_limit_reached()->CopyFrom(solver->getLimitReached(j));
            lm::io::TrajectoryState* finalState = status->mutable_final_state();
            finalState->set_trajectory_id(rwuMsg.part(i+j).initial_state().trajectory_id());
            finalState->set_trajectory_started(true);
            solver->getState(finalState, j);
        }

        // Save the output.
        for (int j=0; j<(int)solver->getSimultaneousTrajectories() && (i+j)<rwuMsg.part_size(); j++)
        {
            lm::message::WorkUnitOutput* output = solver->getOutput(j);
            if (output->has_output())
            {
//                    output->set_condense_output(rwuMsg.output_options().condense_output());
                if (rwuMsg.output_options().has_record_name_prefix())
                {
                    output->set_record_name_prefix(rwuMsg.output_options().record_name_prefix());
                }

                outputParts->AddAllocated(output);
                hasOutput = true;
            }
            else
            {
                delete output;
            }
        }

        PROF_END(PROF_WORK_UNIT_SAVE_PART);
    }

//        // Cache some metadata for later use with condensed output, if needed
//        if (rwuMsg.output_options().condense_output())
//        {
//            condensedOutputAddress = rwuMsg.output_address();

//    #ifdef OPT_CPP11
//            // The unordered set that trajectoryList wraps (in c++11 mode) tends to reverse first and last trajectory, check
//            uint64_t _firstId = rwuMsg.part(0).initial_state().trajectory_id();
//            uint64_t _lastId = rwuMsg.part(rwuMsg.part_size() - 1).initial_state().trajectory_id();

//            firstTrajectoryIds.push_back(std::min(_firstId, _lastId));
//            lastTrajectoryIds.push_back(std::max(_firstId, _lastId));
//    #else
//            firstTrajectoryIds.push_back(rwuMsg.part(0).initial_state().trajectory_id());
//            lastTrajectoryIds.push_back(rwuMsg.part(rwuMsg.part_size() - 1).initial_state().trajectory_id());
//    #endif
//        }

    // Send the output.
    if (hasOutput) communicator->sendMessage(rwuMsg.output_address(), &outputParent);

    // Track the steps and time.
    fwuMsg->set_steps(totalSteps);
    fwuMsg->set_run_time(convertHrToSeconds(totalTime));

    // Tell the supervisor the work unit has finished.
    communicator->sendMessage(rwuMsg.supervisor_address(), &finalStateMsg);

    PROF_END(PROF_WORK_UNIT_RUN);
}

//bool WorkUnitRunner::flushOutputCondensed()
//{
//    bool hasOutput = false;

//    if (solver->workUnitCondenseOutput)
//    {
//        lm::message::Message outputParent;
//        lm::message::ProcessWorkUnitOutput* pwuMsg = outputParent.mutable_process_work_unit_output();
//        pwuMsg->set_work_unit_id(0);
//        google::protobuf::RepeatedPtrField<lm::message::WorkUnitOutput>* outputParts = pwuMsg->mutable_part_output();

//        lm::message::WorkUnitOutput* output = outputParts->Add();
//        if (solver->getOutputCondensed(output, firstTrajectoryIds, lastTrajectoryIds))
//        {
//            // reset the trajectory id vectors
//            firstTrajectoryIds.clear();
//            lastTrajectoryIds.clear();

//            hasOutput = true;

//            output->set_condense_output(solver->workUnitCondenseOutput);
//            output->set_record_name_prefix(solver->workUnitOutputPrefix);
//        }

//        // Send the output.
//        if (hasOutput) communicator->sendMessage(condensedOutputAddress, &outputParent);
//    }

//    return hasOutput;
//}

}
}
