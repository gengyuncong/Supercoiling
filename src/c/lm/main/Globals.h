/*
 * University of Illinois Open Source License
 * Copyright 2012-2016 Roberts Group,
 * All rights reserved.
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
 * - Neither the names of the Roberts Group, Johns Hopkins University,
 * nor the names of its contributors may be used to endorse or
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
#ifndef LM_MAIN_GLOBALS_H
#define LM_MAIN_GLOBALS_H

#include <ctime>
#include <list>
#include <string>
#include <vector>

#include "hrtime.h"
#include "lm/resource/ResourceMap.h"
#include "lm/Types.h"

/**
 * The function being performed.
 */
extern std::string functionOption;

/**
 * The names of the files containing the simulation input.
 */
extern std::vector<std::string> simulationInputFilenames;

/**
 * The name of the file containing the simulation output.
 */
extern std::string simulationOutputFilename;

/**
 * The input reader to use for the simulations.
 */
extern std::string inputReaderClassName;

/**
 * The output writer to use for the simulations.
 */
extern std::string outputWriterClassName;

/**
 * The prefix to use in sfile record names.
 */
extern std::string recordNamePrefixGlobal;



/*
 * Options related to parameter sweep.
 */

/**
 * The parameter sweep filename.
 */
extern std::string parameterSweepFilename;

/**
 * The parameter sweep filename.
 */
extern int parameterSweepParallel;



/*
 * Options related to replicate sampling.
 */

/**
 * The number of replicates of the simulation that should be performed.
 */
extern std::vector<uint64_t> replicates;

/**
  * The number of parts to send in each work unit.
  */
extern uint64_t replicateBatchSize;

/**
 * If messages should be printed about individual replicates.
 */
extern bool replicatePrintMessages;




/**
 * The interval at which the performance statistics should be outputted.
 */
extern time_t performanceOutputInterval;

/**
 * The interval at which the results file should be checkpointed.
 */
extern time_t checkpointInterval;

/**
 * If a global abort signal has been received.
 */
extern volatile bool globalAbort;

/**
 * The communicator to use.
 */
extern std::string communicatorClassName;

/**
 * The hypervisor to use for the job.
 */
extern std::string hypervisorClassName;

/**
 * The supervisor to use for the simulations.
 */
extern std::string supervisorClassName;

/**
 * The master equation solver to use for the simulations.
 */
extern std::string meSolverClassName;

/**
 * The diffusion PDE solver to use for the simulations.
 */
extern std::string diffusionPDESolverClassName;

/**
 * The filename for the resource list.
 */
extern std::string resourceFilename;

/**
 * The format for the resource file.
 */
extern lm::resource::ResourceMap::ResourceFileFormat resourceFileFormat;

/**
 * The number of cpu cores assigned to each process.
 */
extern int cpuCores;

/**
 * The number of cpu cores to assign per runner (can be a fraction, e.g., 1/2, 1/4, etc).
 */
extern double cpuCoresPerRunner;

/**
 * Whether we should use CPU affinity.
 */
extern bool useCPUAffinity;

/**
 * The number gpu devices assigned to each process.
 */
extern int gpuDevices;

/**
 * The number of gpu devices to assign per runner (can be a fraction, e.g., 1/2, 1/4, etc).
 */
extern double gpuDevicesPerRunner;

/**
 * Whether we should print the cuda device capabilities on startup.
 */
extern bool shouldPrintGPUCapabilities;

/**
 * Whether we should reserve a core for the output thread.
 */
extern bool shouldReserveOutputCore;

/**
 * Flag to indicate that we're running a test of the program's input and output
 */
extern bool ioTestFlag;

/*
 * Options related to ffpilot sampling.
 */

/**
 * Flag to whether we should print ffpilot stage messages.
 */
extern bool ffpilotPrintStageMessages;


#endif /* LM_MAIN_GLOBALS_H */
