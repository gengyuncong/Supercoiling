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

#include <string>
#include <list>
#include <vector>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cstdio>
#include <ctime>

#include <sys/stat.h>
#if defined(MACOSX)
#include <sys/sysctl.h>
#elif defined(LINUX)
#include <sys/sysinfo.h>
#endif
#include "hrtime.h"
#include "lm/ClassFactory.h"
#ifdef OPT_CUDA
#include "lm/Cuda.h"
#endif
#include "lm/Exceptions.h"
#include "lm/Print.h"
#include "lm/Types.h"
#include "lm/Version.h"
#include "lm/main/Globals.h"
#include "lm/resource/ResourceMap.h"

using std::string;
using std::vector;

/**
 * The function being performed.
 */
string functionOption = "interpreter";

/**
 * The name of the file containing the simulation input.
 */
vector<string> simulationInputFilenames;

/**
 * The name of the file containing the simulation output.
 */
string simulationOutputFilename;

/**
 * The input reader to use for the simulations.
 */
string inputReaderClassName;

/**
 * The output writer to use for the simulations.
 */
string outputWriterClassName;

/**
 * The prefix to use in sfile record names.
 */
string recordNamePrefixGlobal;


/*
 * Options related to parameter sweep.
 */

/**
 * The parameter sweep filename.
 */
std::string parameterSweepFilename;

/**
 * The parameter sweep filename.
 */
int parameterSweepParallel;


/**
 * The number of replicates of the simulation that should be performed.
 */
vector<uint64_t> replicates;

/**
  * The number of replicates to send in each work unit.
  */
uint64_t replicateBatchSize;

/**
 * If messages should be printed about individual replicates.
 */
bool replicatePrintMessages;

/**
 * The interval at which the performance statistics should be outputted.
 */
time_t performanceOutputInterval;

/**
 * The interval at which the results file should be checkpointed.
 */
time_t checkpointInterval = 0;

/**
 * If a global abort signal has been received.
 */
volatile bool globalAbort = false;

/**
 * The communicator to use.
 */
string communicatorClassName;

/**
 * The hypervisor to use for the job.
 */
string hypervisorClassName;

/**
 * The supervisor to use for the simulations.
 */
string supervisorClassName;

/**
 * The solver to use for the simulations.
 */
string meSolverClassName;

/**
 * The diffusion PDE solver to use for the simulations.
 */
string diffusionPDESolverClassName;

/**
 * The filename for the resource list.
 */
string resourceFilename = "";

/**
 * The format for the resource file.
 */
lm::resource::ResourceMap::ResourceFileFormat resourceFileFormat;

/**
 * The number of cpu cores assigned to each process.
 */
int cpuCores;

/**
 * The number of cpu cores to assign per runner (can be a fraction, e.g., 1/2, 1/4, etc).
 */
double cpuCoresPerRunner;

/**
 * Whether we should use CPU affinity.
 */
bool useCPUAffinity;

/**
 * The number gpu devices assigned to each process.
 */
int gpuDevices;

/**
 * The number of gpu devices to assign per runner (can be a fraction, e.g., 1/2, 1/4, etc).
 */
double gpuDevicesPerRunner;

/**
 * Whether we should print the cuda device capabilities on startup.
 */
bool shouldPrintGPUCapabilities;

/**
 * Whether we should reserve a core for the output thread.
 */
bool shouldReserveOutputCore;

/**
 * Flag to indicate that we're running a test of the program's input and output
 */
bool ioTestFlag;

/*
 * Options related to ffpilot sampling.
 */

/**
 * Flag to whether we should print ffpilot stage messages.
 */
bool ffpilotPrintStageMessages;

