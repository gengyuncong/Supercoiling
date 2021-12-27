/*
 * University of Illinois Open Source License
 * Copyright 2008-2011 Luthey-Schulten Group,
 * Copyright 2012-2014 Roberts Group,
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
#include "lm/String.h"
#include "lm/Types.h"
#include "lm/Version.h"
#include "lm/main/Globals.h"
#include "lm/main/MainArgs.h"

using std::string;
using std::vector;


void printCopyright(int argc, char** argv)
{
    std::cout << "Lattice Microbe ES v" << VERSION_NUM;
    std::cout << " " << BUILD_TYPE << " build " << BUILD_TIMESTAMP;
    std::cout << " in " << (sizeof(uintv_t)*8) << "-bit mode";
    std::cout << " with options";
#ifdef OPT_AVX
    std::cout << " AVX";
#endif
#ifdef OPT_CPP11
    std::cout << " CPP11";
#endif
#ifdef OPT_CUDA
    std::cout << " CUDA";
#endif
#ifdef OPT_FMA
    std::cout << " FMA";
#endif
#ifdef OPT_MPI
    std::cout << " MPI";
#endif
#ifdef OPT_SBML
    std::cout << " SBML";
#endif
#ifdef OPT_SNAPPY
    std::cout << " SNAPPY";
#endif
#ifdef OPT_SVML
    std::cout << " SVML";
#endif
    std::cout << "." << std::endl;
    std::cout << "Copyright (C) " << COPYRIGHT_DATE_JHU << " Roberts Group, Johns Hopkins University." << std::endl;
    std::cout << "Copyright (C) " << COPYRIGHT_DATE << " Luthey-Schulten Group, University of Illinois at Urbana-Champaign." << std::endl << std::endl;
}

/**
 * Parses the command line arguments.
 */
void parseArguments(int argc, char** argv, bool printInfo)
{
    // Set any default options.
#ifdef OPT_MPI
    communicatorClassName = "lm::mpi::MPICommunicator";
#else
    communicatorClassName = "lm::message::LocalCommunicator";
#endif
    cpuCores = -1;
    cpuCoresPerRunner = 1.0;
    gpuDevices = -1;
#ifdef OPT_CUDA
    gpuDevicesPerRunner = 1.0;
#else
    gpuDevicesPerRunner = 0.0;
#endif
    ioTestFlag = false;
    replicateBatchSize = 1;

    performanceOutputInterval = 600;
    checkpointInterval = 0;

    // Parameter sweep options.
    parameterSweepFilename = "";
    parameterSweepParallel = 1;

    // Replicate sampling args.
    replicates.push_back(1);
    replicatePrintMessages = true;

    // FFPilot args.
    ffpilotPrintStageMessages = true;

    shouldPrintGPUCapabilities = true;
    shouldReserveOutputCore = true;
    simulationInputFilenames.clear();
    simulationOutputFilename = "";
    hypervisorClassName = "lm::simulation::SimpleHypervisor";
    supervisorClassName = "lm::replicates::ReplicateSupervisor";
    inputReaderClassName = "lm::input::Input";
    outputWriterClassName = "lm::io::hdf5::Hdf5OutputWriter";
    useCPUAffinity = false;

#ifdef OPT_AVX
    //mesolverClassName = "lm::avx::GillespieDSolverAVX";
    meSolverClassName = "lm::cme::GillespieDSolver";
#else
    meSolverClassName = "lm::cme::GillespieDSolver";
#endif

#ifdef OPT_AVX
    //diffusionPDESolverClassName = "lm::avx::ExplicitFiniteDifferenceSolverAVX";
    diffusionPDESolverClassName = "lm::pde::ExplicitFiniteDifferenceSolver";
#else
    diffusionPDESolverClassName = "lm::pde::ExplicitFiniteDifferenceSolver";
#endif

    /*
     * Simulation presets. Preparsed here, since they may affect the default values of other cmd line args
     */
    for (int i=1; i<argc; i++)
    {
        char *option = argv[i];
        while (*option == ' ') option++;

        if ((strcmp(option, "-sw") == 0 || strcmp(option, "--sweep-params") == 0))
        {
            hypervisorClassName = "lm::paramsweep::ParameterSweepHypervisor";
        }

        // See if the user is trying to perform a replicate sampling simulation.
        else if ((strcmp(option, "-rs") == 0 || strcmp(option, "--replicate-sampling") == 0))
        {
            supervisorClassName = "lm::replicates::ReplicateSupervisor";
        }

        // See if the user is trying to use forward flux pilot sampling.
        else if ((strcmp(option, "-ffp") == 0 || strcmp(option, "--ffpilot") == 0))
        {
            supervisorClassName = "lm::ffpilot::FFPilotSupervisor";
            replicates.clear();
        }

        // See if the user is trying to perform a microenvironment simulation.
        else if ((strcmp(option, "-me") == 0 || strcmp(option, "--microenvironment") == 0))
        {
            supervisorClassName = "lm::microenv::MicroenvironmentSupervisor";
            outputWriterClassName = "lm::io::sfile::SFileOutputWriter";
        }
    }


    // Parse any arguments.
    for (int i=1; i<argc; i++)
    {
        char *option = argv[i];
        while (*option == ' ') option++;
        
        //See if the user is trying to get help.
        if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
        	functionOption = "help";
        	break;
        }
        
        //See if the user is trying to get the version info.
        else if (strcmp(option, "-v") == 0 || strcmp(option, "--version") == 0) {
        	functionOption = "version";
            break;
        }

        //See if the user is trying to perfrom a debug test.
        else if (strcmp(option, "--debug") == 0) {
            functionOption = "debug";
            break;
        }

        //See if the user is trying to execute an iotest.
        else if (strcmp(option, "-iotest") == 0 || strcmp(option, "--input-ouput-test") == 0)
        {
            functionOption = "iotest";

            // Get the filename.
            parseStringListArg(simulationInputFilenames, argv[++i]);
        }

        //See if the user is trying to get the device info.
        else if (strcmp(option, "-l") == 0 || strcmp(option, "--list-devices") == 0) {
            functionOption = "devices";
        }

        // See if the user is trying to execute a simulation.
        else if ((strcmp(option, "-f") == 0 || strcmp(option, "--file") == 0) && i < (argc-1))
        {
            functionOption = "simulation";
            parseStringListArg(simulationInputFilenames, argv[++i], false);
        }
        else if (strncmp(option, "--file=", strlen("--file=")) == 0)
        {
            functionOption = "simulation";
            parseStringListArg(simulationInputFilenames, option+strlen("--file="), false);
        }

        //See if the user is trying to set the output format.
        else if ((strcmp(option, "-ff") == 0 || strcmp(option, "--output-format") == 0) && i < (argc-1))
        {
            outputWriterClassName=parseOutputFormatArg(argv[++i]);
        }
        else if (strncmp(option, "--output-format=", strlen("--output-format=")) == 0)
        {
            outputWriterClassName=parseOutputFormatArg(option+strlen("--output-format="));
        }

        //See if the user is trying to set the output filename.
        else if ((strcmp(option, "-fo") == 0 || strcmp(option, "--output-file") == 0) && i < (argc-1))
        {
            simulationOutputFilename=argv[++i];
        }
        else if (strncmp(option, "--output-file=", strlen("--output-file=")) == 0)
        {
            simulationOutputFilename=option+strlen("--output-file=");
        }

        //See if the user is trying to set the output record prefix.
        else if ((strcmp(option, "-fp") == 0 || strcmp(option, "--output-prefix") == 0) && i < (argc-1))
        {
            recordNamePrefixGlobal=argv[++i];
        }
        else if (strncmp(option, "--output-prefix=", strlen("--output-prefix=")) == 0)
        {
            recordNamePrefixGlobal=option+strlen("--output-prefix=");
        }

        /*
         * Options related to parameter sweep.
         */

        //See if the user is trying to turn on parameter sweep.
        else if ((strcmp(option, "-sw") == 0 || strcmp(option, "--sweep-params") == 0))
        {
            // Do nothing since the preset is taken care of above.
        }

        //See if the user is trying to set the parameter sweep filename.
        else if ((strcmp(option, "-swf") == 0 || strcmp(option, "--sweep-file") == 0) && i < (argc-1))
        {
            parameterSweepFilename=argv[++i];
        }
        else if (strncmp(option, "--sweep-file=", strlen("--sweep-file=")) == 0)
        {
            parameterSweepFilename=option+strlen("--sweep-file=");
        }

        //See if the user is trying to set the number of parallel parameter sweeps.
        else if ((strcmp(option, "-swp") == 0 || strcmp(option, "--sweep-parallel") == 0) && i < (argc-1))
        {
            parameterSweepParallel=atoi(argv[++i]);
        }
        else if (strncmp(option, "--sweep-parallel=", strlen("--sweep-parallel=")) == 0)
        {
            parameterSweepParallel=atoi(option+strlen("--sweep-parallel="));
        }


        /*
         * Options related to replicate sampling.
         */

        //See if the user is trying to set the replicates.
        else if ((strcmp(option, "-r") == 0 || strcmp(option, "--replicates") == 0) && i < (argc-1))
        {
            parseIntListArg(replicates, argv[++i]);
        }
        else if (strncmp(option, "--replicates=", strlen("--replicates=")) == 0)
        {
            parseIntListArg(replicates, option+strlen("--replicates="));
        }

        //See if the user is trying to set the parts per work unit (ie the trajectory batch size).
        else if ((strcmp(option, "-rb") == 0 || strcmp(option, "--replicate-batch-size") == 0) && i < (argc-1))
        {
            replicateBatchSize = uint64_t(atoll(argv[++i]));
        }
        else if (strncmp(option, "--replicate-batch-size=", strlen("--replicate-batch-size=")) == 0)
        {
            replicateBatchSize = uint64_t(atoll(option+strlen("--replicate-batch-size=")));
        }

        //See if the user is trying to turn off printing of completed replicates.
         else if ((strcmp(option, "-rnp") == 0 || strcmp(option, "--replicate-no-print") == 0))
         {
             replicatePrintMessages = false;
         }

        //See if the user is trying to set the performance output interval.
        else if ((strcmp(option, "-pf") == 0 || strcmp(option, "--perf-stats") == 0) && i < (argc-1))
        {
            performanceOutputInterval=parseTimeArg(argv[++i]);
        }
        else if (strncmp(option, "--perf-stats=", strlen("--perf-stats=")) == 0)
        {
            performanceOutputInterval=parseTimeArg(option+strlen("--perf-stats="));
        }

        //See if the user is trying to set the checkpoint interval.
        else if ((strcmp(option, "-ck") == 0 || strcmp(option, "--checkpoint") == 0) && i < (argc-1))
        {
            checkpointInterval=parseTimeArg(argv[++i]);
        }
        else if (strncmp(option, "--checkpoint=", strlen("--checkpoint=")) == 0)
        {
            checkpointInterval=parseTimeArg(option+strlen("--checkpoint="));
        }

        //See if the user is trying to set the solver.
        else if ((strcmp(option, "-sp") == 0 || strcmp(option, "--spatially-resolved") == 0))
        {
            meSolverClassName = "lm::rdme::NextSubvolumeSolver";
        }
        else if ((strcmp(option, "-ws") == 0 || strcmp(option, "--well-stirred") == 0))
        {
            meSolverClassName = "lm::cme::GillespieDSolver";
        }
        else if ((strcmp(option, "-sl") == 0 || strcmp(option, "--solver") == 0) && i < (argc-1))
        {
            meSolverClassName = argv[++i];
        }
        else if (strncmp(option, "--solver=", strlen("--solver=")) == 0)
        {
            meSolverClassName = option+strlen("--solver=");
        }
        else if ((strcmp(option, "-slp") == 0 || strcmp(option, "--pde-solver") == 0) && i < (argc-1))
        {
            diffusionPDESolverClassName = argv[++i];
        }
        else if (strncmp(option, "--pde-solver=", strlen("--pde-solver=")) == 0)
        {
            diffusionPDESolverClassName = option+strlen("--pde-solver=");
        }

        //See if the user is trying to set the node list.
        else if ((strcmp(option, "-n") == 0 || strcmp(option, "--nodelist") == 0) && i < (argc-1))
        {
            resourceFilename=argv[++i];
            resourceFileFormat = lm::resource::ResourceMap::NODELIST;
        }
        else if (strncmp(option, "--nodelist=", strlen("--nodelist=")) == 0)
        {
            resourceFilename=option+strlen("--nodelist=");
            resourceFileFormat = lm::resource::ResourceMap::NODELIST;
        }

        //See if the user is trying to set the resource map.
        else if ((strcmp(option, "-m") == 0 || strcmp(option, "--resource-map") == 0) && i < (argc-1))
        {
            resourceFilename=argv[++i];
            resourceFileFormat = lm::resource::ResourceMap::RESOURCE_MAP;
        }
        else if (strncmp(option, "--resource-map=", strlen("--resource-map=")) == 0)
        {
            resourceFilename=option+strlen("--resource-map=");
            resourceFileFormat = lm::resource::ResourceMap::RESOURCE_MAP;
        }

        //See if the user is trying to set the number of cpus.
        else if ((strcmp(option, "-c") == 0 || strcmp(option, "--cpu") == 0) && i < (argc-1))
        {
            cpuCores=atoi(argv[++i]);
        }
        else if (strncmp(option, "--cpu=", strlen("--cpu=")) == 0)
        {
            cpuCores=atoi(option+strlen("--cpu="));
        }

        //See if the user is trying to set the number of gpu devices per runner.
        else if ((strcmp(option, "-cr") == 0 || strcmp(option, "--cpus-per-runner") == 0 || strcmp(option, "--cpus-per-replicate") == 0) && i < (argc-1))
        {
             cpuCoresPerRunner=parseIntReciprocalArg(argv[++i]);
        }
        else if (strncmp(option, "--cpus-per-runner=", strlen("--cpus-per-runner=")) == 0)
        {
            cpuCoresPerRunner=parseIntReciprocalArg(option+strlen("--cpus-per-runner="));
        }
        else if (strncmp(option, "--cpus-per-replicate=", strlen("--cpus-per-replicate=")) == 0)
        {
             cpuCoresPerRunner=parseIntReciprocalArg(option+strlen("--cpus-per-replicate="));
        }


        //See if the user is trying to turn on cpu affinity.
        else if ((strcmp(option, "-ca") == 0 || strcmp(option, "--cpu-affinity") == 0))
        {
             useCPUAffinity = true;
        }

        //See if the user is trying to set the gpu devices.
        else if ((strcmp(option, "-g") == 0 || strcmp(option, "--gpu") == 0) && i < (argc-1))
        {
             gpuDevices=atoi(argv[++i]);
        }
        else if (strncmp(option, "--gpu=", strlen("--gpu=")) == 0)
        {
             gpuDevices=atoi(option+strlen("--gpu="));
        }

        //See if the user is trying to set the number of gpu devices per runner.
        else if ((strcmp(option, "-gr") == 0 || strcmp(option, "--gpus-per-runner") == 0 || strcmp(option, "--gpus-per-replicate") == 0) && i < (argc-1))
        {
             gpuDevicesPerRunner=parseIntReciprocalArg(argv[++i]);
        }
        else if (strncmp(option, "--gpus-per-runner=", strlen("--gpus-per-runner=")) == 0)
        {
             gpuDevicesPerRunner=parseIntReciprocalArg(option+strlen("--gpus-per-runner="));
        }
        else if (strncmp(option, "--gpus-per-replicate=", strlen("--gpus-per-replicate=")) == 0)
        {
            gpuDevicesPerRunner=parseIntReciprocalArg(option+strlen("--gpus-per-replicate="));
        }

        //See if the user is trying to turn off cuda capability printing.
        else if ((strcmp(option, "-nc") == 0 || strcmp(option, "--no-capabilities") == 0))
        {
             shouldPrintGPUCapabilities = false;
        }

        //See if the user is trying to turn off cuda capability printing.
        else if ((strcmp(option, "-nr") == 0 || strcmp(option, "--no-reserve-core") == 0))
        {
        	 shouldReserveOutputCore = false;
        }

        //See if the user is trying to set the supervisor directly.
        else if ((strcmp(option, "-su") == 0 || strcmp(option, "--supervisor") == 0) && i < (argc-1))
        {
            supervisorClassName = argv[++i];
        }
        else if (strncmp(option, "--supervisor=", strlen("--supervisor=")) == 0)
        {
            supervisorClassName = option+strlen("--supervisor=");
        }
        //See if the user is trying to do an input output test.
        else if ((strcmp(option, "-ioflag") == 0 || strcmp(option, "--do-io-test") == 0))
        {
            ioTestFlag = true;
        }
        //See if the user is trying to set the gpu devices.
        else if ((strcmp(option, "-so") == 0 || strcmp(option, "--shared-libraries") == 0) && i < (argc-1))
        {
            vector<string> sharedLibraries;
            parseStringListArg(sharedLibraries, argv[++i]);
            for (vector<string>::iterator it=sharedLibraries.begin(); it!=sharedLibraries.end(); it++)
                lm::ClassFactory::getInstance().registerClassesFromExternalLibrary(*it);
        }
        else if (strncmp(option, "--shared-libraries=", strlen("--shared-libraries=")) == 0)
        {
            vector<string> sharedLibraries;
            parseStringListArg(sharedLibraries, option+strlen("--shared-libraries="));
            for (vector<string>::iterator it=sharedLibraries.begin(); it!=sharedLibraries.end(); it++)
                lm::ClassFactory::getInstance().registerClassesFromExternalLibrary(*it);
        }

        // See if the user is trying to set the communicator.
        else if (strcmp(option, "--local") == 0)
        {
            communicatorClassName = "lm::message::LocalCommunicator";
        }
        else if (strcmp(option, "--mpi") == 0)
        {
            communicatorClassName = "lm::mpi::MPICommunicator";
        }
        else if (strcmp(option, "--mpi-async") == 0)
        {
            communicatorClassName = "lm::mpi::AsyncMPICommunicator";
        }

        /*
         * Options related to ffpilot.
         */

        //See if the user is trying to turn on ffpilot.
        else if ((strcmp(option, "-ffp") == 0 || strcmp(option, "--ffpilot") == 0))
        {
            // Do nothing since the preset is taken care of above.
        }

        //See if the user is trying to turn off printing of ffpilot messages.
        else if ((strcmp(option, "-ffpnp") == 0 || strcmp(option, "--ffpilot-no-print") == 0))
        {
            ffpilotPrintStageMessages = false;
        }


        /*
         * Simulation type arguments. Parsed earlier, but included here so the args are treated as valid
         */
        else if ((strcmp(option, "-rs") == 0 || strcmp(option, "--replicate-sampling") == 0)) continue;
        else if ((strcmp(option, "-me") == 0 || strcmp(option, "--microenvironment") == 0)) continue;

        // This must be an invalid option.
        else {
            throw lm::CommandLineArgumentException(option);
        }
    }

    // Figure out where to save the simulation output
    if (functionOption=="simulation")
    {
        if (simulationInputFilenames.size()==0)
            throw lm::CommandLineArgumentException("missing simulation input file.");

        if (outputWriterClassName.find("hdf5") != std::string::npos)
        {
            if (simulationOutputFilename=="")
            {
                // If the output name is blank, assume the first input name is the .lm file and set the output to be that .lm file
                simulationOutputFilename = simulationInputFilenames[0];
            }
            else if (simulationOutputFilename != simulationInputFilenames[0])
            {
                throw lm::CommandLineArgumentException("cannot specify separate input and output files with the hdf5 format.");
            }
        }
        else if (outputWriterClassName.find("sfile") != std::string::npos && simulationOutputFilename=="")
        {
            // default SFile output path is the input path with "_-_out.sfile" suffix
            simulationOutputFilename = lm::pathWithSuffix(simulationInputFilenames[0], ".sfile");
        }
        if (printInfo) lm::Print::printf(lm::Print::INFO, "Saving simulation output to %s using %s.", simulationOutputFilename.c_str(), outputWriterClassName.c_str());
    }

    // zero out the CUDA args if CUDA is off. Warn the user if we have to change any arg vals
    #ifndef OPT_CUDA
    // if cuda is off, printInfo the user if they try to set gpuDevices, but then set it to 0 anyway
    if (gpuDevices > 0)
    {
        if (printInfo) lm::Print::printf(lm::Print::WARNING, "attempting to set gpuDevices=%.2f, but CUDA support is turned off", gpuDevices);
        gpuDevices=0;
    }
    // if cuda is off, printInfo the user if they try to set gpuDevicesPerRunner, but set it to 0 anyway
    if (gpuDevicesPerRunner > 0)
    {
        if (printInfo) lm::Print::printf(lm::Print::WARNING, "attempting to set gpuDevicesPerRunner=%.2f, but CUDA support is turned off", gpuDevicesPerRunner);
        gpuDevicesPerRunner=0.0;
    }
    #endif /* OPT_CUDA */
}

string parseOutputFormatArg(char* option)
{
    if (strcmp(option, "hdf5") == 0)
        return "lm::io::hdf5::Hdf5OutputWriter";
    else if (strcmp(option, "sfile") == 0)
        return "lm::io::sfile::SFileOutputWriter";
    else if (strcmp(option, "log") == 0)
        return "lm::io::ConsoleOutputWriter";
    else if (strcmp(option, "null") == 0)
        return "lm::io::NullOutputWriter";
    throw lm::CommandLineArgumentException("Invalid output file format: %s", option);
}

void parseIntListArg(vector<uint64_t> & list, char* arg)
{
    list.clear();
    char * argbuf = new char[strlen(arg)+1];
    strcpy(argbuf,arg);
    char * pch = strtok(argbuf," ,;:\"");
    while (pch != NULL)
    {
        char * rangeDelimiter;
        if ((rangeDelimiter=strstr(pch,"-")) != NULL)
        {
            *rangeDelimiter='\0';
            int begin=atoi(pch);
            int end=atoi(rangeDelimiter+1);
            for (int i=begin; i<=end; i++)
                list.push_back(i);
        }
        else
        {
            if (strlen(pch) > 0) list.push_back(atoi(pch));
        }
        pch = strtok(NULL," ,;:");
    }
    delete[] argbuf;
}

void parseStringListArg(vector<string>& list, char* arg, bool clear)
{
    // Clear the lsit first, if necessary.
    if (clear) list.clear();

    char * argbuf = new char[strlen(arg)+1];
    strcpy(argbuf,arg);
    char * pch = strtok(argbuf," ,;:\"");
    while (pch != NULL)
    {
        if (strlen(pch) > 0) list.push_back(string(pch));
        pch = strtok(NULL," ,;:");
    }
    delete[] argbuf;
}

time_t parseTimeArg(char * arg)
{
    char * argbuf = new char[strlen(arg)+1];
    strcpy(argbuf,arg);
    char * pch = strtok(argbuf,":");

    // Parse the arguments into tokens.
    int tokenNumber=0;
    int tokens[3];
    while (pch != NULL)
    {
        if (tokenNumber < 3)
        {
            tokens[tokenNumber++] = atoi(pch);
        }
        else
        {
            delete[] argbuf;
            throw lm::CommandLineArgumentException(arg);
        }
        pch = strtok(NULL,":");
    }
    delete[] argbuf;

    // Calculate the time from the tokens.
    time_t time=0;
    if (tokenNumber == 1)
        time = tokens[0];
    else if (tokenNumber == 2)
        time = (60*tokens[0])+tokens[1];
    else if (tokenNumber == 3)
        time = (3600*tokens[0])+(60*tokens[1])+tokens[2];

    return time;
}

double parseIntReciprocalArg(char * arg)
{
    if (strlen(arg) >= 3 && arg[0] == '1' && arg[1] == '/')
    {
        return 1.0/(double)atoi(arg+2);
    }
    else
    {
        return (double)atoi(arg);
    }
}

/**
 * Prints the usage for the program.
 */
void printUsage(int argc, char** argv)
{
    std::cout << "Usage: lm (-h|--help)" << std::endl;
    std::cout << "Usage: lm (-v|--version)" << std::endl;
    std::cout << "Usage: lm (-l|--list-devices)" << std::endl;
    std::cout << "Usage: lm [OPTIONS] [SIM_OPTIONS] (-f input_filename_list | --file=input_filename_list)" << std::endl;
    std::cout << std::endl;
    std::cout << "OPTIONS" << std::endl;
    std::cout << "  -ff format        --output-format=format        The file format for the simulation output. Valid values are \"hdf5\" (default)|\"sfile\"|\"log\"|\"null\"." << std::endl;
    std::cout << "  -fo output_file   --output-file=output_filename The file for the simulation output, if different than the input file. Required for sfile, invalid for hdf5." << std::endl;
    std::cout << "  -fp record_prefix --output-prefix=record_prefix The prefix to use for the record names. Optional for sfile or hdf5 output." << std::endl;
    std::cout << "  -n node_file      --nodelist=node_file          A file containing the list of nodes on which to run, one line per available CPU core." << std::endl;
    std::cout << "  -m map_file       --resource-map=map_file       A file containing the map of resources to use: hostname processor_id_list gpu_id_list." << std::endl;
    std::cout << "  -c num_cpus       --cpu=num_cpus                The number of CPUs on which to execute (default all)." << std::endl;
    std::cout << "  -cr num           --cpus-per-runner=num         The number of CPUs (possibly fractional) to assign per runner, e.g. \"2\", \"1/4\" (default 1)." << std::endl;
    std::cout << "  -ca               --cpu-affinity                Turn on CPU affinity." << std::endl;
    std::cout << "  -g num_gpus       --gpu=num_gpus                The number of GPUs on which to execute (default all)." << std::endl;
    std::cout << "  -gr num           --gpus-per-runner=num         The number of GPUs (possibly fractional) to assign per runner, e.g. \"2\", \"1/4\" (default 1)." << std::endl;
    std::cout << "  -nc               --no-capabilities             Don't print the capabilities of the GPU devices." << std::endl;
    std::cout << "  -nr               --no-reserve-core             Don't reserve a CPU core for the output thread." << std::endl;
    std::cout << "  -so               --shared-libraries=libs       A comma delimited list of shared library to load." << std::endl;
    std::cout << "                    --local                       Use a local communicator on only this process." << std::endl;
    std::cout << "                    --mpi                         Use an MPI communicator across multiple processes." << std::endl;
    std::cout << "                    --mpi-async                   Use an anstnchronous MPI communicator across multiple processes." << std::endl;
    std::cout << std::endl;
    std::cout << "SIMULATION_OPTIONS" << std::endl;
    std::cout << "  -sp               --spatially-resolved          The simulations should use the spatially resolved reaction model (default)." << std::endl;
    std::cout << "  -ws               --well-stirred                The simulations should use the well-stirred reaction model." << std::endl;
    std::cout << "  -su supervisor    --supervisor=classname        Perform a simulation using the specified supervisor." << std::endl;
    std::cout << "  -sl solver        --solver=solver               The specific master equation solver class to use for the simulations." << std::endl;
    std::cout << "  -slp solver       --pde-solver=solver           The specific diffusion PDE solver class to use for the simulations." << std::endl;
    std::cout << "  -pf               --perf-stats=interval         Set the interval for writing performance statistics as hh:mm:ss (default 00:10:00)." << std::endl;
    std::cout << "  -ck               --checkpoint=interval         Enable checkpointing with the given interval as hh:mm:ss (default 00:00:00 -- disabled)." << std::endl;
    std::cout << std::endl;
    std::cout << "PARAMETER_SWEEP_OPTIONS" << std::endl;
    std::cout << "  -sw               --sweep-params                Perform a parameter sweep." << std::endl;
    std::cout << "  -swf              --sweep-file=filename         The csv file to read the parameter values from." << std::endl;
    std::cout << "  -swp              --sweep-parallel=num          The number of parameters to run in parallel (default 1)." << std::endl;
    std::cout << std::endl;
    std::cout << "REPLICATE_OPTIONS" << std::endl;
    std::cout << "  -rs               --replicate-sampling          Perform a replicate sampling simulation (default)." << std::endl;
    std::cout << "  -r replicates     --replicates=replicates       A list of replicates to run, e.g. \"0-9\", \"0,11,21\" (default 0)." << std::endl;
    std::cout << "  -rb size          --replicate-batch-size=size   The minimum number of replicates to include in each work unit (default 1)." << std::endl;
    std::cout << "  -rnp              --replicate-no-print          Turn off printing of per replicate messages." << std::endl;
    std::cout << std::endl;
    std::cout << "MICROENVIRONMENT_OPTIONS" << std::endl;
    std::cout << "  -me               --microenvironment            Perform a microenvironment simulation." << std::endl;
    std::cout << "  -r replicates     --replicates=replicates       A list of replicates to run, e.g. \"0-9\", \"0,11,21\" (default 0)." << std::endl;
    std::cout << std::endl;
    std::cout << "FORWARD_FLUX_OPTIONS" << std::endl;
    std::cout << "  -ffp              --ffpilot                    Perform a forward-flux pilot simulation." << std::endl;
    std::cout << "  -ffpnp            --ffpilot-no-print           Turn off printing of ffpilot stage information." << std::endl;
}











