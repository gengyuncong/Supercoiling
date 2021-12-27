/*
 * University of Illinois Open Source License
 * Copyright 2010 Luthey-Schulten Group,
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
 * Author(s): Elijah Roberts
 */

#include <iostream>
#include <string>
#include <map>
#include <vector>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include "lm/Exceptions.h"
#include "lm/Print.h"
#include "lm/String.h"
#include "lm/Version.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/input/ffpilot/FFPilotOptions.pb.h"
#include "lm/io/hdf5/Hdf5InputReader.h"
#include "lm/io/hdf5/SimulationFile.h"
#include "lm/types/OrderParameters.pb.h"
#include "lm/types/Tilings.pb.h"
#include "lptf/Profile.h"
#include "util.h"

using std::string;
using std::map;
using std::vector;
using lm::Print;
using lm::input::ffpilot::FFPilotOptions;
using lm::io::hdf5::Hdf5InputReader;
using lm::io::hdf5::Hdf5File;
using lm::types::OrderParameters;
using lm::types::Tilings;

void processParameter(string key, string range, string value, Hdf5File& file);
bool processOrderParameterOption(string key, string range, string value, OrderParameters* orderParameters);
bool processTilingOption(string key, string range, string value, Tilings* tilings);
bool processFFPilotOption(string key, string range, string value, FFPilotOptions* options, Tilings* tilings);
void printCopyright(int argc, char** argv);
void parseArguments(int argc, char** argv);
void printUsage(int argc, char** argv);

/**
 * The function being performed.
 */
static string function = "";

/**
 * The file to modify.
 */
static string filename = "";

/**
 * If we should print verbose messages.
 */
static bool verbose = false;

/**
 * The number of species.
 */
static uint numberSpecies=0;

/**
 * The parameters to set in the file.
 */
static map<string,string> parameters;

// Allocate the profile space.
PROF_ALLOC;

int main(int argc, char** argv)
{	
    PROF_INIT;
	try
	{
		printCopyright(argc, argv);
		parseArguments(argc, argv);
		
		if (function == "help")
		{
			printUsage(argc, argv);
		}
		else if (function == "version")
		{
		}
		else if (function == "set")
		{
		    // If the file doesn't exist, create it.
            bool newFile=false;
            struct stat fileStats;
		    if (stat(filename.c_str(), &fileStats) != 0)
		    {
		        if (numberSpecies > 0)
		            Hdf5File::create(filename, numberSpecies);
		        else
		            Hdf5File::create(filename);
                newFile = true;
            }

			// Open the file.
		    Hdf5File file(filename);

            // Read the reaction model.
            lm::input::ReactionModel model;
            if (!newFile)
            {
                file.getReactionModel(&model);
                numberSpecies = model.number_species();
            }

            // Read the order parameters list.
            bool setOrderParameterOptions = false;
            lm::types::OrderParameters orderParameters;
            if (!newFile) file.getOrderParameters(&orderParameters);

            // Read the tiling list.
            bool setTilingOptions = false;
            lm::types::Tilings tilings;
            if (!newFile) file.getTilings(&tilings);

            // Read the ffpilot options list.
            bool setFFPilotOptions = false;
            lm::input::ffpilot::FFPilotOptions ffpilotOptions;
            if (!newFile) file.getFFPilotOptions(&ffpilotOptions);

            // Set the parameters.
            Print::printf(Print::INFO, "Updating file %s.\n", filename.c_str());
		    for (std::map<std::string,string>::iterator it=parameters.begin(); it != parameters.end(); it++)
		    {
                // Get the key and value.
                string key=it->first;
                string value=it->second;

                // See if the user specified a range.
                string range("");
                size_t rangeStart=key.find_first_of("(");
                if (rangeStart != string::npos)
                {
                    range = key.substr(rangeStart);
                    key = key.substr(0, rangeStart);
                }

                Print::printf(Print::DEBUG, "Setting parameter: %s = %s for range %s", key.c_str(), value.c_str(), range.c_str());

                // Call the appropriate function based on the key.
                if (strncmp(key.c_str(), "OrderParameter", strlen("OrderParameter")) == 0)
                {
                    setOrderParameterOptions = processOrderParameterOption(key, range, value, &orderParameters);
                }
                else if (strncmp(key.c_str(), "Tiling", strlen("Tiling")) == 0)
                {
                    setTilingOptions = processTilingOption(key, range, value, &tilings);
                }
                else if (strncmp(key.c_str(), "FFPilot", strlen("ffpilot")) == 0)
                {
                    setFFPilotOptions = processFFPilotOption(key, range, value, &ffpilotOptions, &tilings);
                }
                else
                {
                    processParameter(key, range, value, file);
                }
		    }

            // Set the order parameters, if any were changed.
            if (setOrderParameterOptions)
                file.setOrderParameters(&orderParameters);

            // Set the tilings, if any were changed.
            if (setTilingOptions)
                file.setTilings(&tilings);

            // Set the ffpilot options, if any were changed.
            if (setFFPilotOptions)
                file.setFFPilotOptions(&ffpilotOptions);

            // Close the file.
            file.close();

            // Print the model.
            if (verbose)
            {
                Print::printf(Print::INFO, "Updated file %s:", filename.c_str());
                Hdf5InputReader reader;
                lm::input::Input input;
                reader.readFileInto(filename, &input);
                input.PrintDebugString();
            }
            Print::printf(Print::INFO, "Done.");
        }
		else
		{
			throw lm::CommandLineArgumentException("unknown function.");
		}
		return 0;
	}
    catch (lm::CommandLineArgumentException& e)
    {
    	std::cerr << "Invalid command line argument: " << e.what() << std::endl << std::endl;
        printUsage(argc, argv);
    }
    catch (lm::Exception& e)
    {
    	std::cerr << "Exception during execution: " << e.what() << std::endl;
    }
    catch (std::exception& e)
    {
    	std::cerr << "Exception during execution: " << e.what() << std::endl;
    }
    catch (...)
    {
    	std::cerr << "Unknown Exception during execution." << std::endl;
    }
    return -1;
}

void processParameter(string key, string range, string value, Hdf5File& file)
{
    file.setParameter(key, value);
}

bool processOrderParameterOption(string key, string range, string value, OrderParameters* orderParameters)
{
    if (numberSpecies == 0) THROW_EXCEPTION(lm::RuntimeException, "number of species must be specified when setting order parameter options");

    if (key == "OrderParameterCoefficients")
    {
        int32_t opIndex = 0;
        if (range != "")
        {
            vector<int> matrixDims;
            matrixDims.push_back(1);
            vector<int> opIndices = parseIndices(range, matrixDims);
            if (opIndices.size() != 1) THROW_EXCEPTION(lm::RuntimeException, "invalid order parameter indices value: %s", range.c_str());
            opIndex = opIndices[0];
        }

        // Create the indices for the species indices.
        vector<int> matrixDims;
        matrixDims.push_back(static_cast<int>(numberSpecies));
        vector<int> valueIndices = parseIndices("", matrixDims);

        // Parse the values.
        vector<double> values = parseValues(value);

        // Add any new order parameters to make sure we have the specified index.
        while (orderParameters->order_parameter().size() <= opIndex)
            orderParameters->add_order_parameter();

        // Get an reference to the order parameter.
        lm::types::OrderParameter* op = orderParameters->mutable_order_parameter(opIndex);
        op->Clear();

        // Set the type to be LinearCombinationOrderParameter::OPARAM_TYPE, which is 0.
        op->set_type(0);

        // Set the values.
        if (valueIndices.size() == values.size())
        {
            for (uint i=0; i<valueIndices.size(); i++)
            {
                op->add_species_index(uint32_t(valueIndices[i]));
                op->add_species_coefficient(values[i]);
            }
        }
        else
            THROW_EXCEPTION(lm::InvalidArgException, "For order parameter %d the number of order parameter value indices did not equal the number of values: %d, %d", opIndex, valueIndices.size(), values.size());
        return true;
    }

    Print::printf(Print::WARNING, "Unknown order parameter option: %s", key.c_str());
    return false;
}

bool processTilingOption(string key, string range, string value, Tilings* tilings)
{
    if (key == "TilingOrderParameters")
    {
        // Get the tiling index.
        int32_t tilingIndex = 0;
        if (range != "")
        {
            vector<int> matrixDims;
            matrixDims.push_back(1);
            vector<int> tilingIndices = parseIndices(range, matrixDims);
            if (tilingIndices.size() != 1) THROW_EXCEPTION(lm::RuntimeException, "invalid tiling indices value: %s", range.c_str());
            tilingIndex = tilingIndices[0];
        }

        // Get the values.
        vector<int> values = parseIntValues(value);

        // Add any new tilings to make sure we have the specified index.
        while (tilings->tiling().size() <= tilingIndex)
            tilings->add_tiling();

        // Get an reference to the tiling.
        lm::types::Tiling* tiling = tilings->mutable_tiling(tilingIndex);

        // Set the tiling's order parameters.
        ndarray<int32_t> op(utuple(static_cast<uint32_t>(values.size())),values.data());
        NDArraySerializer::serializeInto<int32_t>(tiling->mutable_order_parameter_indices(), op);

        return true;
    }
    else if (key == "TilingEdges")
    {
        // Get the tiling index.
        int32_t tilingIndex = 0;
        if (range != "")
        {
            vector<int> matrixDims;
            matrixDims.push_back(1);
            vector<int> tilingIndices = parseIndices(range, matrixDims);
            if (tilingIndices.size() != 1) THROW_EXCEPTION(lm::RuntimeException, "invalid tiling indices value: %s", range.c_str());
            tilingIndex = tilingIndices[0];
        }

        // Get the values.
        vector<double> values = parseValues(value);

        // Add any new tilings to make sure we have the specified index.
        while (tilings->tiling().size() <= tilingIndex)
            tilings->add_tiling();

        // Get an reference to the tiling.
        lm::types::Tiling* tiling = tilings->mutable_tiling(tilingIndex);

        // Set the tiling's order parameters.
        ndarray<double> edges(utuple(static_cast<uint32_t>(values.size())),values.data());
        NDArraySerializer::serializeInto<double>(tiling->mutable_edges(), edges);

        return true;
    }

    Print::printf(Print::WARNING, "Unknown tiling option: %s", key.c_str());
    return false;
}

bool processFFPilotOption(string key, string range, string value, FFPilotOptions* options, Tilings* tilings)
{
    if (numberSpecies == 0) THROW_EXCEPTION(lm::RuntimeException, "number of species must be specified when setting ffpilot options");

    if (key == "FFPilotBasins")
    {
        vector<int> matrixDims;

        // Parse the indices and values.
        matrixDims.push_back(2);
        matrixDims.push_back(static_cast<int>(numberSpecies));
        vector<utuple> indices = parseIndicesAsTuple(range, matrixDims);
        vector<int> values = parseIntValues(value);
        if (indices.size() != values.size()) THROW_EXCEPTION(lm::RuntimeException, "The number of indices must equal the number of values", indices.size(), values.size());

        // Set the values.
        ndarray<int32_t> counts(utuple(2, static_cast<uint32_t>(numberSpecies)));
        for (uint i=0; i<indices.size(); i++)
            counts[indices[i]] = values[i];

        // Set the counts in the options.
        NDArraySerializer::serializeInto<int32_t>(options->mutable_phase_zero_basins(), counts);

        return true;
    }
    else if (key == "FFPilotErrorGoal")
    {
        options->set_error_goal(atof(value.c_str()));
        return true;
    }
    else if (key == "FFPilotFallbackMethod")
    {
        if (value == "TILING")
            options->set_fallback_method(lm::input::ffpilot::FFPilotOptions::FALLBACK_TILING_EDGE);
        else if (value == "BASIN")
            options->set_fallback_method(lm::input::ffpilot::FFPilotOptions::FALLBACK_BASIN);
        else
            throw lm::CommandLineArgumentException("unknown value for FFPilotFallbackMethod", value.c_str());
        return true;
    }
    else if (key == "FFPilotPilotSkip")
    {
        options->set_pilot_skip(atoi(value.c_str()) != 0);
        return true;
    }
    else if (key == "FFPilotPilotStageCrossings")
    {
        options->set_pilot_stage_crossings(atoi(value.c_str()));
        return true;
    }
    else if (key == "FFPilotProdSkip")
    {
        options->set_prod_skip(atoi(value.c_str()) != 0);
        return true;
    }
    else if (key == "FFPilotProdTrajectoryCounts")
    {
        // Make sure we have a tiling.
        if (options->tiling_index() >= tilings->tiling_size()) THROW_EXCEPTION(lm::RuntimeException, "FFPilotProdTrajectoryCounts option requires a valid tiling index", options->tiling_index(), tilings->tiling_size());

        vector<int> matrixDims;

        // Parse the indices and values.
        int numberEdges = static_cast<int>(tilings->tiling(options->tiling_index()).edges().shape(0));
        matrixDims.push_back(2);
        matrixDims.push_back(numberEdges);
        vector<utuple> indices = parseIndicesAsTuple(range, matrixDims);
        vector<int> values = parseIntValues(value);
        ndarray<uint64_t> counts(utuple(2, static_cast<uint>(numberEdges)));
        if (values.size() == 1)
        {
            // Set all of the values to be the same.
            for (uint i=0; i<indices.size(); i++)
                counts[indices[i]] = static_cast<uint64_t>(values[0]);
        }
        else if(indices.size() == values.size())
        {
            // Set the values.
            for (uint i=0; i<indices.size(); i++)
                counts[indices[i]] = static_cast<uint64_t>(values[i]);
        }
        else
        {
            THROW_EXCEPTION(lm::RuntimeException, "The number of indices must be singular or equal the number of values", indices.size(), values.size());
        }

        // Set the counts in the options.
        NDArraySerializer::serializeInto<uint64_t>(options->mutable_prod_trajectory_counts(), counts);
        return true;
    }
    else if (key == "FFPilotSamplingMultipliers")
    {
        // Make sure we have a tiling.
        if (options->tiling_index() >= tilings->tiling_size()) THROW_EXCEPTION(lm::RuntimeException, "FFPilotSamplingMultiplier option requires a valid tiling index", options->tiling_index(), tilings->tiling_size());

        vector<int> matrixDims;

        // Parse the indices and values.
        int numberEdges = static_cast<int>(tilings->tiling(options->tiling_index()).edges().shape(0));
        matrixDims.push_back(2);
        matrixDims.push_back(numberEdges);
        vector<utuple> indices = parseIndicesAsTuple(range, matrixDims);
        vector<double> values = parseValues(value);
        ndarray<double> multipliers(utuple(2, static_cast<uint>(numberEdges)));
        if (values.size() == 1)
        {
            // Set all of the values to be the same.
            for (uint i=0; i<indices.size(); i++)
                multipliers[indices[i]] = static_cast<double>(values[0]);
        }
        else if(indices.size() == values.size())
        {
            // Set the values.
            for (uint i=0; i<indices.size(); i++)
                multipliers[indices[i]] = static_cast<double>(values[i]);
        }
        else
        {
            THROW_EXCEPTION(lm::RuntimeException, "The number of indices must be singular or equal the number of values: %d,%d", indices.size(), values.size());
        }

        // Set the counts in the options.
        NDArraySerializer::serializeInto<double>(options->mutable_optimize_sampling_multipliers(), multipliers);
        return true;
    }
    else if (key == "FFPilotTrajectoryOutput")
    {
        options->set_trajectory_output(atoi(value.c_str()) != 0);
        return true;
    }

    Print::printf(Print::WARNING, "Unknown ffpilot option: %s", key.c_str());
    return false;
}


/**
 * This function prints the copyright notice.
 */
void printCopyright(int argc, char** argv) {

    std::cout << argv[0] << " v" << VERSION_NUM << " build " << BUILD_TIMESTAMP << std::endl;
    std::cout << "Copyright (C) " << COPYRIGHT_DATE << " Luthey-Schulten Group," << std::endl;
	std::cout << "University of Illinois at Urbana-Champaign." << std::endl;
	std::cout << std::endl;
}

/**
 * Parses the command line arguments.
 */
void parseArguments(int argc, char** argv)
{
    //Parse any arguments.
    for (int i=1; i<argc; i++)
    {
        char *option = argv[i];
        while (*option == ' ') option++;
        
        //See if the user is trying to get help.
        if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
        	function = "help";
        	break;
        }
        
        //See if the user is trying to get the version info.
        else if (strcmp(option, "-v") == 0 || strcmp(option, "--version") == 0) {
        	function = "version";
        	break;
        }
            
        //See if the user is trying to get verbose information.
        else if (strcmp(option, "--verbose") == 0) {
            verbose = true;
        }

        //See if the user is trying to specify a filename.
        else if (i == 1)
        {
        	function = "set";
    		filename = option;
        }
        
        //See if the user is trying to specify a species count.
        else if (i == 2 && strstr(option, "=") == NULL)
        {
            numberSpecies = static_cast<uint32_t>(atoi(option));
        }

        else if (strstr(option, "=") != NULL)
        {
            char * separator=strstr(option, "=");
            if (separator != NULL && separator > option && strlen(separator) > 1)
            {
                string key(option, separator-option);
                string value(separator+1);
                parameters[key] = value;
            }
        }
             
        //This must be an invalid option.
        else {
            throw lm::CommandLineArgumentException(option);
        }
    }
}

/**
 * Prints the usage for the program.
 */
void printUsage(int argc, char** argv)
{
	std::cout << "Usage: " << argv[0] << " (-h|--help)" << std::endl;
	std::cout << "Usage: " << argv[0] << " (-v|--version)" << std::endl;
    std::cout << "Usage: " << argv[0] << " [--verbose] filename [number_species] (key=value)+" << std::endl;
	std::cout << std::endl;
}
