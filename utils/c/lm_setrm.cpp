/*
 * University of Illinois Open Source License
 * Copyright 2008-2011 Luthey-Schulten Group,
 * Copyright 2012-2019 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Luthey-Schulten Group
 *                  University of Illinois at Urbana-Champaign
 *                  http://www.scs.uiuc.edu/~schulten
 *
 * Developed by: Roberts Group
 *                  Johns Hopkins University
 *                  http://biophysics.jhu.edu/roberts/
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

#include <iostream>
#include <list>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>

#include "util.h"
#include "lm/Exceptions.h"
#include "lm/String.h"
#include "lm/Version.h"
#include "lm/cme/CMEOrderParameters.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/io/hdf5/SimulationFile.h"
#include "lptf/Profile.h"
#include "robertslab/pbuf/NDArraySerializer.h"

void printCopyright(int argc, char** argv);
void parseArguments(int argc, char** argv);
void printUsage(int argc, char** argv);

using std::list;
using std::pair;
using std::string;
using std::vector;
using lm::input::ReactionModel;
using lm::io::hdf5::Hdf5File;
using lm::parseUint;
using robertslab::pbuf::NDArraySerializer;

/**
 * The function being performed.
 */
string function = "";

/**
 * The file to modify.
 */
string filename = "";

/**
 * The parameters to set in the file.
 */
list<pair<string,string> > parameters;

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
                Hdf5File::create(filename);
                newFile = true;
		    }

			// Open the file.
		    Hdf5File file(filename);

		    // Read the reaction model.
            ReactionModel model;
		    if (!newFile) file.getReactionModel(&model);

            printf("Loaded model:\n");
            model.PrintDebugString();

            // Set the parameters.
		    for (list<pair<string,string> >::iterator it=parameters.begin(); it != parameters.end(); it++)
		    {
		        string key=it->first;
                string value=it->second;

                // See if the user specified a range.
                vector<int> matrixDims;
		        string range("");
                size_t rangeStart=key.find_first_of("(");
		        if (rangeStart != string::npos)
		        {
		            range = key.substr(rangeStart);
	                key = key.substr(0, rangeStart);
		        }

                // Set the parameter.
                if (key == "numberSpecies")
                {
                	vector<double> values = parseValues(value);
                    if (values.size() != 1) throw lm::Exception("A single value must be specified for the number of species.");
                    model.set_number_species((uint)lround(values[0]));
                }
                else if (key == "numberReactions")
                {
                	vector<double> values = parseValues(value);
                    if (values.size() != 1) throw lm::Exception("A single value must be specified for the number of reactions.");
                    model.set_number_reactions((uint)lround(values[0]));
                    for (uint i=model.reaction_size(); i<model.number_reactions(); i++)
                        model.add_reaction();
                }
                else if (key == "InitialSpeciesCounts")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.initial_species_count_size(); i<model.number_species(); i++)
                        model.add_initial_species_count(0);

                    // Parse the indices.
                	matrixDims.push_back(model.number_species());

    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (indices.size() >= 1 && values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_initial_species_count(indices[i],(uint)lround(values[0]));
                        }
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_initial_species_count(indices[i],(uint)lround(values[i]));
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());
                }
                else if (key == "ReactionTypes")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.reaction_size(); i<model.number_reactions(); i++)
                    	model.add_reaction();

                    // Parse the indices.
                	matrixDims.push_back(model.number_reactions());
    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (indices.size() >= 1 && values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.mutable_reaction(indices[i])->set_type((uint)lround(values[0]));
                        }
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.mutable_reaction(indices[i])->set_type((uint)lround(values[i]));
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());
                }
                else if (key == "ReactionRateConstants")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.reaction_size(); i<model.number_reactions(); i++)
                    	model.add_reaction();

                    // Parse the indices.
                	matrixDims.push_back(model.number_reactions());
                	matrixDims.push_back(10);
    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (indices.size() >= 1 && values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                        	uint rindex=indices[i]/10;
                        	uint kindex=indices[i]%10;
                        	for (uint j=model.reaction(rindex).rate_constant_size(); j<=kindex; j++)
                        		model.mutable_reaction(rindex)->add_rate_constant(0.0);
                            model.mutable_reaction(rindex)->set_rate_constant(kindex, values[0]);
                        }
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                        	uint rindex=indices[i]/10;
                        	uint kindex=indices[i]%10;
                        	for (uint j=model.reaction(rindex).rate_constant_size(); j<=kindex; j++)
                        		model.mutable_reaction(rindex)->add_rate_constant(0.0);
                        	model.mutable_reaction(rindex)->set_rate_constant(kindex, values[i]);
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());
                }
                else if (key == "StoichiometricMatrix")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.stoichiometric_matrix_size(); i<model.number_species()*model.number_reactions(); i++)
                        model.add_stoichiometric_matrix(0);

                    // Parse the indices.
                	matrixDims.push_back(model.number_species());
                	matrixDims.push_back(model.number_reactions());
    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_stoichiometric_matrix(indices[i], (int)lround(values[0]));
                        }
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_stoichiometric_matrix(indices[i], (int)lround(values[i]));
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());
                }
                else if (key == "DependencyMatrix")
                {
                    // Make sure the matrix is the correct size.
                    for (uint i=model.dependency_matrix_size(); i<model.number_species()*model.number_reactions(); i++)
                        model.add_dependency_matrix(0);

                    // Parse the indices.
                	matrixDims.push_back(model.number_species());
                	matrixDims.push_back(model.number_reactions());
    		        vector<int> indices = parseIndices(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Set the values.
                    if (values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_dependency_matrix(indices[i], (uint)lround(values[0]));
                        }
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                        {
                            model.set_dependency_matrix(indices[i], (uint)lround(values[i]));
                        }
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());
                }
                else if (key == "noiseProcesses")
                {
                    vector<double> values = parseValues(value);
                    if (values.size() != 1) throw lm::Exception("A single value must be specified for the number of noise processes.");
                    model.mutable_noise_model()->set_number_processes(static_cast<uint>(lround(values[0])));
                }
                else if (key == "noiseUpdateInterval")
                {
                    vector<double> values = parseValues(value);
                    if (values.size() != 1) throw lm::Exception("A single value must be specified for the process update interval.");
                    model.mutable_noise_model()->set_process_update_interval(values[0]);
                }
                else if (key == "NoiseTypes")
                {
                    // Parse the indices.
                    matrixDims.push_back(model.noise_model().number_processes());
                    vector<utuple> indices = parseIndicesAsTuple(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Create the array.
                    ndarray<uint32_t> A(utuple(model.noise_model().number_processes()));

                    // Keep the any previous values.
                    if (model.noise_model().process_types().shape().size() > 0)
                        A = NDArraySerializer::deserialize<uint32_t>(model.noise_model().process_types());

                    // Set the values.
                    if (indices.size() > 1 && values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                            A[indices[i]] = static_cast<uint>(lround(values[0]));
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                            A[indices[i]] = static_cast<uint>(lround(values[i]));
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());

                    // Store the matrix in the model.
                    NDArraySerializer::serializeInto(model.mutable_noise_model()->mutable_process_types(), A);

                }
                else if (key == "NoiseParameters")
                {
                    // Parse the indices.
                    matrixDims.push_back(model.noise_model().number_processes());
                    matrixDims.push_back(10);
                    vector<utuple> indices = parseIndicesAsTuple(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Create the array.
                    ndarray<double> A(utuple(model.noise_model().number_processes(),10));

                    // Keep the any previous values.
                    if (model.noise_model().process_parameters().shape().size() > 0)
                        A = NDArraySerializer::deserialize<double>(model.noise_model().process_parameters());

                    // Set the values.
                    if (indices.size() > 1 && values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                            A[indices[i]] = values[0];
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                            A[indices[i]] = values[i];
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());

                    // Store the matrix in the model.
                    NDArraySerializer::serializeInto(model.mutable_noise_model()->mutable_process_parameters(), A);

                }
                else if (key == "NoiseDependencies")
                {
                    // Parse the indices.
                    matrixDims.push_back(model.number_reactions());
                    matrixDims.push_back(10);
                    matrixDims.push_back(model.noise_model().number_processes());
                    vector<utuple> indices = parseIndicesAsTuple(range, matrixDims);
                    vector<double> values = parseValues(value);

                    // Create the array.
                    ndarray<uint32_t> A(utuple(model.number_reactions(),10,model.noise_model().number_processes()));

                    // Keep the any previous values.
                    if (model.noise_model().reaction_dependencies().shape().size() > 0)
                        A = NDArraySerializer::deserialize<uint32_t>(model.noise_model().reaction_dependencies());

                    // Set the values.
                    if (indices.size() > 1 && values.size() == 1)
                    {
                        for (uint i=0; i<indices.size(); i++)
                            A[indices[i]] = static_cast<uint>(lround(values[0]));
                    }
                    else if (indices.size() == values.size())
                    {
                        for (uint i=0; i<indices.size(); i++)
                            A[indices[i]] = static_cast<uint>(lround(values[i]));
                    }
                    else
                        throw lm::Exception("The number of indices must equal the number of values", indices.size(), values.size());

                    // Store the matrix in the model.
                    NDArraySerializer::serializeInto(model.mutable_noise_model()->mutable_reaction_dependencies(), A);

                }
                else
                {
                    throw lm::Exception("Unknown key: ",key.c_str());
                }
		    }

		    // Set the model.
		    file.setReactionModel(&model);

            // Close the file.
		    file.close();

            // Print the model.
            printf("\nUpdated model:\n");
            model.PrintDebugString();
            printf("Done.\n");
		}
		else
		{
			throw lm::CommandLineArgumentException("unknown function.");
		}
		return 0;
	}
    catch (lm::CommandLineArgumentException & e)
    {
    	std::cerr << "Invalid command line argument: " << e.what() << std::endl << std::endl;
        printUsage(argc, argv);
    }
    catch (lm::Exception & e)
    {
    	std::cerr << "Exception during execution: " << e.what() << std::endl;
    }
    catch (std::exception & e)
    {
    	std::cerr << "Exception during execution: " << e.what() << std::endl;
    }
    catch (...)
    {
    	std::cerr << "Unknown Exception during execution." << std::endl;
    }
    return -1;
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
            
        //See if the user is trying to specify a filename.
        else if (i == 1)
        {
        	function = "set";
    		filename = option;
        }
        
        //See if the user is trying to specify a key value pair.
        else if (strstr(option, "=") != NULL)
        {
            char * separator=strstr(option, "=");
            if (separator != NULL && separator > option && strlen(separator) > 1)
            {
                string key(option, separator-option);
                string value(separator+1);
                parameters.push_back(pair<string,string>(key,value));
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
	std::cout << "Usage: " << argv[0] << " filename (key=value)+" << std::endl;
	std::cout << std::endl;
}
