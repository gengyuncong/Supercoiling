/*
 * University of Illinois Open Source License
 * Copyright 2008-2012 Luthey-Schulten Group,
 * Copyright 2012-2016 Roberts Group,
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

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <sys/stat.h>
#include <sbml/conversion/ConversionProperties.h>
#include <sbml/math/FormulaFormatter.h>
#include <sbml/SBMLDocument.h>
#include <sbml/SBMLReader.h>
#include <sbml/SBMLTypes.h>
#include <sbml/xml/XMLErrorLog.h>
#include "lm/ClassFactory.h"
#include "lm/Exceptions.h"
#include "lm/Print.h"
#include "lm/Version.h"
#include "lm/cme/CMEPropensityFunctions.h"
#include "lm/me/PropensityFunction.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/io/hdf5/SimulationFile.h"
#include "lptf/Profile.h"
#include "robertslab/sbml/SBMLImporterL3V1.h"
#include "robertslab/sbml/SBMLImporterL3V1COPASI.h"

void printCopyright(int argc, char** argv);
void parseArguments(int argc, char** argv);
void printUsage(int argc, char** argv);

using std::map;
using std::vector;
using std::string;
using lm::Exception;
using lm::input::ReactionModel;
using lm::io::hdf5::Hdf5File;
using lm::Print;

/**
 * The function being performed.
 */
string function = "";

/**
 * The lm file to output the model into.
 */
string outputFilename = "";

/**
 * The sbml file to import the model from.
 */
string inputFilename = "";

/**
 * Whether to ignore errors in the SBML file.
 */
bool ignoreErrors = false;

/**
 * Whether to ignore unmatched reactions.
 */
bool ignoreUnmatchedReactions = false;

/**
 * Whether to print out more verbose messages.
 */
bool verbose = false;
bool reallyVerbose = false;

/**
 * If the file was generated by COPASI, we willprocess it differently.
 */
bool isCopasi = false;

/**
 * If rate constants are given in terms of concentrations.
 */
bool constantsUseConcentrations = false;

/**
 * Whether to ignore non-constant global parameters.
 */
bool ignoreVariableParameters = false;

/**
 * Assignment rules specified by the user.
 */
map<string,string> userRules;

/**
 * Parameters specified by the user.
 */
map<string,double> userParameters;

void importSBMLModel(Hdf5File * lmFile, string sbmlFilename);
bool importSBMLModelL3V1(ReactionModel * lmModel, Model * sbmlModel);
bool importSBMLModelL3V1Kinetics(Reaction * reaction, uint reactionIndex, KineticLaw * kinetics, ndarray<uint>& T, ndarray<double>& K, ndarray<uint>& D, map<string,uint> & speciesIndices, uint numberReactions, map<string,double> & globalParameterValues, map<string,ASTNode_t*>& globalExpressions);
bool matchL3V1KineticsWithPropensityFunction(Reaction * reaction, uint reactionIndex, KineticLaw * kinetics, ndarray<uint>& T, ndarray<double>& K, ndarray<uint>& D, map<string,uint>& speciesIndices, map<string,double>& parameterValues, map<string,ASTNode_t*>& globalExpressions);
bool createReactionModelEntry(ASTNode_t* sourceFormula, ASTNode_t* propensityFormula, uint reactionIndex, ndarray<double>& K, ndarray<uint>& D, map<string,uint>& speciesIndices);

// Allocate the profile space.
PROF_ALLOC;

int main(int argc, char** argv)
{	
	try
	{
		printCopyright(argc, argv);
		parseArguments(argc, argv);

        lm::ClassFactory::getInstance().printRegisteredClasses();
		
		if (function == "help")
		{
			printUsage(argc, argv);
		}
		else if (function == "version")
		{
		}
		else if (function == "import")
		{
            // Make sure the input file exists.
            struct stat fileStats;
            if (stat(inputFilename.c_str(), &fileStats) != 0)
            {
                throw Exception("The specified SBML file does not exist", inputFilename.c_str());
            }

		    // If the output file doesn't exist, create it.
		    if (stat(outputFilename.c_str(), &fileStats) != 0)
		    {
		        Hdf5File::create(outputFilename);
		    }

			// Open the file.
		    Hdf5File outputFile(outputFilename);

		    // Open the sbml file.
		    importSBMLModel(&outputFile, inputFilename);

		    // Close the file.
		    outputFile.close();

            Print::printf(Print::INFO, "Import completed successfully.");
		}
		else
		{
			throw lm::CommandLineArgumentException("unknown function.");
		}
		return 0;
	}
    catch (lm::CommandLineArgumentException& e)
    {
        std::cout << "Invalid command line argument: " << e.what() << std::endl << std::endl;
        printUsage(argc, argv);
    }
    catch (Exception& e)
    {
        std::cout << "Exception during execution: " << e.what() << std::endl;
    }
    catch (std::exception& e)
    {
        std::cout << "Exception during execution: " << e.what() << std::endl;
    }
    catch (...)
    {
        std::cout << "Unknown Exception during execution." << std::endl;
    }
    return -1;
}

void importSBMLModel(Hdf5File * lmFile, string sbmlFilename)
{
    // Read in the SBML document.
    SBMLReader reader;
    std::auto_ptr<SBMLDocument> sbmlDocument(reader.readSBML(sbmlFilename));
    if (sbmlDocument->getNumErrors() > 0)
    {
        bool criticalErrrors = false;
        for (int i=0; i<sbmlDocument->getNumErrors(); i++)
        {
            const SBMLError* error = sbmlDocument->getError(i);
            if (error->getSeverity() >= LIBSBML_SEV_ERROR)
                criticalErrrors = true;
        }


        Print::printf(Print::WARNING,"Problems detected while parsing the SBML file %s",sbmlFilename.c_str());
        Print::printf(Print::WARNING,"-----------------------------------");
        sbmlDocument->printErrors(std::cout);
        Print::printf(Print::WARNING,"-----------------------------------");

        // If there were critical errors and we are not ignoring exceptions, stop.
        if (criticalErrrors && !ignoreErrors)
            throw Exception("There were critical errors detected while parsing the SBML file. Either fix the errors or execute the command again with the --ignore-errors flag set.");
    }

    // Create the appropriate importer instance.
    robertslab::sbml::SBMLImporterL3V1 *importer = NULL;
    if (sbmlDocument->getLevel() == 3 && sbmlDocument->getVersion() == 1)
    {
        if (isCopasi)
            importer = new robertslab::sbml::SBMLImporterL3V1COPASI();
        else
            importer = new robertslab::sbml::SBMLImporterL3V1();
    }
    else
    {
        throw Exception("Unsupported SBML format", sbmlDocument->getLevel(), sbmlDocument->getVersion());
    }

    // Import the document.
    importer->setOptions(constantsUseConcentrations, verbose, reallyVerbose, ignoreErrors, ignoreUnmatchedReactions, ignoreVariableParameters);
    bool success = importer->import(sbmlDocument.get(), userParameters, userRules);
    if (success || ignoreErrors || ignoreUnmatchedReactions)
    {
        if (ignoreErrors || ignoreUnmatchedReactions)
            Print::printf(Print::WARNING, "There were errors during the import process, but the output file was still generated. However, the output file is likely not entirely correct.");
        lmFile->setReactionModel(importer->getReactionModel());
    }
    else
    {
        Print::printf(Print::ERROR, "There were errors during the import process, the output file was not generated.  Either fix the errors or execute the command again with the --ignore-errors flag set.");
    }
}


/**
 * This function prints the copyright notice.
 */
void printCopyright(int argc, char** argv)
{
    std::cout << argv[0] << " v" << VERSION_NUM << " build " << BUILD_TIMESTAMP << std::endl;
    std::cout << "Copyright (C) " << COPYRIGHT_DATE_JHU << " Roberts Group, Johns Hopkins University." << std::endl;
    std::cout << "Copyright (C) " << COPYRIGHT_DATE << " Luthey-Schulten Group, University of Illinois at Urbana-Champaign." << std::endl << std::endl;
	std::cout << std::endl;
}

/**
 * Parses the command line arguments.
 */
void parseArguments(int argc, char** argv)
{
    function = "import";
    outputFilename = "";
    inputFilename = "";

    // Parse any arguments.
    for (int i=1; i<argc; i++)
    {
        char *option = argv[i];
        while (*option == ' ') option++;
        
        // See if the user is trying to get help.
        if (strcmp(option, "-h") == 0 || strcmp(option, "--help") == 0) {
        	function = "help";
        	break;
        }
        
        // See if the user is trying to get the version info.
        else if (strcmp(option, "-v") == 0 || strcmp(option, "--version") == 0) {
        	function = "version";
        	break;
        }

        // See if the user wants to ignore any error and continue processing.
        else if (strcmp(option, "--ignore-errors") == 0)
        {
            ignoreErrors = true;
        }

        // See if the user wants to ignore any unmatched reactions and continue processing.
        else if (strcmp(option, "--ignore-unmatched") == 0)
        {
            ignoreUnmatchedReactions = true;
        }

        // See if the user wants to ignore variable parameters.
        else if (strcmp(option, "--ignore-variable-parameters") == 0)
        {
            ignoreVariableParameters = true;
        }

        // See if the user wants verbose messages.
        else if (strcmp(option, "--verbose") == 0)
        {
            if (!verbose)
                verbose = true;
            else
                reallyVerbose = true;
        }

        // See if the user is trying to import a sbml file that was originally exported by COPASI.
        else if (strcmp(option, "--copasi") == 0)
        {
            isCopasi = true;
            constantsUseConcentrations = true;
        }

        // See if the user is trying to import a sbml file that was originally exported by COPASI.
        else if (strcmp(option, "--constants-use-amounts") == 0)
        {
            constantsUseConcentrations = false;
        }

        // See if the user is trying to import a sbml file that was originally exported by COPASI.
        else if (strcmp(option, "--constants-use-concentration") == 0)
        {
            constantsUseConcentrations = true;
        }

        // See if the user is trying to specify an expression rule key value pair.
        else if (strncmp(option, "--rule:", strlen("--rule:"))==0 && strstr(option, "=") != NULL)
        {
            option += strlen("--rule:");
            char * separator=strstr(option, "=");
            if (separator != NULL && separator > option && strlen(separator) > 1)
            {
                string key(option, separator-option);
                userRules[key] = string(separator+1);
            }
        }

        // See if the user is trying to specify a parameter key value pair.
        else if (strstr(option, "=") != NULL)
        {
            char * separator=strstr(option, "=");
            if (separator != NULL && separator > option && strlen(separator) > 1)
            {
                string key(option, separator-option);
                userParameters[key] = atof(separator+1);
            }
        }

        // See if the user is trying to specify an output filename.
        else if (outputFilename == "")
        {
            outputFilename = option;
        }

        // See if the user is trying to specify an input filename
        else if (inputFilename == "")
        {
            inputFilename = option;
        }

        // This must be an invalid option.
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
    std::cout << "Usage: " << argv[0] << " [OPTIONS] lm_filename sbml_filename" << std::endl; // TODO: uncomment rest of line when userParameterValues is implemented (see below) // (simulation_parameter_key=value)+" << std::endl;
	std::cout << std::endl;
    std::cout << "OPTIONS" << std::endl;
    std::cout << "  key=double_value                Specify a new or override an existing global parameter in the SBML file." << std::endl;
    std::cout << "  --rule:key=string_value         Specify a new or override an existing assignment rule in the SBML file." << std::endl;
    std::cout << "  --ignore-errors                 Use this option to ignore any errors in the SBML file and attempt to import it." << std::endl;
    std::cout << "  --ignore-unmatched              Use this option to ignore any unmatched reactions during the import." << std::endl;
    std::cout << "  --ignore-variable-parameters    Use this option to ignore variable parameters." << std::endl;
    std::cout << "  --verbose                       This option causes more detailed error information to be printed. Use twice for debugging messages." << std::endl;
    std::cout << "  --copasi                        Specify that the SBML file was generated by COPASI and needs to be processed accordingly." << std::endl;
    std::cout << "  --constants-use-amounts         Specify that rate constants are given in terms of particle counts." << std::endl;
    std::cout << "  --constants-use-concentration   Specify that rate constants are given in terms of concentrations." << std::endl;

    std::cout << std::endl;
}

