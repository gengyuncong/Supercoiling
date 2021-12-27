/*
 * Copyright 2017-2019 Johns Hopkins University
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * Developed by: Roberts Group
 *               Johns Hopkins University
 *               http://biophysics.jhu.edu/roberts/
 *
 * Author(s): Elijah Roberts
 */

#include <cmath>
#include <fstream>
#include <iostream>
#include <list>
#include <regex>
#include <set>

#include <sbml/math/ASTNode.h>
#include <sbml/math/L3Parser.h>

#include "lm/Exceptions.h"
#include "lm/Math.h"
#include "lm/Print.h"
#include "lm/Types.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/me/PropensityFunction.h"
#include "robertslab/sbml/ASTHelper.h"
#include "robertslab/bngl/BNGLImporter.h"
#include "robertslab/bngl/InstanceDefinitions.h"
#include "robertslab/bngl/PatternDefinitions.h"
#include "robertslab/bngl/TypeDefinitions.h"
#include "robertslab/graph/Graph.h"

using std::list;
using std::pair;
using std::regex;
using std::regex_token_iterator;
using std::set;
using lm::Exception;
using lm::Print;

using robertslab::sbml::ASTHelper;
using robertslab::graph::GraphMapping;

namespace robertslab {
namespace bngl {

BNGLImporter::BNGLImporter()
:maxRounds(20),propensityFunctions(NULL),constantsUseConcentrations(false),verbose(false),reallyVerbose(false),sorted(false),includeDuplicates(false),allImportStepsSuccessful(true),numberSpecies(0),C(NULL),numberReactions(0),S(NULL),T(NULL),K(NULL),D(NULL)
{
    propensityFunctions = new lm::me::PropensityFunctionFactory();
}

BNGLImporter::~BNGLImporter()
{
    if (propensityFunctions != NULL) delete propensityFunctions; propensityFunctions = NULL;
    if (C != NULL) delete C; C = nullptr;
    if (S != NULL) delete S; S = NULL;
    if (T != NULL) delete T; T = NULL;
    if (K != NULL) delete K; K = NULL;
    if (D != NULL) delete D; D = NULL;
}

void BNGLImporter::setOptions(bool constantsUseConcentrations, bool verbose, bool reallyVerbose, bool sorted, bool includeDuplicates)
{
    this->constantsUseConcentrations = constantsUseConcentrations;
    this->verbose = verbose;
    this->reallyVerbose = reallyVerbose;
    this->sorted = sorted;
    this->includeDuplicates = includeDuplicates;
}

bool BNGLImporter::import(string filename, map<string,double> userParameters)
{
    Print::printf(Print::INFO, "Processing BNGL document: %s", filename.c_str());

    // Initialize the parameters with any user specified parameters.
    for (map<string,double>::iterator it=userParameters.begin(); it!=userParameters.end(); it++)
    {
        if (parameters.count(it->first) == 0)
        {
            parameters[it->first] = it->second;
            Print::printf(Print::INFO, "Added user defined parameter %s: %e", it->first.c_str(), parameters[it->first]);
        }
    }

    // Open the file.
    std::ifstream input(filename, std::ifstream::in);

    regex beginPattern("^begin\\s+(.+)$");
    regex endPattern("^end\\s+(.+)$");
    std::smatch match, match2;
    string section = "";
    list<string> sectionLines;
    int lineNumber=0;
    string continuedLine = "";
    while (!input.eof())
    {
        lineNumber++;
        string line;
        std::getline(input,line);

        // If the line end with a backslash, it is continued on the next line.
        if (line.back() == '\\')
        {
            continuedLine += line.substr(0,line.size()-1);
            continue;
        }
        else if (continuedLine != "")
        {
            line = continuedLine + line;
            continuedLine = "";
        }

        // Strip any comments.
        if (line.find_first_of('#') != string::npos) line = line.substr(0,line.find_first_of('#'));

        // Skip the line if it is blank.
        if (line.find_first_not_of(" \t\r\n") == string::npos) continue;

        // Strip any leading or trailing whitespace.
        size_t start = line.find_first_not_of(" \t\r\n");
        size_t end = line.find_last_not_of(" \t\r\n");
        line = line.substr(start,end-start+1);

        // Strip any leading line numbers.
        regex lineNumberPattern("^(?:\\d+\\s+)?(\\S.*)$");
        regex nullReactionPattern("^0\\s*\\<?\\-\\>?.*$");
        if (std::regex_match(line, match, lineNumberPattern) && match.size() == 2 && !std::regex_match(line, match2, nullReactionPattern))
            line = match[1].str();

        // Strip any leading line labels.
        regex lineLabelPattern("^(?:\\S+\\:\\s+)?(\\S.*)$");
        if (std::regex_match(line, match, lineLabelPattern) && match.size() == 2)
            line = match[1].str();

        // Check for section blocks.
        if (std::regex_match(line, match, beginPattern) && match.size() == 2 && match[1].str() != "model")
        {
            if (section == "")
            {
                section = match[1].str();
                sectionLines.clear();
            }
            else
            {
                Print::printf(Print::ERROR, "invalid begin section on line %d: %s",lineNumber,line.c_str());
                allImportStepsSuccessful = false;
            }
        }
        else if (std::regex_match(line, match, endPattern) && match.size() == 2 && match[1].str() != "model")
        {
            if (section == match[1].str())
            {
                // Process the section.
                if (section == "parameters")
                {
                    parseParameters(sectionLines);
                }
                else if (section == "molecule types")
                {
                    parseMoleculeTypes(sectionLines);
                }
                else if (section == "species" || section == "seed species" )
                {
                    parseSeedSpecies(sectionLines);
                }
                else if (section == "reaction rules")
                {
                    // Add any molecules that are in the species list but not in the molecule types before parsing the reactions.
                    supplementMoleculeTypesFromSeedSpecies();

                    // Parse the reactions.
                    parseReactions(sectionLines);
                }
                else
                {
                    Print::printf(Print::INFO, "Ignoring unsupported block type: %s", section.c_str());
                }

                // Reset the section.
                section = "";
            }
            else
            {
                Print::printf(Print::ERROR, "invalid end section on line %d: %s",lineNumber,line.c_str());
                allImportStepsSuccessful = false;
            }
        }
        else if (section != "" && line != "")
        {
            sectionLines.push_back(line);
        }

    }

    // If we imported the data correctly, process it.
    if (allImportStepsSuccessful)
    {
        Print::printf(Print::INFO, "Processing BNGL model.");
        processModel();
    }

    return allImportStepsSuccessful;
}

lm::input::ReactionModel* BNGLImporter::getReactionModel()
{
    return &lmModel;
}

string BNGLImporter::getDescription()
{
    return "BNGL Importer";
}

void BNGLImporter::parseParameters(list<string>& lines)
{
    Print::printf(Print::INFO, "Parsing parameters block.");

    regex parameterPattern("^(\\S+)\\s+(\\S+)$");
    std::smatch match;
    for (list<string>::iterator it=lines.begin(); it != lines.end(); it++)
    {
        string line = *it;
        if (std::regex_match(line, match, parameterPattern) && match.size() == 3)
        {
            // Get the key and the expression.
            string key = match[1].str();
            string expression = match[2].str();

            // If we got to a numeric expression, save it.
            double value;
            if (evaluteExpression(expression, value))
            {
                if (parameters.count(key) == 0)
                {
                    parameters[key] = value;
                    Print::printf(Print::INFO, "Added parameter %s: %e", match[1].str().c_str(), value);
                }
                else
                {
                    Print::printf(Print::INFO, "Skipped duplicate parameter %s: %e", key.c_str(), value);
                }
            }
            else
            {
                Print::printf(Print::ERROR, "Could not simplify parameter %s: %s", key.c_str(), expression.c_str());
                allImportStepsSuccessful = false;
            }
        }
        else
        {
            Print::printf(Print::WARNING, "Could not parse parameter from block: \"%s\"", line.c_str());
        }
    }
}

void BNGLImporter::parseMoleculeTypes(list<string>& lines)
{
    Print::printf(Print::INFO, "Parsing molecule types block.");

    for (list<string>::iterator it=lines.begin(); it != lines.end(); it++)
    {
        string line = *it;
        MoleculeClass* molecule = new MoleculeClass(line);
        if (molecule->isValid())
        {
            moleculeTypes[molecule->getName()] = molecule;
            Print::printf(Print::INFO, "Added molecule definition: %s", molecule->getString().c_str());
        }
        else
        {
            Print::printf(Print::WARNING, "Could not parse molecule from block: \"%s\"", line.c_str());
        }
    }
}

void BNGLImporter::parseSeedSpecies(list<string>& lines)
{
    Print::printf(Print::INFO, "Parsing species block.");

    // Add the round to the complex species and reactions list.
    complexSpecies.push_back(vector<ComplexInstance*>());
    reactions.push_back(vector<ReactionInstance*>());

    regex parameterPattern("^(\\S+)\\s+(\\S+)$");
    std::smatch match;
    for (list<string>::iterator it=lines.begin(); it != lines.end(); it++)
    {
        string line = *it;
        if (std::regex_match(line, match, parameterPattern) && match.size() == 3)
        {
            // Get the count.
            string countExpression = match[2].str();
            double count=0.0;
            bool isValidCount = evaluteExpression(countExpression,count);

            // Get the complex.
            string complexInstanceString = match[1].str();
            ComplexInstance* complex = new ComplexInstance(complexInstanceString, count);

            // Make sure we processed a valid record.
            if (complex->isValid() && isValidCount)
            {
                complexSpecies[0].push_back(complex);
                Print::printf(Print::INFO, "Added initial count %s", complex->getString(true,true).c_str());
            }
            else if (!complex->isValid())
            {
                Print::printf(Print::ERROR, "Could not simplify complex %s", complexInstanceString.c_str());
                allImportStepsSuccessful = false;
            }
            else if (!isValidCount)
            {
                Print::printf(Print::ERROR, "Could not simplify initial count: \"%s\"", countExpression.c_str());
                allImportStepsSuccessful = false;
            }
        }
        else
        {
            Print::printf(Print::WARNING, "Could not parse initial count from block: \"%s\"", line.c_str());
        }
    }
}

void BNGLImporter::parseReactions(list<string>& lines)
{
    Print::printf(Print::INFO, "Parsing reactions block.");
    regex reversibleReactionPattern("^([^<->]+)\\s*<->\\s*([^<->]+)\\s+(\\S+),\\s*(\\S+)$");
    regex irreversibleReactionPattern("^([^<->]+)\\s*->\\s*([^<->]+)\\s+(\\S+)$");
    std::smatch match;
    for (list<string>::iterator it=lines.begin(); it != lines.end(); it++)
    {
        string line = *it;
        if (std::regex_match(line, match, reversibleReactionPattern) && match.size() == 5)
        {
            // Get the reaction rates.
            string rateFString = match[3].str();
            double rateF=0.0;
            bool isValidRateF = evaluteExpression(rateFString,rateF);
            string rateRString = match[4].str();
            double rateR=0.0;
            bool isValidRateR = evaluteExpression(rateRString,rateR);

            // Get the lhs and rhs of the equation.
            string lhsString = match[1].str();
            string rhsString = match[2].str();
            ReactionPattern* reactionF = new ReactionPattern(lhsString, rhsString, rateF, moleculeTypes);
            ReactionPattern* reactionR = new ReactionPattern(rhsString, lhsString, rateR, moleculeTypes);

            // Make sure we processed a valid record.
            if (isValidRateF && isValidRateR && reactionF->isValid() && reactionR->isValid())
            {
                reactionPatterns.push_back(reactionF);
                reactionPatterns.push_back(reactionR);
                Print::printf(Print::INFO, "Added forward reaction pattern %s", reactionF->getString(true,true).c_str());
                if (reallyVerbose) Print::printf(Print::INFO, "Reaction substrate to product mapping:\n%s", reactionF->getSubstrateToProductMapping().getString().c_str());
                Print::printf(Print::INFO, "Added reverse reaction pattern %s", reactionR->getString(true,true).c_str());
                if (reallyVerbose) Print::printf(Print::INFO, "Reaction substrate to product mapping:\n%s", reactionR->getSubstrateToProductMapping().getString().c_str());
            }
            else if (!reactionF->isValid() || !reactionR->isValid())
            {
                Print::printf(Print::ERROR, "Could not simplify reaction: \"%s\" \"%s\"", lhsString.c_str(), rhsString.c_str());
                allImportStepsSuccessful = false;
            }
            else if (!isValidRateF)
            {
                Print::printf(Print::ERROR, "Could not simplify reaction rate: \"%s\"", rateFString.c_str());
                allImportStepsSuccessful = false;
            }
            else if (!isValidRateR)
            {
                Print::printf(Print::ERROR, "Could not simplify reaction rate: \"%s\"", rateRString.c_str());
                allImportStepsSuccessful = false;
            }
        }
        else if (std::regex_match(line, match, irreversibleReactionPattern) && match.size() == 4)
        {
            // Get the reaction rate.
            string rateString = match[3].str();
            double rate=0.0;
            bool isValidRate = evaluteExpression(rateString,rate);

            // Get the lhs and rhs of the equation.
            string lhsString = match[1].str();
            string rhsString = match[2].str();
            ReactionPattern* reaction = new ReactionPattern(lhsString, rhsString, rate, moleculeTypes);

            // Make sure we processed a valid record.
            if (reaction->isValid() && isValidRate)
            {
                reactionPatterns.push_back(reaction);
                Print::printf(Print::INFO, "Added irreversible reaction pattern %s", reaction->getString(true,true).c_str());
                if (reallyVerbose) Print::printf(Print::INFO, "Reaction substrate to product mapping:\n%s", reaction->getSubstrateToProductMapping().getString().c_str());
            }
            else if (!reaction->isValid())
            {
                Print::printf(Print::ERROR, "Could not simplify left reaction: \"%s\" \"%s\"", lhsString.c_str(), rhsString.c_str());
                allImportStepsSuccessful = false;
            }
            else if (!isValidRate)
            {
                Print::printf(Print::ERROR, "Could not simplify reaction rate: \"%s\"", rateString.c_str());
                allImportStepsSuccessful = false;
            }
        }
        else
        {
            Print::printf(Print::ERROR, "Could not parse reaction from block: \"%s\"", line.c_str());
            allImportStepsSuccessful = false;
        }
    }
}

bool BNGLImporter::evaluteExpression(string expression, double& value)
{
    // Simplify the formula.
    ASTNode_t* formula = SBML_parseL3Formula(expression.c_str());
    ASTHelper::substituteASTParameters(formula, parameters);
    ASTHelper::simplifyASTExpression(formula);

    // If we are printing debug info, print the AST tree.
    if (reallyVerbose)
    {
        Print::printf(Print::INFO, "Simplified expression %s as:", expression.c_str());
        ASTHelper::printASTNode(formula);
    }

    // If we got to a numeric expression, return it.
    if (ASTHelper::isNumeric(formula))
    {
        value = ASTHelper::getNumericValue(formula);
        return true;
    }
    else
    {
        value = ASTHelper::evaluateASTOperator(formula);
        return true;
    }
}

void BNGLImporter::processModel()
{
    // Figure out the list of atomic species that we need.
    //enumerateMoleculeSpecies();

    // Go through the reaction patterns iteratively until the species have converged.
    bool complexSpeciesConverged=false;
    for (int round=1; round<maxRounds; round++)
    {
        if (reallyVerbose) Print::printf(Print::INFO, "** Processing reactions for round %d **", round);
        if (!processReactions(round))
        {
            if (reallyVerbose) Print::printf(Print::INFO, "** Finished all rounds **", round);
            complexSpeciesConverged = true;
            break;
        }
        if (reallyVerbose) Print::printf(Print::INFO, "** Finished round **", round);
    }


    // Combine the species from each round into one list.
    for (size_t i=0; i<complexSpecies.size(); i++)
    {
        for (size_t j=0; j<complexSpecies[i].size(); j++)
        {
            if (sorted)
            {
                bool added = false;
                for (auto it=allComplexSpecies.begin(); it != allComplexSpecies.end(); it++)
                {
//                    Print::printf(Print::INFO, "Comparing %s to %s = %d", complexSpecies[i][j]->getString().c_str(), (*it)->getString());
                    if (complexSpecies[i][j]->getString(false,false).compare((*it)->getString(false,false)) <= 0)
                    {
                        allComplexSpecies.insert(it, complexSpecies[i][j]);
                        added = true;
                        break;
                    }
                }
                if (!added) allComplexSpecies.push_back(complexSpecies[i][j]);
            }
            else
            {
                allComplexSpecies.push_back(complexSpecies[i][j]);
            }
        }
    }

    // Combine the reactions from each round into one list.
    for (size_t i=0; i<reactions.size(); i++)
    {
        for (size_t j=0; j<reactions[i].size(); j++)
        {
            if (sorted)
            {
                bool added = false;
                for (auto it=allReactions.begin(); it != allReactions.end(); it++)
                {
                    if (reactions[i][j]->getSubstrate()->getNumberComplexes() < (*it)->getSubstrate()->getNumberComplexes() || (reactions[i][j]->getSubstrate()->getNumberComplexes() == (*it)->getSubstrate()->getNumberComplexes() && reactions[i][j]->getString(false,false).compare((*it)->getString(false,false)) <= 0))
                    {
                        allReactions.insert(it, reactions[i][j]);
                        added = true;
                        break;
                    }
                }
                if (!added) allReactions.push_back(reactions[i][j]);
            }
            else
            {
                allReactions.push_back(reactions[i][j]);
            }
        }
    }

    Print::printf(Print::INFO, "Found %d total complex species.", allComplexSpecies.size());
    Print::printf(Print::INFO, "Found %d total reactions.", allReactions.size());


    // Print debugging information, if necessary.
    if (verbose)
    {
        for (size_t i=0; i<allComplexSpecies.size(); i++)
            Print::printf(Print::INFO, "Species %d:  %s", i, allComplexSpecies[i]->getString().c_str());
        for (size_t i=0; i<allReactions.size(); i++)
            Print::printf(Print::INFO, "Reaction %d: %s", i, allReactions[i]->getString().c_str());
    }

    // Build the species model.
    buildSpeciesModel();

    // Build the reaction model.
    buildReactionModel();

    // Print the matrices.
    if (verbose)
    {
        Print::printf(Print::INFO, "Species count matrix was:");
        C->print(); printf("\n");
        Print::printf(Print::INFO, "Reaction type matrix was:");
        T->print(); printf("\n");
        Print::printf(Print::INFO, "Rate constant matrix was:");
        K->print(); printf("\n");
        Print::printf(Print::INFO, "Stoichiometry matrix was:");
        S->print(); printf("\n");
        Print::printf(Print::INFO, "Dependency matrix was:");
        D->print(); printf("\n");
    }

    // Build the lm model.
    buildLMModel();

    // If the species did not converge, print a warning at the end.
    if (!complexSpeciesConverged) Print::printf(Print::WARNING, "Complex species did not converge after %d rounds. Either ensure that you have enough species depth or run again with more rounds.");
}

void BNGLImporter::supplementMoleculeTypesFromSeedSpecies()
{
    // Go through the list of complexes with initial counts.
    for (size_t i=0; i<complexSpecies.size(); i++)
    {
        ComplexInstance* complex = complexSpecies[0][i];

        // Go through the molecules in the complex.
        for (int j=0; j<complex->getNumberMolecules(); j++)
        {
            MoleculeInstance* molecule = complex->getMolecule(j);

            // If we don't already have a molecule type with this name, add one.
            if (moleculeTypes.count(molecule->getName()) == 0)
            {
                MoleculeClass* moleculeClass = new MoleculeClass(molecule->getString());
                moleculeTypes[moleculeClass->getName()] = moleculeClass;
                Print::printf(Print::WARNING, "Added inferred molecule definition: %s", moleculeClass->getString().c_str());
            }
        }
    }
}

void BNGLImporter::enumerateMoleculeSpecies()
{
    // We need one species for each combination of states for every molecule.
    for (auto it=moleculeTypes.begin(); it != moleculeTypes.end(); it++)
    {
        MoleculeClass* molecule = it->second;
        vector<string> stateCombinations = molecule->getStateCombinations();
        for (size_t j=0; j<stateCombinations.size(); j++)
        {
            moleculeSpecies.push_back(new MoleculeInstance(molecule->getName()+"("+stateCombinations[j]+")"));
        }
        Print::printf(Print::INFO, "Added species to represent possible states for molecule %s: %d species", molecule->getName().c_str(), stateCombinations.size());
    }
    Print::printf(Print::INFO, "Added %d total molecular species.", moleculeSpecies.size());

    // Print debugging information, if necessary.
    if (reallyVerbose)
    {
        for (size_t i=0; i<moleculeSpecies.size(); i++)
        {
            Print::printf(Print::INFO, "%s", moleculeSpecies[i]->getString().c_str());
        }
    }
}

bool BNGLImporter::processReactions(int round)
{
    // And the round to the complex species and reactions list.
    if (complexSpecies.size() != round) throw Exception("inconsistent complex species list size",round,complexSpecies.size());
    if (reactions.size() != round) throw Exception("inconsistent reaction list size",round,reactions.size());
    complexSpecies.push_back(vector<ComplexInstance*>());
    reactions.push_back(vector<ReactionInstance*>());

    // Loop through each reaction.
    for (size_t i=0; i<reactionPatterns.size(); i++)
    {
        ReactionPattern* reactionPattern = reactionPatterns[i];
        processReaction(round, reactionPattern);
    }

    Print::printf(Print::INFO, "Finished round round %d: %d new complex species, %d new reactions", round, complexSpecies[round].size(), reactions[round].size());

    // Return whether or not we added any new complexes this round.
    return complexSpecies[round].size() != 0;
}

void BNGLImporter::processReaction(int round, ReactionPattern* reactionPattern)
{
    if (reallyVerbose) Print::printf(Print::INFO, "**** Processing reaction %s ****", reactionPattern->getString(false).c_str());

    // Process the reaction according to its order.
    if (reactionPattern->getSubstrates()->getNumberReactants() == 0 || (reactionPattern->getSubstrates()->getNumberReactants() == 1 && reactionPattern->getSubstrates()->getReactant(0)->isNull()))
        processReactionZerothOrder(round, reactionPattern);
    else if (reactionPattern->getSubstrates()->getNumberReactants() == 1)
        processReactionFirstOrder(round, reactionPattern);
    else if (reactionPattern->getSubstrates()->getNumberReactants() == 2)
        processReactionSecondOrder(round, reactionPattern);
    else if (reactionPattern->getSubstrates()->getNumberReactants() == 3)
        processReactionThirdOrder(round, reactionPattern);
    else if (reactionPattern->getSubstrates()->getNumberReactants() == 4)
        processReactionFourthOrder(round, reactionPattern);
    else
        throw Exception("unsupported reaction order",reactionPattern->getSubstrates()->getNumberReactants());

    if (reallyVerbose) Print::printf(Print::INFO, "**** Finished reaction ****");
}

void BNGLImporter::processReactionZerothOrder(int round, ReactionPattern* reactionPattern)
{
    // Only process zeroth order reactions in round 1.
    if (round == 1)
    {
        // Get the products.
        ReactantInstance* substrate = new ReactantInstance();
        ReactantInstance* product = new ReactantInstance();

        // Create any new the products, they must be completely defined.
        ReactantPattern* productPattern = reactionPattern->getProducts();
        for (int m=0; m<productPattern->getNumberReactants(); m++)
        {
            // Skip any null products.
            if (productPattern->getReactant(m)->isNull())
                continue;

            // Find any products from the previous round that match the pattern.
            int productMatchesFound=0;
            ComplexInstance* createdProductComplex=NULL;
            for (size_t o=0; o<complexSpecies[round-1].size(); o++)
            {
                ComplexInstance* complex3 = complexSpecies[round-1][o];
                if (complex3->isIsomorphic(productPattern->getReactant(m)))
                {
                    productMatchesFound++;
                    if (createdProductComplex == NULL)
                        createdProductComplex = new ComplexInstance(*complex3);
                }
            }

            // Make sure we found only a single matching instance.
            if (productMatchesFound == 0)
            {
                Print::printf(Print::ERROR, "Could not find a matching complex instance for product pattern %s in reaction %s", productPattern->getReactant(m)->getString().c_str(), reactionPattern->getString(false).c_str());
                throw Exception("invalid input model");
            }
            else if (productMatchesFound > 1)
            {
                Print::printf(Print::ERROR, "Found mutliple (%d) matching complex instances for product pattern %s in reaction %s", productMatchesFound, productPattern->getReactant(m)->getString().c_str(), reactionPattern->getString(false).c_str());
                throw Exception("invalid input model");
            }
            else
            {
                // See if the product was not used in the substrate to product mapping.
                if (reallyVerbose) Print::printf(Print::INFO, "Product %s created", createdProductComplex->getString().c_str());
                product->addComplex(createdProductComplex);
            }
        }

        // Save the new reaction.
        ReactionInstance* reaction = new ReactionInstance(substrate, product, reactionPattern->getRate());
        if (isNewReaction(reaction))
        {
            reactions[round].push_back(reaction);
            if (reallyVerbose) Print::printf(Print::INFO, "Added new reaction %s", reaction->getString().c_str());
        }

        // If the product contains any new species, add them to the list.
        for (int j=0; j<product->getNumberComplexes(); j++)
        {
            if (isNewComplexSpecies(product->getComplex(j)))
                complexSpecies[round].push_back(product->getComplex(j));
        }
    }
}

void BNGLImporter::processReactionFirstOrder(int round, ReactionPattern* reactionPattern)
{
    // Get the substrate pattern.
    ReactantPattern* substratePattern = reactionPattern->getSubstrates();

    // Go through every complex species from the previous round and see if they match the pattern.
    for (size_t i=0; i<complexSpecies[round-1].size(); i++)
    {
        // Get the substrate complex.
        ComplexInstance* complex1 = new ComplexInstance(*complexSpecies[round-1][i]);

        // Get the substrate.
        ReactantInstance* substrate = new ReactantInstance(complex1);

        // Find all the matches to the substrate pattern.
        list<GraphMapping> matches = complex1->findAllIsomorphicSubgraphs(substratePattern);

        // Go through each match.
        for (auto it=matches.begin(); it != matches.end(); it++)
        {
            GraphMapping substrateMapping(substrate, substratePattern);
            substrateMapping.appendMappings(*it);
            if (reallyVerbose) Print::printf(Print::INFO, "Found matching substrate %s", substrate->getString().c_str());

            // Create any transformed products.
            ReactantInstance* product = rewriteSubstrateToProduct(substrate, substrateMapping, substratePattern, reactionPattern->getSubstrateToProductMapping());
            if (reallyVerbose) Print::printf(Print::INFO, "Substrate %s transformed into product %s", substrate->getString().c_str(), product->getString().c_str());

            // Create any new the products, they must be completely defined.
            ReactantPattern* productPattern = reactionPattern->getProducts();
            for (int m=0; m<productPattern->getNumberReactants(); m++)
            {
                // Skip any null products.
                if (productPattern->getReactant(m)->isNull())
                    continue;

                // Skip any products that were transformed in the reaction.
                if (reactionPattern->isProductReactantTransformed(m))
                    continue;

                // Find a product that matches the pattern.
                int productMatchesFound=0;
                ComplexInstance* createdProductComplex=NULL;
                for (int n=0; n<=round-1; n++)
                {
                    for (size_t o=0; o<complexSpecies[n].size(); o++)
                    {
                        ComplexInstance* complex3 = complexSpecies[n][o];
                        if (complex3->isIsomorphic(productPattern->getReactant(m)))
                        {
                            productMatchesFound++;
                            if (createdProductComplex == NULL)
                                createdProductComplex = new ComplexInstance(*complex3);
                        }
                    }
                }

                // Make sure we found only a single matching instance.
                if (productMatchesFound == 0)
                {
                    Print::printf(Print::ERROR, "Could not find a matching complex instance for product pattern %s in reaction %s", productPattern->getReactant(m)->getString().c_str(), reactionPattern->getString(false).c_str());
                    throw Exception("invalid input model");
                }
                else if (productMatchesFound > 1)
                {
                    Print::printf(Print::ERROR, "Found mutliple (%d) matching complex instances for product pattern %s in reaction %s", productMatchesFound, productPattern->getReactant(m)->getString().c_str(), reactionPattern->getString(false).c_str());
                    throw Exception("invalid input model");
                }
                else
                {
                    // See if the product was not used in the substrate to product mapping.
                    if (reallyVerbose) Print::printf(Print::INFO, "Product %s created", createdProductComplex->getString().c_str());
                    product->addComplex(createdProductComplex);
                }
            }

            // Save the new reaction.
            ReactionInstance* reaction = new ReactionInstance(substrate, product, reactionPattern->getRate());
            if (isNewReaction(reaction))
            {
                reactions[round].push_back(reaction);
                if (reallyVerbose) Print::printf(Print::INFO, "Added new reaction %s", reaction->getString().c_str());
            }

            // If the product contains any new species, add them to the list.
            for (int j=0; j<product->getNumberComplexes(); j++)
            {
                if (isNewComplexSpecies(product->getComplex(j)))
                    complexSpecies[round].push_back(product->getComplex(j));
            }
        }
    }
}

void BNGLImporter::processReactionSecondOrder(int round, ReactionPattern* reactionPattern)
{
    // Get the substrate pattern.
    ReactantPattern* substratePattern = reactionPattern->getSubstrates();

    // Go through every pair of complex species.
    for (int i=0; i<=round-1; i++)
    {
        for (size_t j=0; j<complexSpecies[i].size(); j++)
        {
            for (int k=0; k<=round-1; k++)
            {
                for (size_t l=0; l<complexSpecies[k].size(); l++)
                {
                    // Make sure at least one of the pair is from the previous round.
                    if (i != round-1 && k != round-1) continue;

                    // Get the substrate complexes.
                    ComplexInstance* complex1 = new ComplexInstance(*complexSpecies[i][j]);
                    ComplexInstance* complex2 = new ComplexInstance(*complexSpecies[k][l]);

                    // Get the substrate.
                    ReactantInstance* substrate = new ReactantInstance(complex1, complex2);

                    // Find all the matches to the substrate pattern.
                    list<GraphMapping> matches1 = complex1->findAllIsomorphicSubgraphs(substratePattern->getReactant(0));
                    list<GraphMapping> matches2 = complex2->findAllIsomorphicSubgraphs(substratePattern->getReactant(1));

                    // Go through each match.
                    for (auto it1=matches1.begin(); it1 != matches1.end(); it1++)
                    {
                        for (auto it2=matches2.begin(); it2 != matches2.end(); it2++)
                        {
                            GraphMapping substrateMapping(substrate, substratePattern);
                            substrateMapping.appendMappings(*it1);
                            substrateMapping.appendMappings(*it2);
                            if (reallyVerbose) Print::printf(Print::INFO, "Found matching substrate %s", substrate->getString().c_str());

                            // Create any transformed products.
                            ReactantInstance* product = rewriteSubstrateToProduct(substrate, substrateMapping, substratePattern, reactionPattern->getSubstrateToProductMapping());
                            if (reallyVerbose) Print::printf(Print::INFO, "Substrate %s transformed into product %s", substrate->getString().c_str(), product->getString().c_str());

                            // Create any new the products, they must be completely defined.
                            ReactantPattern* productPattern = reactionPattern->getProducts();
                            for (int m=0; m<productPattern->getNumberReactants(); m++)
                            {
                                // Skip any null products.
                                if (productPattern->getReactant(m)->isNull())
                                    continue;

                                // Skip any products that were transformed in the reaction.
                                if (reactionPattern->isProductReactantTransformed(m))
                                    continue;

                                // Find a product that matches the pattern.
                                int productMatchesFound=0;
                                ComplexInstance* createdProductComplex=NULL;
                                for (int n=0; n<=round-1; n++)
                                {
                                    for (size_t o=0; o<complexSpecies[n].size(); o++)
                                    {
                                        ComplexInstance* complex3 = complexSpecies[n][o];
                                        if (complex3->isIsomorphic(productPattern->getReactant(m)))
                                        {
                                            if (reallyVerbose) Print::printf(Print::INFO, "Complex %s matched product pattern %s", complex3->getString().c_str(), productPattern->getReactant(m)->getString().c_str());
                                            productMatchesFound++;
                                            if (createdProductComplex == NULL)
                                                createdProductComplex = new ComplexInstance(*complex3);
                                        }
                                    }
                                }

                                // Make sure we found only a single matching instance.
                                if (productMatchesFound == 0)
                                {
                                    Print::printf(Print::ERROR, "Could not find a matching complex instance for product pattern %s in reaction %s", productPattern->getReactant(m)->getString().c_str(), reactionPattern->getString(false).c_str());
                                    throw Exception("invalid input model");
                                }
                                else if (productMatchesFound > 1)
                                {
                                    Print::printf(Print::ERROR, "Found mutliple (%d) matching complex instances for product pattern %s in reaction %s", productMatchesFound, productPattern->getReactant(m)->getString().c_str(), reactionPattern->getString(false).c_str());
                                    throw Exception("invalid input model");
                                }
                                else
                                {
                                    // See if the product was not used in the substrate to product mapping.
                                    if (reallyVerbose) Print::printf(Print::INFO, "Product %s created", createdProductComplex->getString().c_str());
                                    product->addComplex(createdProductComplex);
                                }
                            }

                            // Save the new reaction.
                            ReactionInstance* reaction = new ReactionInstance(substrate, product, reactionPattern->getRate());
                            if (isNewReaction(reaction))
                            {
                                reactions[round].push_back(reaction);
                                if (reallyVerbose) Print::printf(Print::INFO, "Added new reaction %s", reaction->getString().c_str());
                            }

                            // If the product contains any new species, add them to the list.
                            for (int m=0; m<product->getNumberComplexes(); m++)
                            {
                                if (isNewComplexSpecies(product->getComplex(m)))
                                    complexSpecies[round].push_back(product->getComplex(m));
                            }
                        }
                    }
                }
            }
        }
    }
}

void BNGLImporter::processReactionThirdOrder(int round, ReactionPattern* reactionPattern)
{
    // Get the substrate pattern.
    ReactantPattern* substratePattern = reactionPattern->getSubstrates();

    // Go through every triplet of complex species.
    for (int i=0; i<=round-1; i++)
    {
        for (size_t j=0; j<complexSpecies[i].size(); j++)
        {
            for (int k=0; k<=round-1; k++)
            {
                for (size_t l=0; l<complexSpecies[k].size(); l++)
                {
                    for (int m=0; m<=round-1; m++)
                    {
                        for (size_t n=0; n<complexSpecies[m].size(); n++)
                        {
                            // Make sure at least one of the set is from the previous round.
                            if (i != round-1 && k != round-1 && m != round-1) continue;

                            // Get the substrate complexes.
                            ComplexInstance* complex1 = new ComplexInstance(*complexSpecies[i][j]);
                            ComplexInstance* complex2 = new ComplexInstance(*complexSpecies[k][l]);
                            ComplexInstance* complex3 = new ComplexInstance(*complexSpecies[m][n]);

                            // Get the substrate.
                            ReactantInstance* substrate = new ReactantInstance(complex1, complex2, complex3);

                            // Find all the matches to the substrate pattern.
                            list<GraphMapping> matches1 = complex1->findAllIsomorphicSubgraphs(substratePattern->getReactant(0));
                            list<GraphMapping> matches2 = complex2->findAllIsomorphicSubgraphs(substratePattern->getReactant(1));
                            list<GraphMapping> matches3 = complex3->findAllIsomorphicSubgraphs(substratePattern->getReactant(2));

                            // Go through each match.
                            for (auto it1=matches1.begin(); it1 != matches1.end(); it1++)
                            {
                                for (auto it2=matches2.begin(); it2 != matches2.end(); it2++)
                                {
                                    for (auto it3=matches3.begin(); it3 != matches3.end(); it3++)
                                    {
                                        GraphMapping substrateMapping(substrate, substratePattern);
                                        substrateMapping.appendMappings(*it1);
                                        substrateMapping.appendMappings(*it2);
                                        substrateMapping.appendMappings(*it3);
                                        if (reallyVerbose) Print::printf(Print::INFO, "Found matching substrate %s", substrate->getString().c_str());

                                        // Create any transformed products.
                                        ReactantInstance* product = rewriteSubstrateToProduct(substrate, substrateMapping, substratePattern, reactionPattern->getSubstrateToProductMapping());
                                        if (reallyVerbose) Print::printf(Print::INFO, "Substrate %s transformed into product %s", substrate->getString().c_str(), product->getString().c_str());

                                        // Create any new the products, they must be completely defined.
                                        ReactantPattern* productPattern = reactionPattern->getProducts();
                                        for (int a=0; a<productPattern->getNumberReactants(); a++)
                                        {
                                            // Skip any null products.
                                            if (productPattern->getReactant(a)->isNull())
                                                continue;

                                            // Skip any products that were transformed in the reaction.
                                            if (reactionPattern->isProductReactantTransformed(a))
                                                continue;

                                            // Find a product that matches the pattern.
                                            int productMatchesFound=0;
                                            ComplexInstance* createdProductComplex=NULL;
                                            for (int b=0; b<=round-1; b++)
                                            {
                                                for (size_t c=0; c<complexSpecies[b].size(); c++)
                                                {
                                                    ComplexInstance* complexToCheck = complexSpecies[b][c];
                                                    if (complexToCheck->isIsomorphic(productPattern->getReactant(a)))
                                                    {
                                                        if (reallyVerbose) Print::printf(Print::INFO, "Complex %s matched product pattern %s", complexToCheck->getString().c_str(), productPattern->getReactant(a)->getString().c_str());
                                                        productMatchesFound++;
                                                        if (createdProductComplex == NULL)
                                                            createdProductComplex = new ComplexInstance(*complexToCheck);
                                                    }
                                                }
                                            }

                                            // Make sure we found only a single matching instance.
                                            if (productMatchesFound == 0)
                                            {
                                                Print::printf(Print::ERROR, "Could not find a matching complex instance for product pattern %s in reaction %s", productPattern->getReactant(a)->getString().c_str(), reactionPattern->getString(false).c_str());
                                                throw Exception("invalid input model");
                                            }
                                            else if (productMatchesFound > 1)
                                            {
                                                Print::printf(Print::ERROR, "Found mutliple (%d) matching complex instances for product pattern %s in reaction %s", productMatchesFound, productPattern->getReactant(a)->getString().c_str(), reactionPattern->getString(false).c_str());
                                                throw Exception("invalid input model");
                                            }
                                            else
                                            {
                                                // See if the product was not used in the substrate to product mapping.
                                                if (reallyVerbose) Print::printf(Print::INFO, "Product %s created", createdProductComplex->getString().c_str());
                                                product->addComplex(createdProductComplex);
                                            }
                                        }

                                        // Save the new reaction.
                                        ReactionInstance* reaction = new ReactionInstance(substrate, product, reactionPattern->getRate());
                                        if (isNewReaction(reaction))
                                        {
                                            reactions[round].push_back(reaction);
                                            if (reallyVerbose) Print::printf(Print::INFO, "Added new reaction %s", reaction->getString().c_str());
                                        }

                                        // If the product contains any new species, add them to the list.
                                        for (int a=0; a<product->getNumberComplexes(); a++)
                                        {
                                            if (isNewComplexSpecies(product->getComplex(a)))
                                                complexSpecies[round].push_back(product->getComplex(a));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void BNGLImporter::processReactionFourthOrder(int round, ReactionPattern* reactionPattern)
{
    // Get the substrate pattern.
    ReactantPattern* substratePattern = reactionPattern->getSubstrates();

    // Go through every quadruplet of complex species.
    for (int i=0; i<=round-1; i++)
    {
        for (size_t j=0; j<complexSpecies[i].size(); j++)
        {
            for (int k=0; k<=round-1; k++)
            {
                for (size_t l=0; l<complexSpecies[k].size(); l++)
                {
                    for (int m=0; m<=round-1; m++)
                    {
                        for (size_t n=0; n<complexSpecies[m].size(); n++)
                        {
                            for (int o=0; o<=round-1; o++)
                            {
                                for (size_t p=0; p<complexSpecies[o].size(); p++)
                                {
                                    // Make sure at least one of the set is from the previous round.
                                    if (i != round-1 && k != round-1 && m != round-1 && o != round-1) continue;

                                    // Get the substrate complexes.
                                    ComplexInstance* complex1 = new ComplexInstance(*complexSpecies[i][j]);
                                    ComplexInstance* complex2 = new ComplexInstance(*complexSpecies[k][l]);
                                    ComplexInstance* complex3 = new ComplexInstance(*complexSpecies[m][n]);
                                    ComplexInstance* complex4 = new ComplexInstance(*complexSpecies[o][p]);

                                    // Get the substrate.
                                    ReactantInstance* substrate = new ReactantInstance(complex1, complex2, complex3, complex4);

                                    // Find all the matches to the substrate pattern.
                                    list<GraphMapping> matches1 = complex1->findAllIsomorphicSubgraphs(substratePattern->getReactant(0));
                                    list<GraphMapping> matches2 = complex2->findAllIsomorphicSubgraphs(substratePattern->getReactant(1));
                                    list<GraphMapping> matches3 = complex3->findAllIsomorphicSubgraphs(substratePattern->getReactant(2));
                                    list<GraphMapping> matches4 = complex4->findAllIsomorphicSubgraphs(substratePattern->getReactant(3));

                                    // Go through each match.
                                    for (auto it1=matches1.begin(); it1 != matches1.end(); it1++)
                                    {
                                        for (auto it2=matches2.begin(); it2 != matches2.end(); it2++)
                                        {
                                            for (auto it3=matches3.begin(); it3 != matches3.end(); it3++)
                                            {
                                                for (auto it4=matches4.begin(); it4 != matches4.end(); it4++)
                                                {
                                                    GraphMapping substrateMapping(substrate, substratePattern);
                                                    substrateMapping.appendMappings(*it1);
                                                    substrateMapping.appendMappings(*it2);
                                                    substrateMapping.appendMappings(*it3);
                                                    substrateMapping.appendMappings(*it4);
                                                    if (reallyVerbose) Print::printf(Print::INFO, "Found matching substrate %s", substrate->getString().c_str());

                                                    // Create any transformed products.
                                                    ReactantInstance* product = rewriteSubstrateToProduct(substrate, substrateMapping, substratePattern, reactionPattern->getSubstrateToProductMapping());
                                                    if (reallyVerbose) Print::printf(Print::INFO, "Substrate %s transformed into product %s", substrate->getString().c_str(), product->getString().c_str());

                                                    // Create any new the products, they must be completely defined.
                                                    ReactantPattern* productPattern = reactionPattern->getProducts();
                                                    for (int a=0; a<productPattern->getNumberReactants(); a++)
                                                    {
                                                        // Skip any null products.
                                                        if (productPattern->getReactant(a)->isNull())
                                                            continue;

                                                        // Skip any products that were transformed in the reaction.
                                                        if (reactionPattern->isProductReactantTransformed(a))
                                                            continue;

                                                        // Find a product that matches the pattern.
                                                        int productMatchesFound=0;
                                                        ComplexInstance* createdProductComplex=NULL;
                                                        for (int b=0; b<=round-1; b++)
                                                        {
                                                            for (size_t c=0; c<complexSpecies[b].size(); c++)
                                                            {
                                                                ComplexInstance* complexToCheck = complexSpecies[b][c];
                                                                if (complexToCheck->isIsomorphic(productPattern->getReactant(a)))
                                                                {
                                                                    if (reallyVerbose) Print::printf(Print::INFO, "Complex %s matched product pattern %s", complexToCheck->getString().c_str(), productPattern->getReactant(a)->getString().c_str());
                                                                    productMatchesFound++;
                                                                    if (createdProductComplex == NULL)
                                                                        createdProductComplex = new ComplexInstance(*complexToCheck);
                                                                }
                                                            }
                                                        }

                                                        // Make sure we found only a single matching instance.
                                                        if (productMatchesFound == 0)
                                                        {
                                                            Print::printf(Print::ERROR, "Could not find a matching complex instance for product pattern %s in reaction %s", productPattern->getReactant(a)->getString().c_str(), reactionPattern->getString(false).c_str());
                                                            throw Exception("invalid input model");
                                                        }
                                                        else if (productMatchesFound > 1)
                                                        {
                                                            Print::printf(Print::ERROR, "Found mutliple (%d) matching complex instances for product pattern %s in reaction %s", productMatchesFound, productPattern->getReactant(a)->getString().c_str(), reactionPattern->getString(false).c_str());
                                                            throw Exception("invalid input model");
                                                        }
                                                        else
                                                        {
                                                            // See if the product was not used in the substrate to product mapping.
                                                            if (reallyVerbose) Print::printf(Print::INFO, "Product %s created", createdProductComplex->getString().c_str());
                                                            product->addComplex(createdProductComplex);
                                                        }
                                                    }

                                                    // Save the new reaction.
                                                    ReactionInstance* reaction = new ReactionInstance(substrate, product, reactionPattern->getRate());
                                                    if (isNewReaction(reaction))
                                                    {
                                                        reactions[round].push_back(reaction);
                                                        if (reallyVerbose) Print::printf(Print::INFO, "Added new reaction %s", reaction->getString().c_str());
                                                    }

                                                    // If the product contains any new species, add them to the list.
                                                    for (int a=0; a<product->getNumberComplexes(); a++)
                                                    {
                                                        if (isNewComplexSpecies(product->getComplex(a)))
                                                            complexSpecies[round].push_back(product->getComplex(a));
                                                    }
                                                }
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

ReactantInstance* BNGLImporter::rewriteSubstrateToProduct(ReactantInstance* substrate, GraphMapping substrateToSubstratePatternMapping, ReactantPattern* substratePattern, GraphMapping substratePatternToProductPatternMapping)
{
    // Create a copy of the substrate to rewrite into the products.
    ReactantInstance* product = new ReactantInstance(*substrate);
    product->clearMark("transformed");

    // Create a copy of the mapping.
    GraphMapping substratePatternToProductMapping(substratePattern, product);
    for (int i=0; i<substrate->getNumberVertices(); i++)
    {
        if (substrateToSubstratePatternMapping.containsSourceVertex(substrate->getVertex(i)))
        {
            Vertex* v1 = substrateToSubstratePatternMapping.getTargetVertex(substrate->getVertex(i));
            Vertex* v2 = product->getVertex(i);
            substratePatternToProductMapping.addMapping(v1,v2);
        }
    }

    // Go through each vertext in the substrate pattern.
    set<pair<Vertex*,Vertex*>> removedEdges;
    set<pair<Vertex*,Vertex*>> addedEdges;
    for (int i=0; i<substratePattern->getNumberVertices(); i++)
    {
        Vertex* substratePatternVertex = substratePattern->getVertex(i);

        // Make sure that this vertext is involved in a substrate->product transformation.
        if (substratePatternToProductPatternMapping.isSourceVertexUsed(substratePatternVertex))
        {
            // Get the matching vertex in the product pattern.
            Vertex* productPatternVertex = substratePatternToProductPatternMapping.getTargetVertex(substratePatternVertex);

            // Get the matching vertex in the product and mark that it is part of the transformation.
            Vertex* productVertex = substratePatternToProductMapping.getTargetVertex(substratePatternVertex);
            productVertex->setMark("transformed");

            // Go through each edge and see if it was changed.
            for (int j=0; j<substratePatternVertex->getMaxNumberEdges(); j++)
            {
                // See if an edge needs to be added.
                if (substratePatternVertex->getEdge(j) == NULL && productPatternVertex->getEdge(j) != NULL)
                {
                    Vertex* v1 = productVertex;
                    int index1 = j;
                    Vertex* v2 = substratePatternToProductMapping.getTargetVertex(substratePatternToProductPatternMapping.getSourceVertex(productPatternVertex->getEdge(j)));
                    if (v2 == NULL) throw Exception("could not add the edge to the product, no second vertex", productPatternVertex->getString().c_str(), productPatternVertex->getEdge(j)->getString().c_str());
                    int index2 = productPatternVertex->getEdge(j)->findEdgeLeadingTo(productPatternVertex);


                    // Make sure we haven't yet removed this edge.
                    if (addedEdges.count(pair<Vertex*,Vertex*>(v1,v2)) == 0 && addedEdges.count(pair<Vertex*,Vertex*>(v2,v1)) == 0)
                    {
                        if (reallyVerbose) Print::printf(Print::INFO, "Adding edge between %s and %s in the product complex", v1->getString().c_str(), v2->getString().c_str());
                        if (!product->addEdge(v1, index1, v2, index2)) throw Exception("could not add the edge to the product", (v1->getString()+":"+std::to_string(index1)).c_str(), (v2->getString()+":"+std::to_string(index2)).c_str(), product->getString().c_str());
                        addedEdges.insert(pair<Vertex*,Vertex*>(v1,v2));
                    }
                }

                // See if an edge needs to be removed.
                else if (substratePatternVertex->getEdge(j) != NULL && productPatternVertex->getEdge(j) == NULL)
                {
                    Vertex* v1 = productVertex;
                    Vertex* v2 = substratePatternToProductMapping.getTargetVertex(substratePatternVertex->getEdge(j));
                    if (v2 == NULL) throw Exception("could not remove edge from the product, no second vertex", substratePatternVertex->getString().c_str(), substratePatternVertex->getEdge(j)->getString().c_str());

                    // Make sure we haven't yet removed this edge.
                    if (removedEdges.count(pair<Vertex*,Vertex*>(v1,v2)) == 0 && removedEdges.count(pair<Vertex*,Vertex*>(v2,v1)) == 0)
                    {
                        if (reallyVerbose) Print::printf(Print::INFO, "Removing edge between %s and %s in the product complex", v1->getString().c_str(), v2->getString().c_str());
                        if (!product->removeEdge(v1, v2)) throw Exception("could not remove the edge from the product", v1->getString().c_str(), v2->getString().c_str(), product->getString().c_str());
                        removedEdges.insert(pair<Vertex*,Vertex*>(v1,v2));
                    }
                }
            }

            // Go through each component and see if its state was changed.
            MoleculePattern* substratePatternMolecule = dynamic_cast<MoleculePattern*>(substratePatternVertex);
            MoleculePattern* productPatternMolecule = dynamic_cast<MoleculePattern*>(productPatternVertex);
            MoleculeInstance* productMolecule = dynamic_cast<MoleculeInstance*>(productVertex);
            if (substratePatternMolecule == NULL) throw std::runtime_error("could not cast substrate pattern to MoleculePattern in BNGLImporter::rewriteSubstrateToProduct");
            if (productPatternMolecule == NULL) throw std::runtime_error("could not cast product pattern to MoleculePattern in BNGLImporter::rewriteSubstrateToProduct");
            if (productMolecule == NULL) throw std::runtime_error("could not cast product to MoleculePattern in BNGLImporter::rewriteSubstrateToProduct");
            for (int j=0; j<substratePatternMolecule->getNumberComponents(); j++)
            {
                // See if the component state is specified in the substrate pattern.
                if (substratePatternMolecule->getComponent(j) != NULL && substratePatternMolecule->getComponent(j)->getState() != "")
                {
                    // See if the state changed between the substrate and the product.
                    if (productPatternMolecule->getComponent(j)->getState() != substratePatternMolecule->getComponent(j)->getState())
                    {
                        if (reallyVerbose) Print::printf(Print::INFO, "Changing state of component %s in the product complex %s to %s", productMolecule->getComponent(j)->getString().c_str(), productMolecule->getString().c_str(), productPatternMolecule->getComponent(j)->getState().c_str());
                        productMolecule->getComponent(j)->setState(productPatternMolecule->getComponent(j)->getState());
                    }
                }
            }
        }
    }

    // Create a new version of the product containing only complexes involved in a transformation.
    ReactantInstance* trimmedProduct = new ReactantInstance();
    trimmedProduct->setValid(true);

    // Get a list of anchor vertices in the product.
    list<Vertex*> anchorVertices = product->findConnectedSubgraphs();

    // Go through each anchor vertex.
    for (auto it=anchorVertices.begin(); it != anchorVertices.end(); it++)
    {
        MoleculeInstance* anchorMolecule = dynamic_cast<MoleculeInstance*>(*it);
        if (anchorMolecule == NULL) throw std::runtime_error("could not cast to MoleculeInstance in ReactantInstance::recreateComplexes");

        // Create a new complex from the connected molecules.
        ComplexInstance* productComplex = new ComplexInstance(anchorMolecule->getConnectedVertices());
        //Print::printf(Print::INFO, "Product complex %s %d", productComplex->getString().c_str(), productComplex->hasMarkOnAny("transformed"));

        // If any vertex in the connected graph was involved in the transformation, add it to the trimmed product.
        if (productComplex->hasMarkOnAny("transformed"))
            trimmedProduct->addComplex(productComplex);
        else
            delete productComplex;
    }

    return trimmedProduct;
}

bool BNGLImporter::isNewComplexSpecies(ComplexInstance* instance)
{
    for (size_t i=0; i<complexSpecies.size(); i++)
    {
        for (size_t j=0; j<complexSpecies[i].size(); j++)
        {
            if (complexSpecies[i][j]->isIsomorphic(instance))
                return false;
        }
    }
    return true;
}

bool BNGLImporter::isNewReaction(ReactionInstance* instance)
{
    for (size_t i=0; i<reactions.size(); i++)
    {
        for (size_t j=0; j<reactions[i].size(); j++)
        {
            if (reactions[i][j]->matches(instance) && !includeDuplicates)
            {
                Print::printf(Print::WARNING, "Reaction %s was excluded as it is a duplicate of one found in round %d step %d: %s. If this is an error, you may need to import again using the --include-duplicates flag.", instance->getString().c_str(), i, j, reactions[i][j]->getString().c_str());
                return false;
            }
        }
    }
    return true;
}

void BNGLImporter::buildSpeciesModel()
{
    // Process the species.
    numberSpecies = allComplexSpecies.size();

    // Initialize the initial species counts matrix.
    C = new ndarray<uint>(utuple(numberSpecies));
    *C=0;

    if (reallyVerbose) Print::printf(Print::INFO, "Processing %d species.", numberSpecies);
    for (uint i=0; i<numberSpecies; i++)
    {
        ComplexInstance* species = allComplexSpecies[i];
        if (!species->isValid()) throw Exception("invalid species", species->getString().c_str(), int(i));

        // Save the species name.
        speciesNames[i] = species->getString();
        (*C)[i] = lround(species->getCount());

        //Print::printf(Print::INFO, "Added species (%d) %s with initial count: %d%s%s", i, species->getId().c_str(), initialSpeciesCount, isSpeciesConst[i]?" (constant)":"", isSpeciesBoundary[i]?" (boundary)":"");
    }
}

void BNGLImporter::buildReactionModel()
{
    // Process the reactions.
    numberReactions = allReactions.size();

    // Initialize the stoichiometry matrix.
    S = new ndarray<int>(utuple(numberSpecies, numberReactions));
    *S=0;

    // Initialize the reaction type matrix.
    T = new ndarray<uint>((utuple(numberReactions)));
    *T=9999;

    // Initialize the rate constant matrix.
    K = new ndarray<double>(utuple(numberReactions,10));
    *K=NAN;

    // Initialize the dependency matrix.
    D = new ndarray<uint>(utuple(numberSpecies, numberReactions));
    *D=0;

    if (reallyVerbose) Print::printf(Print::INFO, "Processing %d reactions.", numberReactions);
    for (uint i=0; i<numberReactions; i++)
    {
        ReactionInstance* reaction = allReactions[i];
        if (!reaction->isValid()) throw Exception("invalid reaction", reaction->getString().c_str(), i);

        // Match the reaction to a propensity function.
        if (!matchPropensityFunction(i, reaction))
            THROW_EXCEPTION(lm::RuntimeException, "could not match reaction %d,%s to a propensity function", i, reaction->getString().c_str());

        // Save the reaction name.
        reactionNames[i] = reaction->getString();

        // Go through each substrate.
        for (int j=0; j<reaction->getSubstrate()->getNumberComplexes(); j++)
        {
            ComplexInstance* s = reaction->getSubstrate()->getComplex(j);

            // Update the S matrix.
            (*S)[utuple(findComplexSpeciesIndex(s),i)] -= 1;

        }

        // Go through each product.
        for (int j=0; j<reaction->getProduct()->getNumberComplexes(); j++)
        {
            ComplexInstance* s = reaction->getProduct()->getComplex(j);

            // Update the S matrix.
            (*S)[utuple(findComplexSpeciesIndex(s),i)] += 1;
        }
    }
}

int BNGLImporter::findComplexSpeciesIndex(ComplexInstance* s)
{
    int index=-1;
    for (uint i=0; i<allComplexSpecies.size(); i++)
    {
        if (s->isIsomorphic(allComplexSpecies[i]))
        {
            if (index == -1)
                index = int(i);
            else
                throw Exception("found a duplicate complex species", index, int(i));
        }
    }

    if (index == -1) throw Exception("could not find index for species");
    return index;
}

string BNGLImporter::buildReactionRateString(ReactionInstance* reaction)
{
    std::stringstream ss;
    ss << reaction->getRate();
    for (int i=0; i<reaction->getSubstrate()->getNumberComplexes(); i++)
        ss << " * " << "S" << findComplexSpeciesIndex(reaction->getSubstrate()->getComplex(i));
    return ss.str();
}

bool BNGLImporter::matchPropensityFunction(int reactionIndex, ReactionInstance* reaction)
{
    if (verbose) Print::printf(Print::INFO, "Matching reaction %d [%s] to a propensity function", reactionIndex, reaction->getString().c_str());

    // Get the kinetic expression.
    ASTNode_t* originalFormula = SBML_parseL3Formula(buildReactionRateString(reaction).c_str());
    if (verbose) Print::printf(Print::INFO, "Rate formula is: %s", SBML_formulaToL3String(originalFormula));
    if (reallyVerbose) ASTHelper::printASTNode(originalFormula);

    // Put the formula into normal form.
    ASTNode_t* normalizedFormula = originalFormula->deepCopy();
    ASTHelper::normalizeASTExpression(normalizedFormula);
    if (verbose) Print::printf(Print::INFO, "Normalized formula is: %s", SBML_formulaToL3String(normalizedFormula));
    if (reallyVerbose) ASTHelper::printASTNode(normalizedFormula);

    // Simplify the formula.
    ASTNode_t* simplifiedFormula = normalizedFormula->deepCopy();
    ASTHelper::simplifyASTExpression(simplifiedFormula);
    ASTHelper::sortASTExpression(simplifiedFormula);
    if (verbose) Print::printf(Print::INFO, "Simplified formula is: %s", SBML_formulaToL3String(simplifiedFormula));
    if (reallyVerbose) ASTHelper::printASTNode(simplifiedFormula);

    // Iterate through each propensity function and see if it matches.
    map<uint,lm::me::PropensityFunctionDefinition> functions = propensityFunctions->getFunctions();
    for (std::map<uint,lm::me::PropensityFunctionDefinition>::const_iterator it=functions.begin(); it != functions.end(); it++)
    {
        uint id = it->first;
        lm::me::PropensityFunctionDefinition p = it->second;
        if (p.expressions.size() > 0)
        {
            for (list<string>::iterator it = p.expressions.begin(); it != p.expressions.end(); it++)
            {
                ASTNode_t* propensityFormula = SBML_parseL3Formula(it->c_str());
                ASTNode_t* normalizedPropensityFormula = propensityFormula->deepCopy();
                ASTHelper::normalizeASTExpression(normalizedPropensityFormula);

                if (ASTHelper::compareASTNodes(simplifiedFormula, normalizedPropensityFormula))
                {
                    // Create the entry for this formula.
                    if (verbose) Print::printf(Print::INFO, "Matched reaction %d to propensity function %s: [%s] == [%s]", reactionIndex, p.name.c_str(), SBML_formulaToL3String(simplifiedFormula), SBML_formulaToL3String(normalizedPropensityFormula));
                    (*T)[utuple(static_cast<uint>(reactionIndex))] = id;
                    return createPropensityFunctionEntry(reactionIndex, simplifiedFormula, normalizedPropensityFormula, p);
                }
                else
                {
                    if (reallyVerbose) Print::printf(Print::INFO, "No match to propensity function %s: [%s]", p.name.c_str(), SBML_formulaToL3String(normalizedPropensityFormula));
                    if (reallyVerbose) ASTHelper::printASTNode(normalizedPropensityFormula);
                }
            }
        }
    }

    // Print out some messages to help the user figure out why there wasn't a match.
    Print::printf(Print::ERROR, "FAILED to match kinetic formula in reaction %d [%s] to a propensity function: [%s] ", reactionIndex, reaction->getString().c_str(), SBML_formulaToL3String(simplifiedFormula));
    if (verbose)
    {
        Print::printf(Print::ERROR, "                                         Normalized formula: [%s]", SBML_formulaToL3String(normalizedFormula));
        Print::printf(Print::ERROR, "                                         Original formula:   [%s]", SBML_formulaToL3String(originalFormula));
        Print::printf(Print::ERROR, "Abstract syntax tree for simplified form was:");
        ASTHelper::printASTNode(simplifiedFormula);
    }

    return false;
}

bool BNGLImporter::createPropensityFunctionEntry(int reactionIndex, ASTNode_t* formula, ASTNode_t* propensityFormula, lm::me::PropensityFunctionDefinition& propensityFunction)
{
    // If this is a number and it matches to a k, store the parameter.
    if (formula->isNumber() && propensityFormula->isName() && propensityFormula->getName()[0] == 'k')
    {
        int parameterIndex = atoi(propensityFormula->getName()+1)-1;
        utuple index = utuple(static_cast<uint>(reactionIndex),static_cast<uint>(parameterIndex));
        if (formula->getType() == AST_INTEGER)
        {
            double value = static_cast<double>(formula->getInteger());
            if (std::isnan((*K)[index]))
            {
                (*K)[index] = value;
                if (verbose) Print::printf(Print::INFO, "    Added parameter for reaction %d parameter %d: %e", reactionIndex, parameterIndex, (*K)[index]);
                return true;
            }
            else if (almostEqual((*K)[index],value))
            {
                return true;
            }
            else
            {
                Print::printf(Print::ERROR, "FAILED to create entry for reaction %d parameter %d: the specified value did not match the previous value for this constant, %e != %e.", reactionIndex, parameterIndex, value, (*K)[index]);
                return false;
            }
        }
        else if (formula->getType() == AST_REAL || formula->getType() == AST_REAL_E)
        {
            double value = formula->getReal();
            if (std::isnan((*K)[index]))
            {
                (*K)[index] = value;
                if (verbose) Print::printf(Print::INFO, "    Added parameter for reaction %d parameter %d: %e", reactionIndex, parameterIndex, (*K)[index]);
                return true;
            }
            else if (almostEqual((*K)[index],value))
            {
                return true;
            }
            else
            {
                Print::printf(Print::ERROR, "FAILED to create entry for reaction %d parameter %d: the specified value did not match the previous value for this constant, %e != %e.", reactionIndex, parameterIndex, value, (*K)[index]);
                return false;
            }
        }
        Print::printf(Print::ERROR, "FAILED to create entry for reaction %d parameter %d: the value did not match a known numeric type.", reactionIndex, parameterIndex);
        return false;
    }

    // If this is species index and it matches to an x, store the dependency.
    if (formula->isName() && formula->getName()[0] == 'S' && propensityFormula->isName() && propensityFormula->getName()[0] == 'x')
    {
        // Get the species index and the variable number.
        int speciesIndex = atoi(formula->getName()+1);
        int dependencyNumber = atoi(propensityFormula->getName()+1);

        // Add an entry in the D matrix.
        utuple index = utuple(static_cast<uint>(speciesIndex),static_cast<uint>(reactionIndex));
        if ((*D)[index] == 0)
        {
            (*D)[index] = static_cast<uint>(dependencyNumber);
            if (verbose) Print::printf(Print::INFO, "    Added dependency for reaction %d on species %d (%s): %d", reactionIndex, speciesIndex, formula->getName(), dependencyNumber);
            return true;
        }
        else if ((*D)[index] == static_cast<uint>(dependencyNumber))
        {
            return true;
        }
        else
        {
            Print::printf(Print::ERROR, "FAILED to create entry for reaction %d species %d: the specified dependency did not match the previous dependency for this species, %d != %d.", reactionIndex, speciesIndex, dependencyNumber, (*D)[index]);
            return false;
        }
    }

    // Go through all of the children, if any failed return false.
    bool ret = true;
    for (uint i=0; i<formula->getNumChildren(); i++)
        if (!createPropensityFunctionEntry(reactionIndex,formula->getChild(i), propensityFormula->getChild(i), propensityFunction))
            ret = false;
    return ret;
}

void BNGLImporter::buildLMModel()
{
    lmModel.set_number_species(numberSpecies);
    for (uint i=0; i<numberSpecies; i++)
    {
        // Add the species to the model.
        lmModel.add_species_name(speciesNames[i]);
        lmModel.add_initial_species_count((*C)[i]);
    }
    lmModel.set_number_reactions(numberReactions);

    // Fill in the S and D matrices.
    for (uint i=0; i<numberSpecies; i++)
    {
        for (uint j=0; j<numberReactions; j++)
        {
            lmModel.add_stoichiometric_matrix((*S)[utuple(i,j)]);
            lmModel.add_dependency_matrix((*D)[utuple(i,j)]);
        }
    }

    // Fill in the T and K matrices.
    for (uint j=0; j<numberReactions; j++)
    {
        lm::input::ReactionModel_Reaction* reaction = lmModel.add_reaction();
        reaction->set_name(reactionNames[j]);
        reaction->set_type((*T)[utuple(j)]);
        for (uint k=0; k<10 && !std::isnan((*K)[utuple(j,k)]); k++)
            reaction->add_rate_constant((*K)[utuple(j,k)]);
    }
}

}
}
