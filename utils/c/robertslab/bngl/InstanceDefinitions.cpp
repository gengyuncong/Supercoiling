/*
 * Copyright 2017 Johns Hopkins University
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
#include <regex>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "robertslab/bngl/InstanceDefinitions.h"
#include "robertslab/bngl/PatternDefinitions.h"

using std::regex;
using std::regex_token_iterator;
using std::set;
using std::string;
using std::vector;

namespace robertslab {
namespace bngl {

ComponentInstance::ComponentInstance()
:molecule(NULL),valid(false),bond(NULL)
{
}

ComponentInstance::ComponentInstance(MoleculeInstance* molecule, string definition)
:molecule(molecule),valid(false),bond(NULL)
{
    // Parse the component definition.
    regex statePattern("^([^~!]+)(?:~([^~!]+))?(?:!([^~!]+))?$");
    std::smatch match;
    if (std::regex_match(definition, match, statePattern) && match.size() == 4)
    {
        name = match[1].str();
        state = match[2].str();
        bondName = match[3].str();
        valid = true;
    }
}

ComponentInstance::ComponentInstance(MoleculeInstance* molecule, const ComponentInstance& other)
:molecule(molecule),valid(other.valid),name(other.name),state(other.state),bondName(other.bondName),bond(NULL)
{
}

bool ComponentInstance::isValid()
{
    return valid;
}

void ComponentInstance::setState(string newState)
{
    state = newState;
}

string ComponentInstance::getString(bool includeBondNames)
{
    std::stringstream ss;
    ss << name;
    if (state != "") ss << "~" << state;
    if (includeBondNames && bondName != "") ss << "!" << bondName;
    return ss.str();
}

MoleculeInstance::MoleculeInstance()
:valid(false)
{
}

MoleculeInstance::MoleculeInstance(string definition)
:valid(false)
{
    regex moleculePattern("^(\\w+)\\(([^\\)]*)\\)$");
    std::smatch match;
    if (std::regex_match(definition, match, moleculePattern) && match.size() == 3)
    {
        valid = true;
        name = match[1].str();

        // Parse the components.
        string componentString = match[2].str();
        regex componentPattern("([^ ,]+)");
        regex_token_iterator<string::iterator> endOfTokens;
        regex_token_iterator<std::string::iterator> tokens(componentString.begin(), componentString.end(), componentPattern);
        while (tokens != endOfTokens)
        {
            string componentString = *tokens++;
            ComponentInstance* component = new ComponentInstance(this, componentString);
            components.push_back(component);
            if (!component->isValid()) valid = false;
        }
    }
}

MoleculeInstance::MoleculeInstance(const MoleculeInstance& other)
:valid(other.valid),name(other.name)
{
    for (size_t i=0; i<other.components.size(); i++)
        components.push_back(new ComponentInstance(this, *other.components[i]));
}

bool MoleculeInstance::isValid()
{
    return valid;
}

string MoleculeInstance::getName()
{
    return name;
}

ComponentInstance* MoleculeInstance::getComponent(int index)
{
    return components[index];
}

string MoleculeInstance::getString()
{
    return getString(true);
}

string MoleculeInstance::getString(bool includeBondNames)
{
    std::stringstream ss;
    ss << name << "(";
    for (size_t i=0; i<components.size(); i++)
    {
        ss << (i==0?"":",") << components[i]->getString(includeBondNames);
    }
    ss << ")";

    return ss.str();
}

bool MoleculeInstance::matches(MoleculeInstance* instance)
{
    if (name != instance->name) return false;
    if (components.size() != instance->components.size()) return false;

    // Check the bonding state of the components.
    for (size_t i=0; i<components.size(); i++)
    {
        // Make sure the components have the same name.
        if (components[i]->name != instance->components[i]->name) return false;

        // See if the source has a bond.
        if (components[i]->bond != NULL)
        {
            // The source has a bond, so make sure the target does to.
            if (instance->components[i]->bond == NULL) return false;

            // TODO: make sure the bond is the same.
        }
        else
        {
            // Otherwise the source didn't have a bond, so make sure the target doesn't either.
            if (instance->components[i]->bond != NULL) return false;
        }

        // Make sure the component state is the same.
        if (components[i]->state != instance->components[i]->state) return false;
    }

    return true;
}

int MoleculeInstance::getMaxNumberEdges()
{
    return components.size();
}

Vertex* MoleculeInstance::getEdge(int i)
{
    if (components[i]->bond != NULL)
        return components[i]->bond->molecule;
    return NULL;
}

void MoleculeInstance::addEdge(int sourceIndex, Vertex* dest, int destIndex)
{
    MoleculeInstance* destMolecule = dynamic_cast<MoleculeInstance*>(dest);
    if (destMolecule == NULL) throw std::runtime_error("could not cast to MoleculeInstance in MoleculeInstance::addEdge");
    components[sourceIndex]->bondName = "*";
    components[sourceIndex]->bond = destMolecule->components[destIndex];
}

void MoleculeInstance::removeEdge(int index)
{
    components[index]->bondName = "";
    components[index]->bond = NULL;
}

bool MoleculeInstance::matches(Vertex* comp)
{
    if (dynamic_cast<MoleculeInstance*>(comp))
    {
        return matches((MoleculeInstance*)comp);
    }
    else if (dynamic_cast<MoleculePattern*>(comp))
    {
        return ((MoleculePattern*)comp)->matches(this);
    }

    return false;
}


ComplexInstance::ComplexInstance()
:valid(false),count(0.0)
{
}

ComplexInstance::ComplexInstance(string definition, double count)
:valid(false),count(count)
{
    valid = true;

    // Parse the molecule definitions.
    regex moleculePattern("([^\\.]+)");
    regex_token_iterator<string::iterator> endOfTokens;
    regex_token_iterator<std::string::iterator> tokens(definition.begin(), definition.end(), moleculePattern);
    while (tokens != endOfTokens)
    {
        string moleculeString = *tokens++;
        MoleculeInstance* molecule = new MoleculeInstance(moleculeString);
        molecules.push_back(molecule);
        if (!molecule->isValid()) valid = false;
    }

    // Connect all of the bonds.
    connectBonds();
}

ComplexInstance::ComplexInstance(set<Vertex*> connectedMolecules)
:valid(false),count(0.0)
{
    valid = true;

    // Add all of the moelcules.
    for (auto it=connectedMolecules.begin(); it != connectedMolecules.end(); it++)
    {
        MoleculeInstance* molecule = dynamic_cast<MoleculeInstance*>(*it);
        if (molecule == NULL) throw std::runtime_error("could not cast to MoleculeInstance in ComplexInstance::ComplexInstance");
        molecules.push_back(molecule);
        if (!molecule->isValid()) valid = false;
    }

    // Construct the bond names.
    constructBondNames();
}


ComplexInstance::ComplexInstance(const ComplexInstance& other)
:valid(other.valid),count(0.0)
{
    // Create the new molecules.
    for (size_t i=0; i<other.molecules.size(); i++)
        molecules.push_back(new MoleculeInstance(*other.molecules[i]));

    // Connect all of the bonds.
    connectBonds();
}

void ComplexInstance::connectBonds()
{
    // Go through and establish the connectivity using the bond names.
    for (size_t m1=0; m1<molecules.size(); m1++)
    {
        for (size_t c1=0; c1<molecules[m1]->components.size(); c1++)
        {
            if (molecules[m1]->components[c1]->bondName != "")
            {
                int matches=0;
                for (size_t m2=0; m2<molecules.size(); m2++)
                {
                    for (size_t c2=0; c2<molecules[m2]->components.size(); c2++)
                    {
                        if ((m1 != m2 || c1 != c2) && molecules[m1]->components[c1]->bondName == molecules[m2]->components[c2]->bondName)
                        {
                            molecules[m1]->components[c1]->bond = molecules[m2]->components[c2];
                            matches++;
                        }
                    }
                }

                if (matches != 1) throw std::invalid_argument("inconsistent number of bond names in ComplexInstance::ComplexInstance");
            }
            else
            {
                molecules[m1]->components[c1]->bond = NULL;
            }
        }
    }
}

void ComplexInstance::constructBondNames()
{
    // Go through and reset all of the bond names.
    for (size_t i=0; i<molecules.size(); i++)
    {
        for (size_t j=0; j<molecules[i]->components.size(); j++)
        {
            if (molecules[i]->components[j] != NULL)
            {
                molecules[i]->components[j]->bondName = "";
            }
        }
    }

    // Go through and specify all of the bond names.
    int nextBond=1;
    for (size_t i=0; i<molecules.size(); i++)
    {
        for (size_t j=0; j<molecules[i]->components.size(); j++)
        {
            if (molecules[i]->components[j] != NULL)
            {
                if (molecules[i]->components[j]->bond != NULL && molecules[i]->components[j]->bondName == "")
                {
                    // Set the name of this bond and the target.
                    molecules[i]->components[j]->bondName = std::to_string(nextBond);
                    molecules[i]->components[j]->bond->bondName = std::to_string(nextBond);
                    nextBond++;
                }
            }
        }
    }
}

bool ComplexInstance::isValid() const
{
    return valid;
}

double ComplexInstance::getCount()
{
    return count;
}


string ComplexInstance::getString()
{
    return getString(false,true);
}

string ComplexInstance::getString(bool includeCounts, bool includeBondNames)
{
    // Get the alphbetical order for the molecules.
    vector<MoleculeInstance*> sortedMolecules;
    for (auto it=molecules.begin(); it != molecules.end(); it++)
    {
        bool added=false;
        for (auto it2=sortedMolecules.begin(); it2 != sortedMolecules.end(); it2++)
        {
            if ((*it)->getString(false).compare((*it2)->getString(false)) <= 0)
            {
                sortedMolecules.insert(it2, *it);
                added = true;
                break;
            }
        }
        if (!added) sortedMolecules.push_back(*it);
    }

    // Construct the string.
    std::stringstream ss;
    for (size_t i=0; i<sortedMolecules.size(); i++)
    {
        ss << (i==0?"":".") << sortedMolecules[i]->getString(includeBondNames);
    }
    if (includeCounts) ss << " " << count;
    return ss.str();
}


int ComplexInstance::getNumberMolecules()
{
    return molecules.size();
}

MoleculeInstance* ComplexInstance::getMolecule(int i)
{
    return molecules[i];
}

/*bool ComplexInstance::matches(ComplexInstance* instance)
{
    printf("Comparing %s to %s\n",getString().c_str(), instance->getString().c_str());
    if (molecules.size() != instance->molecules.size()) return false;

    // Make sure we can find every molecule, even though they may not be in the same order.
    for (int i=0; i<molecules.size(); i++)
    {
        bool foundMatch=false;
        for (int j=0; j<instance->molecules.size(); j++)
        {
            if (molecules[i]->matches(instance->molecules[j]))
            {
                foundMatch = true;
                break;
            }
        }
        if (!foundMatch) return false;
    }
    printf("yes\n");
    return true;
}*/


int ComplexInstance::getNumberVertices()
{
    return molecules.size();
}

Vertex* ComplexInstance::getVertex(int i)
{
    return molecules[i];
}

ReactantInstance::ReactantInstance()
:valid(false)
{
}

ReactantInstance::ReactantInstance(ComplexInstance* complex1)
{
    valid = true;
    if (!complex1->isValid()) valid = false;
    complexes.push_back(complex1);
}

ReactantInstance::ReactantInstance(ComplexInstance* complex1, ComplexInstance* complex2)
{
    valid = true;
    if (!complex1->isValid()) valid = false;
    if (!complex2->isValid()) valid = false;
    complexes.push_back(complex1);
    complexes.push_back(complex2);
}

ReactantInstance::ReactantInstance(ComplexInstance* complex1, ComplexInstance* complex2, ComplexInstance* complex3)
{
    valid = true;
    if (!complex1->isValid()) valid = false;
    if (!complex2->isValid()) valid = false;
    if (!complex3->isValid()) valid = false;
    complexes.push_back(complex1);
    complexes.push_back(complex2);
    complexes.push_back(complex3);
}

ReactantInstance::ReactantInstance(ComplexInstance* complex1, ComplexInstance* complex2, ComplexInstance* complex3, ComplexInstance* complex4)
{
    valid = true;
    if (!complex1->isValid()) valid = false;
    if (!complex2->isValid()) valid = false;
    if (!complex3->isValid()) valid = false;
    if (!complex4->isValid()) valid = false;
    complexes.push_back(complex1);
    complexes.push_back(complex2);
    complexes.push_back(complex3);
    complexes.push_back(complex4);
}

ReactantInstance::ReactantInstance(const ReactantInstance& other)
:valid(other.valid)
{
    // Create the new complexes.
    for (size_t i=0; i<other.complexes.size(); i++)
        complexes.push_back(new ComplexInstance(*other.complexes[i]));

}

bool ReactantInstance::isValid()
{
    return valid;
}

void ReactantInstance::setValid(bool valid)
{
    this->valid = valid;
}

string ReactantInstance::getString()
{
    return getString(true);
}

string ReactantInstance::getString(bool includeBondNames=true)
{
    // Get the alphbetical order for the complexes.
    vector<ComplexInstance*> sortedComplexes;
    for (auto it=complexes.begin(); it != complexes.end(); it++)
    {
        bool added=false;
        for (auto it2=sortedComplexes.begin(); it2 != sortedComplexes.end(); it2++)
        {
            if ((*it)->getString(false,false).compare((*it2)->getString(false,false)) <= 0)
            {
                sortedComplexes.insert(it2, *it);
                added = true;
                break;
            }
        }
        if (!added) sortedComplexes.push_back(*it);
    }

    // Construct the string.
    std::stringstream ss;
    for (size_t i=0; i<sortedComplexes.size(); i++)
        ss << (i==0?"":" + ") << sortedComplexes[i]->getString(false,includeBondNames);
    return ss.str();
}

int ReactantInstance::getNumberComplexes()
{
    return complexes.size();
}

ComplexInstance* ReactantInstance::getComplex(int index)
{
    return complexes[index];
}

void ReactantInstance::addComplex(ComplexInstance* complex)
{
    if (!complex->isValid()) valid = false;
    complexes.push_back(complex);
}

int ReactantInstance::getNumberVertices()
{
    int ret=0;
    for (size_t i=0; i<complexes.size(); i++)
        ret += complexes[i]->getNumberVertices();
    return ret;
}

Vertex* ReactantInstance::getVertex(int index)
{
    for (size_t i=0; i<complexes.size(); i++)
    {
        if (index <  complexes[i]->getNumberVertices())
            return complexes[i]->getVertex(index);
        index -= complexes[i]->getNumberVertices();
    }
    throw std::out_of_range("index out of range in call to ReactantInstance::getVertex");
}

void ReactantInstance::recreateComplexes()
{
    // Get a list of vertices that anchor a set of connected subgraphs.
    list<Vertex*> anchorVertices = findConnectedSubgraphs();

    // Clear the complexes list.
    complexes.clear();

    // Go through each anchor vertex.
    for (auto it=anchorVertices.begin(); it != anchorVertices.end(); it++)
    {
        MoleculeInstance* anchorMolecule = dynamic_cast<MoleculeInstance*>(*it);
        if (anchorMolecule == NULL) throw std::runtime_error("could not cast to MoleculeInstance in ReactantInstance::recreateComplexes");
        complexes.push_back(new ComplexInstance(anchorMolecule->getConnectedVertices()));
    }

}

ReactionInstance::ReactionInstance()
:valid(false),substrate(NULL),product(NULL),rate(0.0)
{
}

ReactionInstance::ReactionInstance(ReactantInstance* substrate, ReactantInstance* product, double rate)
:valid(true),substrate(substrate),product(product),rate(rate)
{
}

bool ReactionInstance::isValid()
{
    return valid;
}

string ReactionInstance::getString(bool includeRate, bool includeBondNames)
{
    std::stringstream ss;
    ss << substrate->getString(includeBondNames);
    ss << " -> ";
    ss << product->getString(includeBondNames);
    if (includeRate) ss << " " << rate;
    return ss.str();
}

ReactantInstance* ReactionInstance::getSubstrate()
{
    return substrate;
}

ReactantInstance* ReactionInstance::getProduct()
{
    return product;
}

double ReactionInstance::getRate()
{
    return rate;
}

bool ReactionInstance::matches(ReactionInstance* target)
{
    // Make sure the number of substrates and products are the same.
    if (substrate->getNumberComplexes() != target->substrate->getNumberComplexes()) return false;
    if (product->getNumberComplexes() != target->product->getNumberComplexes()) return false;

    // Make sure the substratea are the same.
    for (size_t i=0; i<substrate->getNumberComplexes(); i++)
    {
        if (!substrate->getComplex(i)->isIsomorphic(target->substrate->getComplex(i)))
            return false;
    }

    // Make sure the products are the same.
    for (size_t i=0; i<product->getNumberComplexes(); i++)
    {
        if (!product->getComplex(i)->isIsomorphic(target->product->getComplex(i)))
            return false;
    }

    // Make sure the rates are the same.
    if (fabs((rate-target->rate)/rate) > 1e-6) return false;

    return true;
}

}
}
