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

#ifndef ROBERTSLAB_BNGL_INSTANCEDEFINITIONS_H
#define ROBERTSLAB_BNGL_INSTANCEDEFINITIONS_H

#include <string>
#include <vector>

#include "robertslab/graph/Graph.h"

using std::string;
using std::vector;
using robertslab::graph::Graph;
using robertslab::graph::Vertex;

namespace robertslab {
namespace bngl {

class ComponentClass;
class MoleculeClass;
class MoleculeInstance;
class MoleculePattern;
class ReactantInstance;

class ComponentInstance
{
public:
    ComponentInstance();
    ComponentInstance(MoleculeInstance* molecule, string definition);
    ComponentInstance(MoleculeInstance* molecule, const ComponentInstance& other);
    virtual ~ComponentInstance() {}
    bool isValid();
    void setState(string newState);
    string getString(bool includeBondNames=true);
    ComponentInstance* getBond();

protected:
    MoleculeInstance* molecule;
    bool valid;
    string name;
    string state;
    string bondName;
    ComponentInstance* bond;

    friend class ComplexInstance;
    friend class MoleculeInstance;
    friend class ComponentClass;
    friend class MoleculePattern;
};

class MoleculeInstance : public Vertex
{
public:
    MoleculeInstance();
    MoleculeInstance(string definition);
    MoleculeInstance(const MoleculeInstance& other);
    virtual ~MoleculeInstance() {}
    bool isValid();
    string getName();
    ComponentInstance* getComponent(int index);
    bool matches(MoleculeInstance* instance);
    string getString(bool includeBondNames);

public:
    virtual int getMaxNumberEdges();
    virtual Vertex* getEdge(int i);
    virtual void addEdge(int sourceIndex, Vertex* dest, int destIndex);
    virtual void removeEdge(int i);
    virtual bool matches(Vertex* v2);
    virtual string getString();

protected:
    bool valid;
    string name;
    vector<ComponentInstance*> components;

    friend class ComplexInstance;
    friend class MoleculeClass;
    friend class MoleculePattern;
};

class ComplexInstance : public Graph
{
public:
    ComplexInstance();
    ComplexInstance(string definition, double count);
    ComplexInstance(set<Vertex*> molecules);
    ComplexInstance(const ComplexInstance& other);
    virtual ~ComplexInstance() {}
    bool isValid() const;
    double getCount();
    int getNumberMolecules();
    MoleculeInstance* getMolecule(int i);
    string getString(bool includeCounts, bool includeBondNames);

public:
    virtual int getNumberVertices();
    virtual Vertex* getVertex(int i);
    virtual string getString();

protected:
    void connectBonds();
    void constructBondNames();

protected:
    bool valid;
    double count;
    vector<MoleculeInstance*> molecules;
};

class ReactantInstance : public Graph
{
public:
    ReactantInstance();
    ReactantInstance(ComplexInstance* complex1);
    ReactantInstance(ComplexInstance* complex1, ComplexInstance* complex2);
    ReactantInstance(ComplexInstance* complex1, ComplexInstance* complex2, ComplexInstance* complex3);
    ReactantInstance(ComplexInstance* complex1, ComplexInstance* complex2, ComplexInstance* complex3, ComplexInstance* complex4);
    ReactantInstance(const ReactantInstance& other);
    virtual ~ReactantInstance() {}
    bool isValid();
    void setValid(bool valid);
    int getNumberComplexes();
    ComplexInstance* getComplex(int index);
    void addComplex(ComplexInstance* complex);
    void recreateComplexes();
    string getString(bool includeBondNames);

public:
    virtual int getNumberVertices();
    virtual Vertex* getVertex(int i);
    virtual string getString();

protected:
    bool valid;
    vector<ComplexInstance*> complexes;
};

class ReactionInstance
{
public:
    ReactionInstance();
    ReactionInstance(ReactantInstance* substrate, ReactantInstance* product, double rate);
    virtual ~ReactionInstance() {}
    bool isValid();
    string getString(bool includeRate=false, bool includeBondNames=true);
    ReactantInstance* getSubstrate();
    ReactantInstance* getProduct();
    double getRate();
    bool matches(ReactionInstance* target);

protected:
    bool valid;
    ReactantInstance* substrate;
    ReactantInstance* product;
    double rate;
};

}
}

#endif // MOLECULE_H
