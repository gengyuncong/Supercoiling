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

#ifndef ROBERTSLAB_BNGL_BNGLIMPORTER_H
#define ROBERTSLAB_BNGL_BNGLIMPORTER_H

#include <map>
#include <string>
#include <vector>
#include <sbml/math/ASTNode.h>

#include "lm/Types.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/me/PropensityFunction.h"
#include "robertslab/bngl/InstanceDefinitions.h"
#include "robertslab/bngl/PatternDefinitions.h"
#include "robertslab/bngl/TypeDefinitions.h"

using std::map;
using std::string;
using std::vector;

namespace robertslab {
namespace bngl {

class BNGLImporter
{
public:
    BNGLImporter();
    virtual ~BNGLImporter();
    virtual string getDescription();
    virtual void setOptions(bool constantsUseConcentrations, bool verbose, bool reallyVerbose, bool sorted, bool includeDuplicates);
    virtual bool import(string filename, map<string,double> userParameters);
    virtual lm::input::ReactionModel* getReactionModel();

protected:
    virtual void parseParameters(list<string>& lines);
    virtual void parseMoleculeTypes(list<string>& lines);
    virtual void parseSeedSpecies(list<string>& lines);
    virtual void parseReactions(list<string>& lines);
    bool evaluteExpression(string expression, double& value);

protected:
    virtual void processModel();
    virtual void supplementMoleculeTypesFromSeedSpecies();
    virtual void enumerateMoleculeSpecies();
    virtual bool processReactions(int round);
    virtual void processReaction(int round, ReactionPattern* reaction);
    virtual void processReactionZerothOrder(int round, ReactionPattern* reaction);
    virtual void processReactionFirstOrder(int round, ReactionPattern* reaction);
    virtual void processReactionSecondOrder(int round, ReactionPattern* reaction);
    virtual void processReactionThirdOrder(int round, ReactionPattern* reaction);
    virtual void processReactionFourthOrder(int round, ReactionPattern* reaction);

    ReactantInstance* rewriteSubstrateToProduct(ReactantInstance* substrate, GraphMapping substrateToSubstratePatternMapping, ReactantPattern* substratePattern, GraphMapping substratePatternToProductPatternMapping);
    bool isNewComplexSpecies(ComplexInstance* instance);
    bool isNewReaction(ReactionInstance* instance);

    void buildSpeciesModel();
    void buildReactionModel();
    void buildLMModel();
    int findComplexSpeciesIndex(ComplexInstance* s);
    string buildReactionRateString(ReactionInstance* reaction);

    bool matchPropensityFunction(int reactionIndex, ReactionInstance* reaction);
    bool createPropensityFunctionEntry(int reactionIndex, ASTNode_t* formula, ASTNode_t* propensityFormula, lm::me::PropensityFunctionDefinition& propensityFunction);


    /*virtual double convertPropensityConstantUnits(string constant, double value, string desiredUnits);
    virtual void convertUnits(ASTNode_t* units);
    virtual double convertVolumeToLiters(double size, string units="");
    virtual double convertSubstanceToParticles(double value, string units="");
    virtual double convertTimeToSeconds(double value, string units="");*/

protected:
    int maxRounds;
    lm::me::PropensityFunctionFactory *propensityFunctions;
    bool constantsUseConcentrations;
    bool verbose, reallyVerbose, sorted, includeDuplicates;
    bool allImportStepsSuccessful;
    map<string,double> parameters;
    map<string,MoleculeClass*> moleculeTypes;
    vector<MoleculeInstance*> moleculeSpecies;
    vector<ReactionPattern*> reactionPatterns;
    vector<vector<ComplexInstance*>> complexSpecies;
    vector<vector<ReactionInstance*>> reactions;
    vector<ComplexInstance*> allComplexSpecies;
    vector<ReactionInstance*> allReactions;


    lm::input::ReactionModel lmModel;

    uint numberSpecies;
    map<uint,string> speciesNames;
    ndarray<uint>* C; // Initial species counts.

    uint numberReactions;
    map<uint,string> reactionNames;
    ndarray<int>* S;
    ndarray<uint>* T;
    ndarray<double>* K;
    ndarray<uint>* D;
};

}
}

#endif // BNGLImporter_H
