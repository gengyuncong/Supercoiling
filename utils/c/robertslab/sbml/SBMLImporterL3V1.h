#ifndef SBMLIMPORTERL3V1_H
#define SBMLIMPORTERL3V1_H

#include <map>
#include <string>
#include <vector>
#include <sbml/SBMLDocument.h>

#include "lm/Types.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/me/PropensityFunction.h"

using std::map;
using std::string;
using std::vector;

namespace robertslab {
namespace sbml {

class SBMLImporterL3V1
{
public:
    SBMLImporterL3V1();
    virtual ~SBMLImporterL3V1();
    virtual void setOptions(bool constantsUseConcentrations, bool verbose, bool reallyVerbose, bool ignoreErrors, bool ignoreUnmatchedReactions, bool ignoreVariableParameters);
    virtual bool import(SBMLDocument* document, map<string,double> userParameters, map<string,string> userExpressions);
    virtual lm::input::ReactionModel* getReactionModel();

protected:
    virtual string getDescription();
    virtual void expandFunctionDefinitions();
    virtual void importGlobalExpressions();
    virtual void importGlobalParameters();
    virtual void importCompartments();
    virtual void importSpecies();
    virtual void importReactions();
    virtual bool importKinetics(Reaction* reaction, int reactionIndex, KineticLaw* kinetics);
    virtual bool importPropensityFunction(Reaction* reaction, int reactionIndex, KineticLaw* kinetics, map<string,double>& parameterValues);
    virtual bool assignPropensityFunction(Reaction* reaction, int reactionIndex, KineticLaw* kinetics, string propensityFunctionName, map<string,double>& parameterValues);
    virtual bool matchPropensityFunction(Reaction* reaction, int reactionIndex, KineticLaw* kinetics, map<string,double>& parameterValues);
    virtual bool createPropensityFunctionEntry(int reactionIndex, ASTNode_t* formula, ASTNode_t* propensityFormula, lm::me::PropensityFunctionDefinition& propensityFunction);

protected:
    virtual ASTNode_t* filterKineticExpression(ASTNode_t* expression);

protected:
    virtual double convertPropensityConstantUnits(string constant, double value, string desiredUnits);
    virtual void convertUnits(ASTNode_t* units);
    virtual double convertVolumeToLiters(double size, string units="");
    virtual double convertSubstanceToParticles(double value, string units="");
    virtual double convertTimeToSeconds(double value, string units="");

protected:
    lm::me::PropensityFunctionFactory *propensityFunctions;
    SBMLDocument* sbmlDocument;
    Model* sbmlModel;
    bool constantsUseConcentrations;
    bool verbose, reallyVerbose;
    bool stopOnError, stopOnUnmatchedReactions;
    bool ignoreVariableParameters;
    map<string,double> userParameters;
    map<string,string> userExpressions;
    bool allImportStepsSuccessful;
    map<string,double> globalParameters;
    map<string,ASTNode_t*> globalExpressions;
    vector<string> compartments;
    map<string,double> compartmentSizes;
    lm::input::ReactionModel reactionModel;
    int numberSpecies;
    map<string,int> speciesIndices;
    map<int,bool> isSpeciesConst;
    map<int,bool> isSpeciesBoundary;
    int numberReactions;
    ndarray<int> *S;
    ndarray<int> *T;
    ndarray<double> *K;
    ndarray<int> *D;
};

}
}

#endif // SBMLIMPORTERL3V1_H
