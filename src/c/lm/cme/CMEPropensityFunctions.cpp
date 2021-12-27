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

#include <list>
#include <map>
#include <string>
#include <vector>
#include "math.h"

#include "lm/ClassFactory.h"
#include "lm/cme/CMEPropensityFunctions.h"
#include "lm/me/PropensityFunction.h"
#include "lm/Types.h"

using std::list;
using std::map;
using std::string;
using std::vector;

namespace lm {
namespace cme {

bool CMEPropensityFunctions::registered=CMEPropensityFunctions::registerClass();

bool CMEPropensityFunctions::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::me::PropensityFunctionCollection", "lm::cme::CMEPropensityFunctions", &CMEPropensityFunctions::allocateObject);
    return true;
}

void* CMEPropensityFunctions::allocateObject()
{
    return new CMEPropensityFunctions();
}

CMEPropensityFunctions::CMEPropensityFunctions()
{
}

CMEPropensityFunctions::~CMEPropensityFunctions()
{
}

class ZerothOrderPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 0;

    ZerothOrderPropensity(double k) :PropensityFunction(REACTION_TYPE,0),k(k) {}
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    virtual void changeVolume(double volumeMultiplier) {k*=volumeMultiplier;}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        return k;
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return _mm256_set1_pd(k);
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple dependencies = getDependencies(reactionIndex, D);
        if (dependencies.len != 0) throw InvalidArgException("D", "zeroth order propensity had invalid number of dependencies",dependencies.len);

        // Find the rate costant.
        if (k.len < 1)  throw InvalidArgException("k", "zeroth order propensity needs one rate constant",k.len);

        return new ZerothOrderPropensity(k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* unitsForConstants[] = {"item/second", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "ZerothOrderPropensity", "k1", unitsForConstants, &create);
    }
};

class FirstOrderPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 1;

    FirstOrderPropensity(uint s, double k) :PropensityFunction(REACTION_TYPE,1),s(s),k(k) {}
    uint s;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        return k * double(speciesCounts[s]);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_load_pd(&speciesCounts[s*DOUBLES_PER_AVX]));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple dependencies = getDependencies(reactionIndex, D);
        if (dependencies.len != 1) throw InvalidArgException("D", "first order propensity had invalid number of dependencies",dependencies.len);

        // Find the rate costant.
        if (k.len < 1)  throw InvalidArgException("k", "first order propensity needs one rate constant",k.len);

        return new FirstOrderPropensity(dependencies[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* unitsForConstants[] = {"1/second", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FirstOrderPropensity", "k1 * x1", unitsForConstants, &create);
    }
};

class SecondOrderPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 2;

    SecondOrderPropensity(uint s1, uint s2, double k) :PropensityFunction(REACTION_TYPE,2),s1(s1),s2(s2),k(k) {}
    uint s1,s2;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=volumeMultiplier;}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        return k * double(speciesCounts[s1]*speciesCounts[s2]);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        return _mm256_mul_pd(_mm256_set1_pd(k),_mm256_mul_pd(_mm256_load_pd(&speciesCounts[s1*DOUBLES_PER_AVX]), _mm256_load_pd(&speciesCounts[s2*DOUBLES_PER_AVX])));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple dependencies = getDependencies(reactionIndex, D);
        if (dependencies.len != 2) throw InvalidArgException("D", "second order propensity had invalid number of dependencies",dependencies.len);

        // Find the rate costant.
        if (k.len < 1)  throw InvalidArgException("k", "second order propensity needs one rate constant",k.len);

        return new SecondOrderPropensity(dependencies[0],dependencies[1],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* unitsForConstants[] = {"1/(item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "SecondOrderPropensity", "k1 * x1 * x2", unitsForConstants, &create);
    }
};

class DimerizationPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 3;

    DimerizationPropensity(uint si, double k) :PropensityFunction(REACTION_TYPE,2),si(si),k(k) {}
    uint si;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=volumeMultiplier;}
    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double s = double(speciesCounts[si]);
        if (s < 2.0) return 0.0;
        return k * s * (s-1.0);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd s = _mm256_load_pd(&speciesCounts[si*DOUBLES_PER_AVX]);
        avxd sm1 = _mm256_sub_pd(s, _mm256_set1_pd(1.0));
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(s,sm1));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple dependencies = getDependencies(reactionIndex, D);
        if (dependencies.len != 1) throw InvalidArgException("D", "dimerization propensity had invalid number of dependencies",dependencies.len);

        // Find the rate costant.
        if (k.len < 1)  throw InvalidArgException("k", "dimerization propensity needs one rate constant",k.len);

        return new DimerizationPropensity(dependencies[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x1", "k1 * x1^2", "k1 * x1 * (x1-1)", NULL};
        const char* unitsForConstants[] = {"1/(item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "DimerizationPropensity", expressions, unitsForConstants, &create);
    }
};

class TrimerizationPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 4;

    TrimerizationPropensity(uint si, double k) :PropensityFunction(REACTION_TYPE,3),si(si),k(k) {}
    uint si;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=pow(volumeMultiplier,2);}

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double s = double(speciesCounts[si]);
        if (s < 3.0) return 0.0;
        return k * s * (s-1.0) * (s-2.0);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd s = _mm256_load_pd(&speciesCounts[si*DOUBLES_PER_AVX]);
        avxd sm1 = _mm256_sub_pd(s, _mm256_set1_pd(1.0));
        avxd sm2 = _mm256_sub_pd(s, _mm256_set1_pd(2.0));
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(s,_mm256_mul_pd(sm1,sm2)));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple dependencies = getDependencies(reactionIndex, D);
        if (dependencies.len != 1) throw InvalidArgException("D", "trimerization propensity had invalid number of dependencies",dependencies.len);

        // Find the rate costant.
        if (k.len < 1)  throw InvalidArgException("k", "trimerization propensity needs one rate constant",k.len);

        return new TrimerizationPropensity(dependencies[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x1 * x1", "k1 * x1^3", "k1 * x1 * (x1-1) * (x1-2)", NULL};
        const char* unitsForConstants[] = {"1/(item*item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "TrimerizationPropensity", expressions, unitsForConstants, &create);
    }
};

class TetramerizationPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 5;

    TetramerizationPropensity(uint si, double k) :PropensityFunction(REACTION_TYPE,4),si(si),k(k) {}
    uint si;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=pow(volumeMultiplier,3);}

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double s = double(speciesCounts[si]);
        if (s < 4.0) return 0.0;
        return k * s * (s-1.0) * (s-2.0) * (s-3.0);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd s = _mm256_load_pd(&speciesCounts[si*DOUBLES_PER_AVX]);
        avxd sm1 = _mm256_sub_pd(s, _mm256_set1_pd(1.0));
        avxd sm2 = _mm256_sub_pd(s, _mm256_set1_pd(2.0));
        avxd sm3 = _mm256_sub_pd(s, _mm256_set1_pd(3.0));
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(s,_mm256_mul_pd(sm1,_mm256_mul_pd(sm2,sm3))));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple dependencies = getDependencies(reactionIndex, D);
        if (dependencies.len != 1) throw InvalidArgException("D", "tetramerization propensity had invalid number of dependencies",dependencies.len);

        // Find the rate costant.
        if (k.len < 1)  throw InvalidArgException("k", "tetramerization propensity needs one rate constant",k.len);

        return new TetramerizationPropensity(dependencies[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x1 * x1 * x1", "k1 * x1^4", "k1 * x1 * (x1-1) * (x1-2) * (x1-3)", NULL};
        const char* unitsForConstants[] = {"1/(item*item*item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "TetramerizationPropensity", expressions, unitsForConstants, &create);
    }
};

class PentamerizationPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 6;

    PentamerizationPropensity(uint si, double k) :PropensityFunction(REACTION_TYPE,5),si(si),k(k) {}
    uint si;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=pow(volumeMultiplier,4);}

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double s = double(speciesCounts[si]);
        if (s < 5.0) return 0.0;
        return k * s * (s-1.0) * (s-2.0) * (s-3.0) * (s-4.0);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd s = _mm256_load_pd(&speciesCounts[si*DOUBLES_PER_AVX]);
        avxd sm1 = _mm256_sub_pd(s, _mm256_set1_pd(1.0));
        avxd sm2 = _mm256_sub_pd(s, _mm256_set1_pd(2.0));
        avxd sm3 = _mm256_sub_pd(s, _mm256_set1_pd(3.0));
        avxd sm4 = _mm256_sub_pd(s, _mm256_set1_pd(4.0));
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(s,_mm256_mul_pd(sm1,_mm256_mul_pd(sm2,_mm256_mul_pd(sm3,sm4)))));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple dependencies = getDependencies(reactionIndex, D);
        if (dependencies.len != 1) throw InvalidArgException("D", "pentamerization propensity had invalid number of dependencies",dependencies.len);

        // Find the rate costant.
        if (k.len < 1)  throw InvalidArgException("k", "pentamerization propensity needs one rate constant",k.len);

        return new PentamerizationPropensity(dependencies[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x1 * x1 * x1 * x1", "k1 * x1^5", "k1 * x1 * (x1-1) * (x1-2) * (x1-3) * (x1-4)", NULL};
        const char* unitsForConstants[] = {"1/(item*item*item*item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "PentamerizationPropensity", expressions, unitsForConstants, &create);
    }
};

class HexamerizationPropensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 7;

    HexamerizationPropensity(uint si, double k) :PropensityFunction(REACTION_TYPE,6),si(si),k(k) {}
    uint si;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=pow(volumeMultiplier,5);}

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double s = double(speciesCounts[si]);
        if (s < 6.0) return 0.0;
        return k * s * (s-1.0) * (s-2.0) * (s-3.0) * (s-4.0) * (s-5.0);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd s = _mm256_load_pd(&speciesCounts[si*DOUBLES_PER_AVX]);
        avxd sm1 = _mm256_sub_pd(s, _mm256_set1_pd(1.0));
        avxd sm2 = _mm256_sub_pd(s, _mm256_set1_pd(2.0));
        avxd sm3 = _mm256_sub_pd(s, _mm256_set1_pd(3.0));
        avxd sm4 = _mm256_sub_pd(s, _mm256_set1_pd(4.0));
        avxd sm5 = _mm256_sub_pd(s, _mm256_set1_pd(5.0));
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(s,_mm256_mul_pd(sm1,_mm256_mul_pd(sm2,_mm256_mul_pd(sm3,_mm256_mul_pd(sm4,sm5))))));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple dependencies = getDependencies(reactionIndex, D);
        if (dependencies.len != 1) throw InvalidArgException("D", "hexamerization propensity had invalid number of dependencies",dependencies.len);

        // Find the rate costant.
        if (k.len < 1)  throw InvalidArgException("k", "hexamerization propensity needs one rate constant",k.len);

        return new HexamerizationPropensity(dependencies[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x1 * x1 * x1 * x1 * x1", "k1 * x1^6", "k1 * x1 * (x1-1) * (x1-2) * (x1-3) * (x1-4) * (x1-5)", NULL};
        const char* unitsForConstants[] = {"1/(item*item*item*item*item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "HexamerizationPropensity", expressions, unitsForConstants, &create);
    }
};


class TwoSpecies21Propensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 103;

    TwoSpecies21Propensity(uint s1, uint s2, double k) :PropensityFunction(REACTION_TYPE,3),s1(s1),s2(s2),k(k) {}
    uint s1, s2;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=(pow(volumeMultiplier,2));}

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double x1 = double(speciesCounts[s1]);
        double x2 = double(speciesCounts[s2]);
        if (x1 < 2.0 || x2 < 1.0) return 0.0;

        return k * x1 * (x1-1.0) * x2;
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd x1 = _mm256_load_pd(&speciesCounts[s1*DOUBLES_PER_AVX]);
        avxd x1m1 = _mm256_sub_pd(x1, _mm256_set1_pd(1.0));
        avxd x2 = _mm256_load_pd(&speciesCounts[s2*DOUBLES_PER_AVX]);
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(_mm256_mul_pd(x1,x1m1),x2));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple d1 = getSpecificDependencies(reactionIndex, D, 1);
        if (d1.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one first species dependency, had %d", d1.len);
        utuple d2 = getSpecificDependencies(reactionIndex, D, 2);
        if (d2.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one second species dependency, had %d", d2.len);

        // Find the rate costant.
        if (k.len < 1)  THROW_EXCEPTION(InvalidArgException, "propensity function needs one rate constant, had %d", k.len);

        return new TwoSpecies21Propensity(d1[0],d2[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x1 * x2", "k1 * x1 * x2 * x1", "k1 * x2 * x1 * x1", "k1 * x1^2 * x2", "k1 * x1 * (x1-1) * x2", NULL};
        const char* unitsForConstants[] = {"1/(item*item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "TwoSpecies21Propensity", expressions, unitsForConstants, &create);
    }
};

class TwoSpecies22Propensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 100;

    TwoSpecies22Propensity(uint s1, uint s2, double k) :PropensityFunction(REACTION_TYPE,4),s1(s1),s2(s2),k(k) {}
    uint s1, s2;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=(pow(volumeMultiplier,3));}

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double x1 = double(speciesCounts[s1]);
        double x2 = double(speciesCounts[s2]);
        if (x1 < 2.0 || x2 < 2.0) return 0.0;

        return k * x1 * (x1-1.0) * x2 * (x2-1.0);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd x1 = _mm256_load_pd(&speciesCounts[s1*DOUBLES_PER_AVX]);
        avxd x1m1 = _mm256_sub_pd(x1, _mm256_set1_pd(1.0));
        avxd x2 = _mm256_load_pd(&speciesCounts[s2*DOUBLES_PER_AVX]);
        avxd x2m1 = _mm256_sub_pd(x2, _mm256_set1_pd(1.0));
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(_mm256_mul_pd(x1,x1m1),_mm256_mul_pd(x2,x2m1)));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple d1 = getSpecificDependencies(reactionIndex, D, 1);
        if (d1.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one first species dependency, had %d", d1.len);
        utuple d2 = getSpecificDependencies(reactionIndex, D, 2);
        if (d2.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one second species dependency, had %d", d2.len);

        // Find the rate costant.
        if (k.len < 1)  THROW_EXCEPTION(InvalidArgException, "propensity function needs one rate constant, had %d", k.len);

        return new TwoSpecies22Propensity(d1[0],d2[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x1 * x2 * x2", "k1 * x1 * x2 * x1 * x2", "k1 * x1 * x2 * x2 * x1", "k1 * x1^2 * x2^2", "k1 * x1 * (x1-1) * x2 * (x2-1)", NULL};
        const char* unitsForConstants[] = {"1/(item*item*item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "TwoSpecies22Propensity", expressions, unitsForConstants, &create);
    }
};


class TwoSpecies32Propensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 101;

    TwoSpecies32Propensity(uint s1, uint s2, double k) :PropensityFunction(REACTION_TYPE,5),s1(s1),s2(s2),k(k) {}
    uint s1, s2;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=(pow(volumeMultiplier,4));}

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double x1 = double(speciesCounts[s1]);
        double x2 = double(speciesCounts[s2]);
        if (x1 < 3.0 || x2 < 2.0) return 0.0;

        return k * x1 * (x1-1.0) * (x1-2.0) * x2 * (x2-1.0);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd x1 = _mm256_load_pd(&speciesCounts[s1*DOUBLES_PER_AVX]);
        avxd x1m1 = _mm256_sub_pd(x1, _mm256_set1_pd(1.0));
        avxd x2 = _mm256_load_pd(&speciesCounts[s2*DOUBLES_PER_AVX]);
        avxd x2m1 = _mm256_sub_pd(x2, _mm256_set1_pd(1.0));
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(_mm256_mul_pd(x1,x1m1),_mm256_mul_pd(x2,x2m1)));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple d1 = getSpecificDependencies(reactionIndex, D, 1);
        if (d1.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one first species dependency, had %d", d1.len);
        utuple d2 = getSpecificDependencies(reactionIndex, D, 2);
        if (d2.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one second species dependency, had %d", d2.len);

        // Find the rate costant.
        if (k.len < 1)  THROW_EXCEPTION(InvalidArgException, "propensity function needs one rate constant, had %d", k.len);

        return new TwoSpecies32Propensity(d1[0],d2[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x1 * x1 * x2 * x2", "k1 * x1^3 * x2^2", "k1 * x1 * (x1-1) * (x1-2) * x2 * (x2-1)", NULL};
        const char* unitsForConstants[] = {"1/(item*item*item*item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "TwoSpecies32Propensity", expressions, unitsForConstants, &create);
    }
};


class TwoSpecies42Propensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 102;

    TwoSpecies42Propensity(uint s1, uint s2, double k) :PropensityFunction(REACTION_TYPE,6),s1(s1),s2(s2),k(k) {}
    uint s1, s2;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=(pow(volumeMultiplier,5));}

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double x1 = double(speciesCounts[s1]);
        double x2 = double(speciesCounts[s2]);
        if (x1 < 4.0 || x2 < 2.0) return 0.0;

        return k * x1 * (x1-1.0) * (x1-2.0) * (x1-3.0) * x2 * (x2-1.0);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd x1 = _mm256_load_pd(&speciesCounts[s1*DOUBLES_PER_AVX]);
        avxd x1m1 = _mm256_sub_pd(x1, _mm256_set1_pd(1.0));
        avxd x2 = _mm256_load_pd(&speciesCounts[s2*DOUBLES_PER_AVX]);
        avxd x2m1 = _mm256_sub_pd(x2, _mm256_set1_pd(1.0));
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(_mm256_mul_pd(x1,x1m1),_mm256_mul_pd(x2,x2m1)));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple d1 = getSpecificDependencies(reactionIndex, D, 1);
        if (d1.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one first species dependency, had %d", d1.len);
        utuple d2 = getSpecificDependencies(reactionIndex, D, 2);
        if (d2.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one second species dependency, had %d", d2.len);

        // Find the rate costant.
        if (k.len < 1)  THROW_EXCEPTION(InvalidArgException, "propensity function needs one rate constant, had %d", k.len);

        return new TwoSpecies42Propensity(d1[0],d2[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x1 * x1 * x1 * x2 * x2", "k1 * x1^4 * x2^2", "k1 * x1 * (x1-1) * (x1-2) * (x1-3) * x2 * (x2-1)", NULL};
        const char* unitsForConstants[] = {"1/(item*item*item*item*item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "TwoSpecies42Propensity", expressions, unitsForConstants, &create);
    }
};

class ThreeSpecies111Propensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 201;

    ThreeSpecies111Propensity(uint s1, uint s2, uint s3, double k) :PropensityFunction(REACTION_TYPE,3),s1(s1),s2(s2),s3(s3),k(k) {}
    uint s1, s2, s3;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=(pow(volumeMultiplier,2));}

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double x1 = double(speciesCounts[s1]);
        double x2 = double(speciesCounts[s2]);
        double x3 = double(speciesCounts[s3]);
        if (x1 < 1.0 || x2 < 1.0 || x3 < 1.0) return 0.0;

        return k * x1 * x2 * x3;
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd x1 = _mm256_load_pd(&speciesCounts[s1*DOUBLES_PER_AVX]);
        avxd x2 = _mm256_load_pd(&speciesCounts[s2*DOUBLES_PER_AVX]);
        avxd x3 = _mm256_load_pd(&speciesCounts[s3*DOUBLES_PER_AVX]);
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(_mm256_mul_pd(x1,x2),x3));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple d1 = getSpecificDependencies(reactionIndex, D, 1);
        if (d1.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one first species dependency, had %d", d1.len);
        utuple d2 = getSpecificDependencies(reactionIndex, D, 2);
        if (d2.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one second species dependency, had %d", d2.len);
        utuple d3 = getSpecificDependencies(reactionIndex, D, 3);
        if (d3.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one third species dependency, had %d", d3.len);

        // Find the rate costant.
        if (k.len < 1)  THROW_EXCEPTION(InvalidArgException, "propensity function needs one rate constant, had %d", k.len);

        return new ThreeSpecies111Propensity(d1[0],d2[0],d3[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x2 * x3", NULL};
        const char* unitsForConstants[] = {"1/(item*item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "ThreeSpecies111Propensity", expressions, unitsForConstants, &create);
    }
};

class ThreeSpecies211Propensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 200;

    ThreeSpecies211Propensity(uint s1, uint s2, uint s3, double k) :PropensityFunction(REACTION_TYPE,4),s1(s1),s2(s2),s3(s3),k(k) {}
    uint s1, s2, s3;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=(pow(volumeMultiplier,3));}

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double x1 = double(speciesCounts[s1]);
        double x2 = double(speciesCounts[s2]);
        double x3 = double(speciesCounts[s3]);
        if (x1 < 2.0 || x2 < 1.0 || x3 < 1.0) return 0.0;

        return k * x1 * (x1-1.0) * x2 * x3;
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd x1 = _mm256_load_pd(&speciesCounts[s1*DOUBLES_PER_AVX]);
        avxd x1m1 = _mm256_sub_pd(x1, _mm256_set1_pd(1.0));
        avxd x2 = _mm256_load_pd(&speciesCounts[s2*DOUBLES_PER_AVX]);
        avxd x3 = _mm256_load_pd(&speciesCounts[s3*DOUBLES_PER_AVX]);
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(_mm256_mul_pd(x1,x1m1),_mm256_mul_pd(x2,x3)));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple d1 = getSpecificDependencies(reactionIndex, D, 1);
        if (d1.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one first species dependency, had %d", d1.len);
        utuple d2 = getSpecificDependencies(reactionIndex, D, 2);
        if (d2.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one second species dependency, had %d", d2.len);
        utuple d3 = getSpecificDependencies(reactionIndex, D, 3);
        if (d3.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one third species dependency, had %d", d3.len);

        // Find the rate costant.
        if (k.len < 1)  THROW_EXCEPTION(InvalidArgException, "propensity function needs one rate constant, had %d", k.len);

        return new ThreeSpecies211Propensity(d1[0],d2[0],d3[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x1 * x2 * x3", "k1 * x2 * x1 * x1 * x3", "k1 * x2 * x3 * x1 * x1", "k1 * x1^2 * x2 * x3", "k1 * x1 * (x1-1) * x2 * x3", NULL};
        const char* unitsForConstants[] = {"1/(item*item*item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "ThreeSpecies211Propensity", expressions, unitsForConstants, &create);
    }
};

class FourSpecies1111Propensity : public lm::me::PropensityFunction
{
public:
    static const uint REACTION_TYPE = 400;

    FourSpecies1111Propensity(uint s1, uint s2, uint s3, uint s4, double k) :PropensityFunction(REACTION_TYPE,4),s1(s1),s2(s2),s3(s3),s4(s4),k(k) {}
    uint s1, s2, s3, s4;
    double k;

    virtual tuple<double> getConstants() {return tuple<double>(k);}
    virtual void setConstants(tuple<double> constants) {k=constants[0];}
    virtual void setConstant(uint index, double constant) {if (index == 0) k=constant; else THROW_EXCEPTION(lm::RuntimeException, "invalid constant index: %d", index);}
    void changeVolume(double volumeMultiplier) {k/=(pow(volumeMultiplier,3));}

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double x1 = double(speciesCounts[s1]);
        double x2 = double(speciesCounts[s2]);
        double x3 = double(speciesCounts[s3]);
        double x4 = double(speciesCounts[s4]);
        if (x1 < 1.0 || x2 < 1.0 || x3 < 1.0 || x4 < 1.0) return 0.0;

        return k * x1 * x2 * x3 * x4;
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd x1 = _mm256_load_pd(&speciesCounts[s1*DOUBLES_PER_AVX]);
        avxd x2 = _mm256_load_pd(&speciesCounts[s2*DOUBLES_PER_AVX]);
        avxd x3 = _mm256_load_pd(&speciesCounts[s3*DOUBLES_PER_AVX]);
        avxd x4 = _mm256_load_pd(&speciesCounts[s4*DOUBLES_PER_AVX]);
        return _mm256_mul_pd(_mm256_set1_pd(k), _mm256_mul_pd(_mm256_mul_pd(x1,x2),_mm256_mul_pd(x3,x4)));
    }
#endif

    static PropensityFunction* create(const uint reactionIndex, const ndarray<int> S, const ndarray<uint> D, const tuple<double>k)
    {
        // Find the species dependencies.
        utuple d1 = getSpecificDependencies(reactionIndex, D, 1);
        if (d1.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one first species dependency, had %d", d1.len);
        utuple d2 = getSpecificDependencies(reactionIndex, D, 2);
        if (d2.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one second species dependency, had %d", d2.len);
        utuple d3 = getSpecificDependencies(reactionIndex, D, 3);
        if (d3.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one third species dependency, had %d", d3.len);
        utuple d4 = getSpecificDependencies(reactionIndex, D, 4);
        if (d4.len != 1) THROW_EXCEPTION(InvalidArgException, "propensity function needs one fourth species dependency, had %d", d4.len);

        // Find the rate costant.
        if (k.len < 1)  THROW_EXCEPTION(InvalidArgException, "propensity function needs one rate constant, had %d", k.len);

        return new FourSpecies1111Propensity(d1[0],d2[0],d3[0],d4[0],k[0]);
    }

    static lm::me::PropensityFunctionDefinition registerFunction()
    {
        const char* expressions[] = {"k1 * x1 * x2 * x3 * x4", NULL};
        const char* unitsForConstants[] = {"1/(item*item*item*second)", NULL};
        return lm::me::PropensityFunctionDefinition(REACTION_TYPE, "FourSpecies1111Propensity", expressions, unitsForConstants, &create);
    }
};


list<lm::me::PropensityFunctionDefinition> CMEPropensityFunctions::getPropensityFunctionDefinitions()
{
    list<lm::me::PropensityFunctionDefinition> defs;
    defs.push_back(ZerothOrderPropensity::registerFunction());
    defs.push_back(FirstOrderPropensity::registerFunction());
    defs.push_back(SecondOrderPropensity::registerFunction());
    defs.push_back(DimerizationPropensity::registerFunction());
    defs.push_back(TrimerizationPropensity::registerFunction());
    defs.push_back(TetramerizationPropensity::registerFunction());
    defs.push_back(PentamerizationPropensity::registerFunction());
    defs.push_back(HexamerizationPropensity::registerFunction());
    defs.push_back(TwoSpecies21Propensity::registerFunction());
    defs.push_back(TwoSpecies22Propensity::registerFunction());
    defs.push_back(TwoSpecies32Propensity::registerFunction());
    defs.push_back(TwoSpecies42Propensity::registerFunction());
    defs.push_back(ThreeSpecies111Propensity::registerFunction());
    defs.push_back(ThreeSpecies211Propensity::registerFunction());
    defs.push_back(FourSpecies1111Propensity::registerFunction());

    return defs;
}


}
}

