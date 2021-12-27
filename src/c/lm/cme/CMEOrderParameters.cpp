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

#include <cmath>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "lm/ClassFactory.h"
#include "lm/Types.h"
#include "lm/cme/CMEOrderParameters.h"
#include "lm/oparam/OrderParameterFunction.h"
#include "lm/types/OrderParameters.pb.h"

using std::list;
using std::map;
using std::string;
using std::vector;

namespace lm {
namespace cme {

bool CMEOrderParameters::registered=CMEOrderParameters::registerClass();

bool CMEOrderParameters::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::oparam::OrderParameterFunctionCollection", "lm::cme::CMEOrderParameters", &CMEOrderParameters::allocateObject);
    return true;
}

void* CMEOrderParameters::allocateObject()
{
    return new CMEOrderParameters();
}

CMEOrderParameters::CMEOrderParameters()
{
}

CMEOrderParameters::~CMEOrderParameters()
{
}


class LinearCombinationOrderParameter : public lm::oparam::OrderParameterFunction
{
public:
    static const uint OPARAM_TYPE = 0;

    LinearCombinationOrderParameter(size_t size, uint* speciesIndex, double* speciesCoefficient, uint oparamType=OPARAM_TYPE)
    :OrderParameterFunction(oparamType),size(size),speciesIndex(speciesIndex),speciesCoefficient(speciesCoefficient) {}
    ~LinearCombinationOrderParameter()
    {
        if (speciesIndex != NULL) delete[] speciesIndex; speciesIndex = NULL;
        if (speciesCoefficient != NULL) delete[] speciesCoefficient; speciesCoefficient = NULL;
    }

    size_t size;
    uint* speciesIndex;
    double* speciesCoefficient;

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double value=0.0;
        for (size_t i=0; i<size; i++)
            value += speciesCoefficient[i]*double(speciesCounts[speciesIndex[i]]);
        return value;
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd value = _mm256_set1_pd(0.0);
        for (size_t i=0; i<size; i++)
            value = _mm256_add_pd(value, _mm256_mul_pd(_mm256_set1_pd(speciesCoefficient[i]), _mm256_load_pd(&speciesCounts[speciesIndex[i]*DOUBLES_PER_AVX])));
        return value;
    }
#endif

    static OrderParameterFunction* create(const lm::types::OrderParameter& op)
    {
        if (op.type() != OPARAM_TYPE)
            throw lm::InvalidArgException("op.type", "Mismatch of types during creation of linear order parameter function",op.type(), OPARAM_TYPE);
        if (op.species_index_size() != op.species_coefficient_size())
            throw lm::InvalidArgException("op.size", "Mismatch of sizes during creation of linear order parameter function",op.species_index_size(), op.species_coefficient_size());

        size_t size = op.species_index_size();
        uint* speciesIndex = new uint[size];
        double* speciesCoefficient = new double[size];
        for (size_t i=0; i<size; i++)
        {
            speciesIndex[i] = op.species_index(i);
            speciesCoefficient[i] = op.species_coefficient(i);
        }

        return new LinearCombinationOrderParameter(size, speciesIndex, speciesCoefficient);
    }

    static lm::oparam::OrderParameterFunctionDefinition registerFunction()
    {
        return lm::oparam::OrderParameterFunctionDefinition(OPARAM_TYPE, &create);
    }
};

class TwoSpeciesOrderParameter : public lm::oparam::OrderParameterFunction
{
public:
    static const uint OPARAM_TYPE = 2;

    TwoSpeciesOrderParameter(uint s1, uint s2, double k1, double k2):OrderParameterFunction(OPARAM_TYPE),s1(s1),s2(s2),k1(k1),k2(k2) {}
    uint s1, s2;
    double k1, k2;

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        return k1*double(speciesCounts[s1]) + k2*double(speciesCounts[s2]);
    }

#ifdef OPT_AVX
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd value = _mm256_mul_pd(_mm256_set1_pd(k1), _mm256_load_pd(&speciesCounts[s1*DOUBLES_PER_AVX]));
        return _mm256_fmadd_pd(_mm256_set1_pd(k2), _mm256_load_pd(&speciesCounts[s2*DOUBLES_PER_AVX]), value);
    }
#endif

    static OrderParameterFunction* create(const lm::types::OrderParameter& op)
    {
        if (op.type() != OPARAM_TYPE)
            throw lm::InvalidArgException("op.type", "Mismatch of types during creation of two species order parameter function",op.type(), OPARAM_TYPE);
        if (op.species_index_size() != 2 && op.species_coefficient_size() != 2)
            throw lm::InvalidArgException("op.size", "Mismatch of sizes during creation of two species order parameter function",op.species_index_size(), op.species_coefficient_size());

        return new TwoSpeciesOrderParameter(op.species_index(0), op.species_index(1), op.species_coefficient(0), op.species_coefficient(1));
    }

    static lm::oparam::OrderParameterFunctionDefinition registerFunction()
    {
        return lm::oparam::OrderParameterFunctionDefinition(OPARAM_TYPE, &create);
    }
};


class PolynomialOrderParameter : public lm::oparam::OrderParameterFunction
{
public:
    static const uint OPARAM_TYPE = 100;

    PolynomialOrderParameter(size_t size, uint* speciesIndex, double* speciesCoefficient, double* speciesExponent)
    :OrderParameterFunction(OPARAM_TYPE),size(size),speciesIndex(speciesIndex),speciesCoefficient(speciesCoefficient),
     speciesExponent(speciesExponent)
    {
    }

    ~PolynomialOrderParameter()
    {
        if (speciesIndex != NULL) delete[] speciesIndex; speciesIndex = NULL;
        if (speciesCoefficient != NULL) delete[] speciesCoefficient; speciesCoefficient = NULL;
        if (speciesExponent != NULL) delete[] speciesExponent; speciesExponent = NULL;
    }

    double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const
    {
        double value=0.0;
        for (size_t i=0; i<size; i++)
            value += speciesCoefficient[i]*pow(double(speciesCounts[speciesIndex[i]]), speciesExponent[i]);
        return value;
    }

#if defined(OPT_AVX) && !defined(OPT_SVML)
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd value = _mm256_set1_pd(0.0);
        for (size_t i=0; i<size; i++)
        {
            size_t specIx = speciesIndex[i];
            double specExp = speciesExponent[i];

            // exponentiate. NB: Order is reversed by the AVX _set commands
            avxd xexp = _mm256_set_pd(pow(speciesCounts[specIx*DOUBLES_PER_AVX+3], specExp),
                                      pow(speciesCounts[specIx*DOUBLES_PER_AVX+2], specExp),
                                      pow(speciesCounts[specIx*DOUBLES_PER_AVX+1], specExp),
                                      pow(speciesCounts[specIx*DOUBLES_PER_AVX], specExp));

            // multiply and accumulate in value
            value = _mm256_add_pd(
                value, _mm256_mul_pd(
                    _mm256_set1_pd(speciesCoefficient[i]), xexp
                )
            );

        }
        return value;
    }
#elif defined(OPT_AVX) && defined(OPT_SVML)
    avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const
    {
        avxd value = _mm256_set1_pd(0.0);
        for (size_t i=0; i<size; i++)
        {
            // exponentiate, multiply, and accumulate in value
            value = _mm256_add_pd(
                value, _mm256_mul_pd(
                    _mm256_set1_pd(speciesCoefficient[i]), _mm256_pow_pd(
                        _mm256_load_pd(&speciesCounts[speciesIndex[i]*DOUBLES_PER_AVX]), _mm256_set1_pd(speciesExponent[i])
                    )
                )
            );

        }
        return value;
    }
#endif

    static OrderParameterFunction* create(const lm::types::OrderParameter& op)
    {
        if (op.type() != OPARAM_TYPE)
            throw lm::InvalidArgException("op.type", "Mismatch of types during creation of polynomial order parameter function", op.type(), OPARAM_TYPE);
        if (op.species_index_size() != op.species_coefficient_size() or op.species_index_size() != op.species_exponent_size())
            throw lm::InvalidArgException("op.size", "Mismatch of sizes during creation of polynomial order parameter function", op.species_index_size(), op.species_coefficient_size(), op.species_exponent_size());

        size_t size = op.species_index_size();
        uint* speciesIndex = new uint[size];
        double* speciesCoefficient = new double[size];
        double* speciesExponent = new double[size];

        for (size_t i=0; i<size; i++)
        {
            speciesIndex[i] = op.species_index(i);
            speciesCoefficient[i] = op.species_coefficient(i);
            speciesExponent[i] = op.species_exponent(i);
        }

        return new PolynomialOrderParameter(size, speciesIndex, speciesCoefficient, speciesExponent);
    }

    static lm::oparam::OrderParameterFunctionDefinition registerFunction()
    {
        return lm::oparam::OrderParameterFunctionDefinition(OPARAM_TYPE, &create);
    }

    size_t size;
    uint* speciesIndex;
    double* speciesCoefficient;
    double* speciesExponent;
};

list<lm::oparam::OrderParameterFunctionDefinition> CMEOrderParameters::getOrderParameterFunctionDefinitions()
{
    list<lm::oparam::OrderParameterFunctionDefinition> defs;
    defs.push_back(LinearCombinationOrderParameter::registerFunction());
    defs.push_back(TwoSpeciesOrderParameter::registerFunction());
    defs.push_back(PolynomialOrderParameter::registerFunction());
    return defs;
}



}
}
