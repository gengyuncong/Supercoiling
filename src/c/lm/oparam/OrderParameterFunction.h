/*
 * Copyright 2012-2019 Johns Hopkins University
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
 * Author(s): Elijah Roberts, Max Klein
 */

#ifndef LM_OPARAM_ORDERPARAMETERFUNCTION_H
#define LM_OPARAM_ORDERPARAMETERFUNCTION_H

#include <limits>
#include <list>
#include <map>
#include <string>
#include <vector>

#include "lm/Types.h"
#include "lm/types/OrderParameters.pb.h"

using std::list;
using std::map;
using std::string;
using std::vector;


namespace lm {
namespace oparam {

// The base class for any propensity function.
class OrderParameterFunction
{
public:
    OrderParameterFunction(const uint type):type(type){}
    virtual ~OrderParameterFunction() {}
    uint getType() const {return type;}
    virtual double calculate(const ndarray<int32_t>& speciesCounts) const;
    virtual double calculate(const double time, const ndarray<int32_t>& speciesCounts) const;
    virtual double calculate(const int* speciesCounts, const uint numberSpecies) const;
    virtual double calculate(const double time, const int* speciesCounts, const uint numberSpecies) const=0;
#ifdef OPT_AVX
    virtual avxd calculateAvx(const avxd time, const double* speciesCounts, const uint numberSpecies) const;
#endif

protected:
    const uint type;
};

// The type definition for a function to create the propensity function.
typedef OrderParameterFunction* (*OrderParameterFunctionCreator)(const lm::types::OrderParameter& msg);

struct OrderParameterFunctionDefinition
{
    OrderParameterFunctionDefinition():type(std::numeric_limits<uint>::max()),create(NULL){}
    OrderParameterFunctionDefinition(uint type, OrderParameterFunctionCreator create):type(type),create(create){}
    OrderParameterFunctionDefinition(const OrderParameterFunctionDefinition& p):type(p.type),create(p.create){}
    uint type;
    OrderParameterFunctionCreator create;
};

class OrderParameterFunctionFactory
{
public:
    OrderParameterFunctionFactory();
    ~OrderParameterFunctionFactory();
    OrderParameterFunction* createOrderParameterFunction(const lm::types::OrderParameter& msg);

private:
    map<uint,OrderParameterFunctionDefinition> functions;
};

// The base class for a collection of propensity functions.
class OrderParameterFunctionCollection
{
public:
    OrderParameterFunctionCollection();
    virtual ~OrderParameterFunctionCollection();
    virtual list<OrderParameterFunctionDefinition> getOrderParameterFunctionDefinitions()=0;
};

}
}

#endif // LM_OPARAM_ORDERPARAMETERFUNCTION_H
