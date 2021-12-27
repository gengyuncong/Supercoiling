/*
 * University of Illinois Open Source License
 * Copyright 2012-2014 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Roberts Group
 *               Johns Hopkins University
 *               http://biophysics.jhu.edu/roberts/
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
//#include "lm/ClassFactory.h"
//#include "lm/types/OrderParameters.pb.h"
//#include "lm/io/TrajectoryState.pb.h"
//#include "lm/oparam/OParam.h"

//namespace lm {
//namespace oparam {

//// base class OParam methods
//OParam::OParam(): op(NULL), opFunc(NULL), val(0), prevVal(0)
//{
//}

//OParam::~OParam()
//{
//    if (op!=NULL) delete op; op = NULL;
//    if (opFunc!=NULL) delete opFunc; opFunc = NULL;
//}

//void OParam::init(const lm::types::OrderParameter& opRef)
//{
//    op = new lm::types::OrderParameter(opRef);

//    OrderParameterFunctionFactory fs;
//    opFunc = fs.createOrderParameterFunction(*op);
//}

//void OParam::initValues(int* speciesCounts, double time)
//{
//    val = calc(speciesCounts, time);
//    prevVal = val;
//}

//double OParam::calcAndStore(int* speciesCounts, double time)
//{
//    prevVal = val;
//    val = calc(speciesCounts, time);
//    return val;
//}

//double OParam::calc(const lm::io::TrajectoryState& state) const
//{
//    const lm::io::SpeciesCounts& sc(state.cme_state().species_counts());

//    // offset the species_count pointer to ensure that we only get the last "row" of values
//    int offset = (sc.number_entries() - 1)*(sc.number_species());
//    return calc(sc.species_count().data() + offset, sc.time(sc.number_entries() - 1));
//}

//// derived class methods
//bool OParamLinear::registered=OParamLinear::registerClass();
//bool OParamLinear::registerClass()
//{
//    lm::ClassFactory::getInstance().registerClass("lm::oparam::OParam","lm::oparam::OParamLinear",&OParamLinear::allocateObject);
//    return true;
//}
//void* OParamLinear::allocateObject()
//{
//    return new OParamLinear();
//}

//OParamLinear::OParamLinear(): OParam(), size(), speciesID(), speciesCoefficient() {}

//void OParamLinear::init(const lm::types::OrderParameter& opRef)
//{
//    // call parent method
//    OParam::init(opRef);
//    size = op->species_index_size();
//    speciesID = op->species_index().data();
//    speciesCoefficient = op->species_coefficient().data();
//}

//bool OParamTwoSpecies::registered=OParamTwoSpecies::registerClass();
//bool OParamTwoSpecies::registerClass()
//{
//    lm::ClassFactory::getInstance().registerClass("lm::oparam::OParam","lm::oparam::OParamTwoSpecies",&OParamTwoSpecies::allocateObject);
//    return true;
//}
//void* OParamTwoSpecies::allocateObject()
//{
//    return new OParamTwoSpecies();
//}

//OParamTwoSpecies::OParamTwoSpecies(): OParam() {}

//void OParamTwoSpecies::init(const lm::types::OrderParameter& opRef)
//{
//    // call parent method
//    OParam::init(opRef);
//    s1 = op->species_index(0);
//    s2 = op->species_index(1);
//    k1 = op->species_coefficient(0);
//    k2 = op->species_coefficient(1);
//}

//bool OParamPolynomial::registered=OParamPolynomial::registerClass();
//bool OParamPolynomial::registerClass()
//{
//    lm::ClassFactory::getInstance().registerClass("lm::oparam::OParam","lm::oparam::OParamPolynomial",&OParamPolynomial::allocateObject);
//    return true;
//}
//void* OParamPolynomial::allocateObject()
//{
//    return new OParamLinear();
//}

//OParamPolynomial::OParamPolynomial(): OParam(), size(), speciesID(), speciesCoefficient(), speciesExponents() {}

//void OParamPolynomial::init(const lm::types::OrderParameter& opRef)
//{
//    // call parent method
//    OParam::init(opRef);
//    size = op->species_index_size();
//    speciesID = op->species_index().data();
//    speciesCoefficient = op->species_coefficient().data();
//    speciesExponents = op->species_exponent().data();
//}

//bool OParamBasins::registered=OParamBasins::registerClass();
//bool OParamBasins::registerClass()
//{
//    lm::ClassFactory::getInstance().registerClass("lm::oparam::OParam","lm::oparam::OParamBasins",&OParamBasins::allocateObject);
//    return true;
//}
//void* OParamBasins::allocateObject()
//{
//    return new OParamBasins();
//}

//OParamBasins::OParamBasins(): OParam(), size(), speciesID(), speciesCoefficient() {}

//void OParamBasins::init(const lm::types::OrderParameter& opRef)
//{
//    // call parent method
//    OParam::init(opRef);
//    size = op->basins(0).species_count_size();
//    for (int i=0; i<size; i++)
//    {
//        op->add_species_index(i);
//        op->add_species_coefficient(op->basins(1).species_count(i) - op->basins(0).species_count(i));
//    }
//    speciesID = op->species_index().data();
//    speciesCoefficient = op->species_coefficient().data();
//}

//}
//}
