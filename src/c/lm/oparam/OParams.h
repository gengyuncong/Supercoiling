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
#ifndef LM_OPARAM_OPARAMLIST
#define LM_OPARAM_OPARAMLIST

//#include <google/protobuf/repeated_field.h>
//#include <iterator>
//#include <map>
//#include <string>
//#include <vector>

//#include "lm/types/OrderParameters.pb.h"
//#include "lm/oparam/OParam.h"

//namespace lm {
//namespace oparam {

//typedef std::map<uint,std::string> OPClassMap;
//typedef google::protobuf::RepeatedPtrField<lm::types::OrderParameter>::const_iterator OPIterator;
//typedef std::map<uint,lm::oparam::OParam*> OPMap;

//class OParams
//{
//public:
//    // constructors/destructors/initializers
//    OParams();
//    ~OParams();
//    void clearOPMap();
//    void init(const lm::types::OrderParameters& oparams);
//    void init();
//    void initOParam(const lm::types::OrderParameter& oparam);
//    void initValues(int* speciesCounts, double time);

//    // operators
//    lm::oparam::OParam* operator[](uint i) {return opMap[i];}
//    OPMap::iterator begin() {return opMap.begin();}
//    OPMap::iterator end() {return opMap.end();}

//    // accessors
//    const lm::oparam::OParam* at(const uint& key) const {return opMap.at(key);}
//    uint64_t size() const {return opMap.size();}

//    // mutators
//    void calcAndStore(int* speciesCounts, double time);
//    lm::types::OrderParameters* getOParamsBuf() {return &oparamsBuf;}
//    void setOParamsBuf(const lm::types::OrderParameters& newOParamsBuf) {*getOParamsBuf() = newOParamsBuf;}

//    // static methods
//    static OPClassMap opClassMap;
//    static OPClassMap makeOPClassMap()
//    {
//        std::map<uint,std::string> m;
//        m[0] = "lm::oparam::OParamLinear";
//        m[2] = "lm::oparam::OParamTwoSpecies";
//        m[100] = "lm::oparam::OParamPolynomial";
//        m[1000] = "lm::oparam::OParamBasins";
//        // m[9999...] = "lm::oparam::OParamTranscendental";
//        return m;
//    }
//protected:
//    lm::types::OrderParameters oparamsBuf;
//    OPMap opMap;
//};

//}
//}

#endif /* LM_OPARAM_OPARAMLIST */
