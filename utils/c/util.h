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

#ifndef LM_UTIL_H_
#define LM_UTIL_H_

#include <string>
#include <utility>
#include <vector>

#include "robertslab/Types.h"

using std::pair;
using std::string;
using std::vector;


vector<int> parseIndices(string s, vector<int> matrixDims);
vector<utuple> parseIndicesAsTuple(string s, vector<int> matrixDims);
pair<int,int> parseRange(char *, int maxDim);
vector<double> parseValues(string s);
vector<int32_t> parseIntValues(string s);
vector<uint32_t> parseUintValues(string s);

#endif

