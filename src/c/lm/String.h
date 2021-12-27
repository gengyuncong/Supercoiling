/*
 * Copyright 2019 Johns Hopkins University
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

#ifndef LM_STRING_H_
#define LM_STRING_H_

#include <cstdio>
#include <google/protobuf/message.h>
#include <google/protobuf/repeated_field.h>
#include <ostream>
#include <string>
#include <sstream>
#include <vector>

#include "robertslab/Types.h"

namespace lm {

// get the first n characters of a string in a new string
std::string head(const std::string& source, size_t length);

// get the last n characters of a string in a new string (see http://stackoverflow.com/a/7597469/425458)
std::string tail(const std::string& source, size_t length);

// See if a string ends with a given suffix.
bool endsWith(std::string const& str, std::string const& suffix);

// Trim a string of whitespace characters.
std::string rtrim(std::string s, const char* t=" \t\n\r");
std::string ltrim(std::string s, const char* t=" \t\n\r");
std::string trim(std::string s, const char* t=" \t\n\r");

// join a vector of path elements into a "/" delineated path.
// If absolute, ensures that there is exactly one "/" at the beginning of the path, otherwise any leading "/" are stripped
std::string pathJoin(const std::vector<std::string>& pathElements, bool absolute=true);

// convenience overloads for pathJoin
std::string pathJoin(const std::string& elem0, const std::string& elem1, bool absolute=true);// {std::vector<std::string> elems; elems.push_back(elem0); elems.push_back(elem1); return pathJoin(elems, absolute);}

// get just the final element from a path.
// example: pathName("foo/bar/re")=="re"
std::string pathName(const std::string& path);

// remove the suffix from the final element of a path, then append the newSuffix.
// example: pathWithSuffix("rey/far/foo.bar", ".rab")=="rey/far/foo.rab", pathWithSuffix("rey/far.bar/foo", ".rab")=="rey/far.bar/foo.rab"
std::string pathWithSuffix(const std::string& path, const std::string& newSuffix);

// strip any "/" from the left and right of str
std::string strip(const std::string& str);

// Parse an int.
int32_t parseInt(const std::string& arg);

// Parse a uint.
uint32_t parseUint(const std::string& arg);

// Parse a list of strings into a vector.
std::vector<std::string> parseStringList(const std::string& arg, const char* delimiters=" ,;:");

// Parse a list of signed integers or ranges of unsigned integers into a vector.
std::vector<int32_t> parseIntList(const std::string& arg, const char* delimiters=" ,;:");

// Parse a list of unsigned integers or ranges of unsigned integers into a vector.
std::vector<uint32_t> parseUintList(const std::string& arg, const char* delimiters=" ,;:");

utuple parseUintRange(const std::string& str, uint32_t maxDim=0xFFFFFFFF);

std::vector<utuple> parseUintIndices(const std::string& str);

std::vector<std::string> parseStringValues(const std::string& s);

}

#endif /* LM_PRINT_H_ */
