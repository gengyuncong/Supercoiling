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

#include <algorithm>
#include <cstdarg>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "lm/Iterator.h"
#include "lm/String.h"
#include "lm/Types.h"

using std::pair;
using std::string;
using std::stringstream;
using std::vector;

namespace lm {


// get the first n characters of a string in a new string
string head(const string& source, size_t length)
{
    return source.substr(0, length);
}

// get the last n characters of a string in a new string (see http://stackoverflow.com/a/7597469/425458)
string tail(const string& source, size_t length)
{
    if (length>=source.size())
    {
        return source;
    }
    return source.substr(source.size() - length);
}

bool endsWith(const string& str, const string & suffix)
{
    return suffix.size() <= str.size() && str.find(suffix, str.size() - suffix.size()) != str.npos;
}

std::string rtrim(std::string s, const char* t)
{
    s.erase(s.find_last_not_of(t) + 1);
    return s;
}

std::string ltrim(std::string s, const char* t)
{
    s.erase(0, s.find_first_not_of(t));
    return s;
}

std::string trim(std::string s, const char* t)
{
    return ltrim(rtrim(s, t), t);
}

string pathJoin(const vector<string>& pathElements, bool absolute)
{
    bool isAbsolute = false;
    stringstream ss;
    vector<string>::const_iterator it=pathElements.begin();

    // find the first non-zero length string in the vector
    while ((it!=pathElements.end()) and (it->size()==0)) it++;

    // check to see if the first non-zero length element is already an absolute path. skip if no non-zero legnth elements (otherwise, segfault due to EXC_BAD_ACCESS)
    if (it!=pathElements.end())
    {
        isAbsolute = ((*it)[0]=='/');
    }

    for (;it!=pathElements.end()--;it++)
    {
//        // see http://stackoverflow.com/a/9359324/425458
//        // By ending at the right iterator, we will do the equivalent of the rstrip operation...
//        string::const_iterator right = std::find_if(it->rbegin(), it->rend(), isNotSlash).base();
//
//        // ...and by starting at the left iterator, we will do the equivalent of the lstrip operation.
//        string::const_iterator left = std::find_if(it->begin(), right, isNotSlash);
//
//        ss << string(left, right);

        // strip any slashes off the left and right sides of the path element
        ss << strip(*it);

        // add a forward slash to the end of the path element if it is not the final one
        if (not isLast(it, pathElements))
        {
            ss << "/";
        }
    }

    string joinedPath(ss.str());
    // the strip ops will have removed any leading "/", so if we want one add it back now
    if ((absolute or isAbsolute) and joinedPath.size() > 0 and joinedPath[0]!='/') joinedPath.insert(0, "/");

    return joinedPath;
}

string pathJoin(const string& elem0, const string& elem1, bool absolute)
{
    std::vector<std::string> elems;
    elems.push_back(elem0);
    elems.push_back(elem1);

    return pathJoin(elems, absolute);
}

string pathName(const string& path)
{
    std::size_t found = path.find_last_of('/');

    if (found != string::npos)
    {
        // at least one slash was found, return everything after the last slash
        return path.substr(found + 1);
    }
    else
    {
        // no slashes found, just return the original path
        return path;
    }
}

string pathWithSuffix(const string& path, const string& suffix)
{
    // find the start of the path suffix, if any, but only in the final element (ie behave like "/foo/bar.re/foobar" doesn't have a suffix)
    std::size_t slashFound, suffixFound;
    slashFound = path.find_last_of('/');

    if (slashFound!=string::npos)
    {
        // at least one slash was found, search for the suffix within the final path element
        suffixFound = path.substr(slashFound + 1).find_last_of('.');

        // convert from pos in the substring to pos in the original path
        if (suffixFound!=string::npos)
        {
            suffixFound+=slashFound + 1;
        }
    }
    else
    {
        // no slashes found, just search the whole path for a suffix
        suffixFound = path.find_last_of('.');
    }

    string newPath;
    if (suffixFound!=string::npos)
    {
        // existing suffix was found, strip it off
        newPath = path.substr(0, suffixFound);
    }
    else
    {
        // no existing suffix, just copy path
        newPath = path;
    }

    // add on the new suffix
    newPath.append(suffix);
    return newPath;
}

// "lambda" function needed for strip
bool _isNotSlash(const char& c) {return c!='/';}

string strip(const string& str)
{
    // see http://stackoverflow.com/a/9359324/425458
    // By ending at the right iterator, we will do the equivalent of the rstrip operation...
    string::const_iterator right = std::find_if(str.rbegin(), str.rend(), _isNotSlash).base();

    // ...and by starting at the left iterator, we will do the equivalent of the lstrip operation.
    string::const_iterator left = std::find_if(str.begin(), right, _isNotSlash);

    return string(left, right);
}

int32_t parseInt(const string& arg)
{
    return static_cast<int32_t>(atoi(arg.c_str()));
}

uint32_t parseUint(const string& arg)
{
    return static_cast<uint32_t>(atoi(arg.c_str()));
}

vector<string> parseStringList(const string& arg, const char* delimiters)
{
    vector<string> list;
    char * argbuf = new char[strlen(arg.c_str())+1];
    strcpy(argbuf,arg.c_str());
    char * pch = strtok(argbuf,delimiters);
    while (pch != NULL)
    {
        if (strlen(pch) > 0) list.push_back(string(pch));
        pch = strtok(NULL,delimiters);
    }
    delete[] argbuf;
    return list;
}

vector<int32_t> parseIntList(const string& arg, const char* delimiters)
{
    vector<int32_t> list;
    char * argbuf = new char[strlen(arg.c_str())+1];
    strcpy(argbuf,arg.c_str());
    char * pch = strtok(argbuf,delimiters);
    while (pch != NULL)
    {
        char * rangeDelimiter;
        if ((rangeDelimiter=strstr(pch,"-")) != NULL)
        {
            *rangeDelimiter='\0';
            int32_t begin=int32_t(atoi(pch));
            int32_t end=int32_t(atoi(rangeDelimiter+1));
            for (int32_t i=begin; i<=end; i++)
                list.push_back(i);
        }
        else
        {
            if (strlen(pch) > 0) list.push_back(int32_t(atoi(pch)));
        }
        pch = strtok(NULL,delimiters);
    }
    delete[] argbuf;
    return list;
}

vector<uint32_t> parseUintList(const string& arg, const char* delimiters)
{
    vector<uint32_t> list;
    char * argbuf = new char[strlen(arg.c_str())+1];
    strcpy(argbuf,arg.c_str());
    char * pch = strtok(argbuf,delimiters);
    while (pch != NULL)
    {
        char * rangeDelimiter;
        if ((rangeDelimiter=strstr(pch,"-")) != NULL)
        {
            *rangeDelimiter='\0';
            uint32_t begin=uint32_t(atoi(pch));
            uint32_t end=uint32_t(atoi(rangeDelimiter+1));
            for (uint32_t i=begin; i<=end; i++)
                list.push_back(i);
        }
        else
        {
            if (strlen(pch) > 0) list.push_back(uint32_t(atoi(pch)));
        }
        pch = strtok(NULL,delimiters);
    }
    delete[] argbuf;
    return list;
}

utuple parseUintRange(const string& str, uint32_t maxDim)
{
    utuple ret(0U,0U);

    char* s = new char[str.size()+1];
    strcpy(s, str.c_str());

    if (strlen(s) == 1 && s[0] == ':')
    {
        ret = utuple(0,maxDim);
    }
    else if (strlen(s) >= 2 && s[0] == ':')
    {
        ret = utuple(0,static_cast<uint32_t>(atoi(&(s[1]))));
    }
    else if (strlen(s) >= 2 && s[strlen(s)-1] == ':')
    {
        s[strlen(s)-1] = 0;
        ret = utuple(static_cast<uint32_t>(atoi(&(s[1]))), maxDim);
    }
    else if (strlen(s) >= 1)
    {
        char * split = strstr(s, ":");
        if (split > s && split < s+strlen(s)-1)
        {
            char * start = s;
            *split = 0;
            char * end = split+1;
            ret = utuple(static_cast<uint32_t>(atoi(start)),static_cast<uint32_t>(atoi(end)));
        }
        else
        {
            ret = utuple(static_cast<uint32_t>(atoi(s)),static_cast<uint32_t>(atoi(s)));
        }
    }

    delete[] s;
    return ret;
}

vector<utuple> parseUintIndices(const string& str)
{
    vector<utuple> indices;

    char * sbuffer = new char[str.size()+1];
    strcpy(sbuffer, str.c_str());
    int sbufferStart=0;

    // Remove any leading/trailing parentheses.
    if (sbuffer[0] == '(') sbufferStart++;
    if (sbuffer[strlen(sbuffer)-1] == ')') sbuffer[strlen(sbuffer)-1] = '\0';

    // Parse each dimension.
    char * pch = strtok(&(sbuffer[sbufferStart]), ",");
    while (pch != NULL)
    {
        indices.push_back(parseUintRange(string(pch)));
        pch = strtok(NULL, ",");
    }
    delete [] sbuffer;
    return indices;
}

vector<string> parseStringValues(const string& s)
{
    char * valueToks = new char[s.size()+1];
    strcpy(valueToks, s.c_str());
    int valueToksStart=0;

    // Remove any leading/trailing brackets.
    if (valueToks[0] == '[') valueToksStart++;
    if (valueToks[s.size()-1] == ']') valueToks[s.size()-1] = '\0';

    // Parse the values.
    vector<string> values;
    char * value = strtok(valueToks+valueToksStart, ";,");
    while (value != NULL)
    {
        string v = value;
        values.push_back(v);
        value = strtok(NULL, ";,");
    }
    delete [] valueToks;
    return values;
}


}
