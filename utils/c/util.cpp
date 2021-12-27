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

#include <list>
#include <string>
#include <utility>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include "lm/Exceptions.h"
#include "robertslab/Types.h"

#include "util.h"

using std::list;
using std::pair;
using std::string;
using std::vector;

vector<int> parseIndices(string s, vector<int> matrixDims)
{
	vector<int> indices;

	if (s.size() == 0)
	{
		if (matrixDims.size() == 1)
		{
			for (int i=0; i<matrixDims[0]; i++)
				indices.push_back(i);
		}
		else if (matrixDims.size() == 2)
		{
			for (int i=0; i<matrixDims[0]; i++)
				for (int j=0; j<matrixDims[1]; j++)
					indices.push_back(i*matrixDims[1]+j);
		}
		else if (matrixDims.size() == 3)
		{
			for (int i=0; i<matrixDims[0]; i++)
				for (int j=0; j<matrixDims[1]; j++)
					for (int k=0; k<matrixDims[2]; k++)
						indices.push_back(i*matrixDims[1]*matrixDims[2]+j*matrixDims[2]+k);
		}
		else
			throw lm::Exception("Unsupported matrix dimension", matrixDims.size());
	}
	else
	{
		char * sbuffer = new char[s.size()+1];
		strcpy(sbuffer, s.c_str());
		int sbufferStart=0;

		// Remove any leading/trailing parentheses.
		if (sbuffer[0] == '(') sbufferStart++;
		if (sbuffer[s.size()-1] == ')') sbuffer[s.size()-1] = '\0';

		if (matrixDims.size() == 1)
		{
			pair<int,int> range = parseRange(&(sbuffer[sbufferStart]), matrixDims[0]);
			for (int i=range.first; i<=range.second; i++)
				indices.push_back(i);
		}
		else if (matrixDims.size() == 2)
		{
			char * istr = &(sbuffer[sbufferStart]);
			char * jstr = (char *)"0";
			char * split = strstr(istr, ",");
			if (split != NULL)
			{
				*split = '\0';
				jstr = split+1;
			}
			pair<int,int> irange = parseRange(istr, matrixDims[0]);
			pair<int,int> jrange = parseRange(jstr, matrixDims[1]);
			for (int i=irange.first; i<=irange.second; i++)
				for (int j=jrange.first; j<=jrange.second; j++)
					indices.push_back(i*matrixDims[1]+j);
		}
		else if (matrixDims.size() == 3)
		{
			char * istr = &(sbuffer[sbufferStart]);
			char * jstr = (char *)"0";
			char * kstr = (char *)"0";
			char * split = strstr(istr, ",");
			if (split != NULL)
			{
				*split = '\0';
				jstr = split+1;
				split = strstr(jstr, ",");
				if (split != NULL)
				{
					*split = '\0';
					kstr = split+1;
				}
			}
			pair<int,int> irange = parseRange(istr, matrixDims[0]);
			pair<int,int> jrange = parseRange(jstr, matrixDims[1]);
			pair<int,int> krange = parseRange(kstr, matrixDims[2]);
			for (int i=irange.first; i<=irange.second; i++)
				for (int j=jrange.first; j<=jrange.second; j++)
					for (int k=krange.first; k<=krange.second; k++)
						indices.push_back(i*matrixDims[1]*matrixDims[2]+j*matrixDims[2]+k);
		}
		else
			throw lm::Exception("Unsupported matrix dimension", matrixDims.size());

		delete [] sbuffer;
	}

	return indices;
}

vector<utuple> parseIndicesAsTuple(string s, vector<int> matrixDims)
{
    vector<utuple> indices;

    if (s.size() == 0)
    {
        if (matrixDims.size() == 1)
        {
            for (int i=0; i<matrixDims[0]; i++)
                indices.push_back(utuple(static_cast<uint>(i)));
        }
        else if (matrixDims.size() == 2)
        {
            for (int i=0; i<matrixDims[0]; i++)
                for (int j=0; j<matrixDims[1]; j++)
                    indices.push_back(utuple(static_cast<uint>(i),static_cast<uint>(j)));
        }
        else if (matrixDims.size() == 3)
        {
            for (int i=0; i<matrixDims[0]; i++)
                for (int j=0; j<matrixDims[1]; j++)
                    for (int k=0; k<matrixDims[2]; k++)
                        indices.push_back(utuple(static_cast<uint>(i),static_cast<uint>(j),static_cast<uint>(k)));
        }
        else
            throw lm::Exception("Unsupported matrix dimension", static_cast<int>(matrixDims.size()));
    }
    else
    {
        char * sbuffer = new char[s.size()+1];
        strcpy(sbuffer, s.c_str());
        int sbufferStart=0;

        // Remove any leading/trailing parentheses.
        if (sbuffer[0] == '(') sbufferStart++;
        if (sbuffer[s.size()-1] == ')') sbuffer[s.size()-1] = '\0';

        if (matrixDims.size() == 1)
        {
            pair<int,int> range = parseRange(&(sbuffer[sbufferStart]), matrixDims[0]);
            for (int i=range.first; i<=range.second; i++)
                indices.push_back(utuple(static_cast<uint>(i)));
        }
        else if (matrixDims.size() == 2)
        {
            char * istr = &(sbuffer[sbufferStart]);
            char * jstr = (char *)"0";
            char * split = strstr(istr, ",");
            if (split != NULL)
            {
                *split = '\0';
                jstr = split+1;
            }
            pair<int,int> irange = parseRange(istr, matrixDims[0]);
            pair<int,int> jrange = parseRange(jstr, matrixDims[1]);
            for (int i=irange.first; i<=irange.second; i++)
                for (int j=jrange.first; j<=jrange.second; j++)
                    indices.push_back(utuple(static_cast<uint>(i),static_cast<uint>(j)));
        }
        else if (matrixDims.size() == 3)
        {
            char * istr = &(sbuffer[sbufferStart]);
            char * jstr = (char *)"0";
            char * kstr = (char *)"0";
            char * split = strstr(istr, ",");
            if (split != NULL)
            {
                *split = '\0';
                jstr = split+1;
                split = strstr(jstr, ",");
                if (split != NULL)
                {
                    *split = '\0';
                    kstr = split+1;
                }
            }
            pair<int,int> irange = parseRange(istr, matrixDims[0]);
            pair<int,int> jrange = parseRange(jstr, matrixDims[1]);
            pair<int,int> krange = parseRange(kstr, matrixDims[2]);
            for (int i=irange.first; i<=irange.second; i++)
                for (int j=jrange.first; j<=jrange.second; j++)
                    for (int k=krange.first; k<=krange.second; k++)
                        indices.push_back(utuple(static_cast<uint>(i),static_cast<uint>(j),static_cast<uint>(k)));
        }
        else
            throw lm::Exception("Unsupported matrix dimension", matrixDims.size());

        delete [] sbuffer;
    }

    return indices;
}

pair<int,int> parseRange(char *rangeStr, int maxDim)
{
	if (strlen(rangeStr) == 1 && rangeStr[0] == ':')
		return pair<int,int>(0,maxDim-1);

	char * split = strstr(rangeStr, ":");
	if (split > rangeStr && split < rangeStr+strlen(rangeStr)-1)
	{
		char * start = rangeStr;
		*split = '\0';
		char * end = split+1;
		return pair<int,int>(atoi(start),atoi(end));
	}
	else
	{
		return pair<int,int>(atoi(rangeStr),atoi(rangeStr));
	}
}

vector<double> parseValues(string s)
{
    char * valueToks = new char[s.size()+1];
    strcpy(valueToks, s.c_str());
    int valueToksStart=0;

    // Remove any leading/trailing brackets.
    if (valueToks[0] == '[') valueToksStart++;
    if (valueToks[s.size()-1] == ']') valueToks[s.size()-1] = '\0';

    // Parse the values.
    vector<double> values;
    char * value = strtok(valueToks+valueToksStart, ";,");
    while (value != NULL)
    {
        double d = atof(value);
        values.push_back(d);
        value = strtok(NULL, ";,");
    }
    delete [] valueToks;
    return values;
}

vector<int32_t> parseIntValues(string s)
{
    char * valueToks = new char[s.size()+1];
    strcpy(valueToks, s.c_str());
    int valueToksStart=0;

    // Remove any leading/trailing brackets.
    if (valueToks[0] == '[') valueToksStart++;
    if (valueToks[s.size()-1] == ']') valueToks[s.size()-1] = '\0';

    // Parse the values.
    vector<int32_t> values;
    char * value = strtok(valueToks+valueToksStart, ";,");
    while (value != NULL)
    {
        int32_t i = atoi(value);
        values.push_back(i);
        value = strtok(NULL, ";,");
    }
    delete [] valueToks;
    return values;
}

vector<uint32_t> parseUintValues(string s)
{
    char * valueToks = new char[s.size()+1];
    strcpy(valueToks, s.c_str());
    int valueToksStart=0;

    // Remove any leading/trailing brackets.
    if (valueToks[0] == '[') valueToksStart++;
    if (valueToks[s.size()-1] == ']') valueToks[s.size()-1] = '\0';

    // Parse the values.
    vector<uint32_t> values;
    char * value = strtok(valueToks+valueToksStart, ";,");
    while (value != NULL)
    {
        uint32_t i = static_cast<uint32_t>(atoi(value));
        values.push_back(i);
        value = strtok(NULL, ";,");
    }
    delete [] valueToks;
    return values;
}
