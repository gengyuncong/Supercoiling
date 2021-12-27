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
 * Author(s): Max Klein
 */
#include "lm/gtest.h"

#include <stdio.h>
#include <unistd.h>
#include "gtest/gtest.h"
#include "gmock/gmock.h"

#include "lm/input/Input.h"
#include "lm/input/OrderParameters.pb.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

const char* filenamesLiteralOParams[] = {TESTDATA_ROOT "/oparamsInput.sfile"};
const vector<string> filenamesOParams(filenamesLiteralOParams, filenamesLiteralOParams+1);

const int speciesCountsLiteral[] = {0,1,2,3,4,5,6};
const vector<int> speciesCounts(speciesCountsLiteral, speciesCountsLiteral+7);

const double timesLiteral[] = {0.0};
const vector<double> times(timesLiteral, timesLiteral+1);

class OParamsFixture : public ::testing::Test
{
public:
    OParamsFixture(): input(filenamesOParams)
    {
//        printCWD();
    }

    lm::input::Input input;
};

TEST_F(OParamsFixture, readSFileInput_test)
{
    // make some references to the innards of the FFPilotSimulationInput msg for easy access
    const lm::oparam::OParams& oparams = input.getOrderParameters();

    double time = 0.0;
    const int state[] = {3, 8, 1, 4, 8, 0, 0};

    double absolute_tolerance = 1e-10;

    EXPECT_NEAR(-1, oparams.at(0)->calc(state, time), absolute_tolerance);
    EXPECT_NEAR(5, oparams.at(1)->calc(state, time), absolute_tolerance);
    EXPECT_NEAR(3, oparams.at(2)->calc(state, time), absolute_tolerance);
};
