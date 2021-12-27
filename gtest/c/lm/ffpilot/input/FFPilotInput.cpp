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

#include "lm/ffpilot/input/FFPilotInput.h"
#include "lm/ffpilot/input/FFPilotStage.pb.h"

#include <string>
#include <vector>

using std::string;
using std::vector;

const char* filenamesLiteralFFPilotInputs[] = {TESTDATA_ROOT "/ffpilotSimulationInput.sfile"};
const vector<string> filenamesFFPilotInputs(filenamesLiteralFFPilotInputs, filenamesLiteralFFPilotInputs+1);

const int speciesCountsLiteral[] = {0,1,2,3,4,5,6};
const vector<int> speciesCounts(speciesCountsLiteral, speciesCountsLiteral+7);

const double timesLiteral[] = {0.0};
const vector<double> times(timesLiteral, timesLiteral+1);

class FFPilotInputFixture : public ::testing::Test
{
public:
    FFPilotInputFixture(): ffpilotInput(filenamesFFPilotInputs)
    {
//        printCWD();
    }

    lm::ffpilot::input::FFPilotInput ffpilotInput;
};

TEST_F(FFPilotInputFixture, readSFileInput_test)
{
    // make some references to the innards of the FFPilotSimulationInput msg for easy access
    const lm::ffpilot::input::FFPilotStage& ffpilotStage = ffpilotInput.getFFPilotStages().Get(0);
    const lm::ffpilot::input::FFPilotPhase& ffpilotPhase = ffpilotStage.ffpilot_phases(0);
    const lm::ffpilot::input::FFPilotPhaseLimit&  ffpilotPhaseLimit = ffpilotPhase.ffpilot_phase_limit();
    const lm::ffpilot::io::EndPoint& ffpilotPhaseStartPoint = ffpilotPhase.start_points(0);

    // test some scalar values in ffpilotStage
    EXPECT_EQ(1, ffpilotStage.tiling_id());
    EXPECT_EQ(1, ffpilotStage.basin_id());

    // test some scalar values in ffpilotPhase
    EXPECT_EQ(1, ffpilotPhase.tiling_id());
    EXPECT_EQ(1, ffpilotPhase.basin_id());
    EXPECT_EQ(4, ffpilotPhase.phase_id());

    // test some scalar values in ffpilotPhaseStartPoint
    EXPECT_EQ(1, ffpilotPhaseStartPoint.count());

    // test some vector/repeated values in ffpilotPhaseStartPoint
    int countSize = ffpilotPhaseStartPoint.species_coordinates_size();
    for (int i=0; i<countSize; i++)
    {
        int valExpected = speciesCounts[i];
        int valActual = ffpilotPhaseStartPoint.species_coordinates(i);

        EXPECT_EQ(valExpected, valActual);
    }

    int timesSize = ffpilotPhaseStartPoint.times_size();
    for (int i=0; i<timesSize; i++)
    {
        double valExpected = times[i];
        double valActual = ffpilotPhaseStartPoint.times(i);

        EXPECT_NEAR(valExpected, valActual, absolute_tolerance);
    }
}
