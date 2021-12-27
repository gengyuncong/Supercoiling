/*
 * Copyright 2016-2019 Johns Hopkins University
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
 * Author(s): Elijah Roberts
 */

#ifndef LM_FFPILOT_FFPILOTMATH_H_
#define LM_FFPILOT_FFPILOTMATH_H_

#include <algorithm>
#include <cmath>
#include <iterator>
#include <limits>
#include <vector>

#include "lm/Math.h"
#include "lm/Stats.h"
#include "lm/Types.h"
#include "lm/VectorMath.h"



// Calculates the variance of a Bernouli random variable based on its probability.
inline double bernouliVariance(const double probability)
{
    return probability * (1 - probability);
}

// Get the lower bound of a confidence interval for a Bernouli random variable based on its probability of number of trails.
// see https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval#Agresti-Coull_Interval for more details.
inline double bernouliConfidenceIntervalLowerBound(const double probability, const double trials, double confidence)
{
    double z = normalZ(confidence);
    double z2 = pow(z, 2);
    double n_twiddle = trials + z2;
    double p_twiddle = (trials*probability + 0.5*z2)/n_twiddle;

    // Return the lower bound of the confidence interval.
    return p_twiddle - z*sqrt((p_twiddle*(1 - p_twiddle)/n_twiddle));
}

// implementation of inverse of eq B1 from the "Automatic error control during forward flux sampling of rare events in master equation models" paper
inline uint64_t npilot(double errorGoal, double errorGoalConfidence, double probability=0)
{
    double nzsq = pow(normalZ(errorGoalConfidence), 2);

    return (uint64_t)((nzsq - probability*nzsq)/pow(errorGoal, 2));
}

// version of npilot that slightly overestimates when errorGoal >= npivot (but underestimates otherwise)
inline uint64_t _npilotConservative(double errorGoal, double errorGoalConfidence, double probability, double npivot)
{
    // round the z score up in order to get a nice round return value (some of the time, anyway)
    double nzsq = ceil(pow(normalZ(errorGoalConfidence), 2));

    return (uint64_t)((nzsq - probability*nzsq)/(npivot*errorGoal));
}

inline uint64_t npilotConservative(double errorGoal, double errorGoalConfidence, double probability=0, double npivot=.02)
{
    // return the max of our two approaches for calculating npilot
    return std::max(npilot(errorGoal, errorGoalConfidence), _npilotConservative(errorGoal, errorGoalConfidence, probability, npivot));
}

// exact formula for calculating the variance of a product of random variables. Based on their individual expected values and variances
inline double productVarianceExact(const std::vector<double>& expected, const std::vector<double>& variance)
{
    std::vector<double> expected2 = pow(expected, 2.0);

    return prod(variance + expected2) - prod(expected2);
}

#endif /* LM_FFLUX_MATH_H_ */
