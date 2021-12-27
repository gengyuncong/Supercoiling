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
 * Author(s): Max Klein, Elijah Roberts
 */

#ifndef LM_STATS_H_
#define LM_STATS_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>
#include <vector>

#include "lm/Types.h"

// constants used in erfinv (ie the inverse error function)
#define erfinv_a3 -0.140543331
#define erfinv_a2 0.914624893
#define erfinv_a1 -1.645349621
#define erfinv_a0 0.886226899

#define erfinv_b4 0.012229801
#define erfinv_b3 -0.329097515
#define erfinv_b2 1.442710462
#define erfinv_b1 -2.118377725
#define erfinv_b0 1

#define erfinv_c3 1.641345311
#define erfinv_c2 3.429567803
#define erfinv_c1 -1.62490649
#define erfinv_c0 -1.970840454

#define erfinv_d2 1.637067800
#define erfinv_d1 3.543889200
#define erfinv_d0 1

/** - Inverse error function. useful for calculating certain values related to the normal distribution.
  *     - Calculates a concrete value via a Taylor series approximation.
  *     - Code modified from libit, found at http://libit.sourceforge.net/math_8c-source.html.
  */
 inline double erfinv (double x)
{
    double x2, r, y;
    int  sign_x;

    if (x < -1.0 || x > 1.0) return std::numeric_limits<double>::quiet_NaN();

    if (x == 0.0) return 0;

    if (x > 0.0)
    {
        sign_x = 1;
    }
    else
    {
        sign_x = -1;
        x = -x;
    }

    if (x <= 0.7)
    {
        x2 = x * x;
        r = x * (((erfinv_a3 * x2 + erfinv_a2) * x2 + erfinv_a1) * x2 + erfinv_a0);
        r /= (((erfinv_b4 * x2 + erfinv_b3) * x2 + erfinv_b2) * x2 + erfinv_b1) * x2 + erfinv_b0;
    }
    else
    {
        y = sqrt (-log ((1 - x) / 2));
        r = (((erfinv_c3 * y + erfinv_c2) * y + erfinv_c1) * y + erfinv_c0);
        r /= ((erfinv_d2 * y + erfinv_d1) * y + erfinv_d0);
    }

    r = r * sign_x;
    x = x * sign_x;

    r -= (erf(r) - x) / (2 / sqrt(PI) * exp (-r * r));
    r -= (erf(r) - x) / (2 / sqrt(PI) * exp (-r * r));

    return r;
}

#undef erfinv_a3
#undef erfinv_a2
#undef erfinv_a1
#undef erfinv_a0

#undef erfinv_b4
#undef erfinv_b3
#undef erfinv_b2
#undef erfinv_b1
#undef erfinv_b0

#undef erfinv_c3
#undef erfinv_c2
#undef erfinv_c1
#undef erfinv_c0

#undef erfinv_d2
#undef erfinv_d1
#undef erfinv_d0

/** - The magic function for which normalZ(.95) == 1.96 (approximately).
  *     - Calculates how many standard deviations from the mean the cut-lines are for a given percentile (also centered on the mean) of the normal distribution.
  *     - Used in calculating confidence intervals.
  *     - Signature based on scipy.stats.norm.ppf, see http://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.norm.html for more details
  */
 inline double normalZ(double percentile, double mean=0.0, double std=1.0)
{
    return sqrt(2)*erfinv(percentile)*std + mean;
}

/** - Gives the index (of a sorted array) of the smallest member of the 100*<percentile> percentile.
  */
inline size_t percentileIndex(double percentile, size_t size)
{
    return static_cast<size_t>(ceil(percentile*size));
}

/** - classs that implements an iterative, numerically stable algorithm for calculating sample variance of a stream of samples, one at a time
  *     - has the advantage of not needing to keep every sample value in memory
  * - used as:
  *     StreamingVariance sv;
  *
  *     sv.push_back(17.0);
  *     sv.push_back(19.0);
  *     sv.push_back(24.0);
  *
  *     double variance = sv.var();
  *
  * - modified from https://www.johndcook.com/blog/standard_deviation/ and https://www.johndcook.com/blog/skewness_kurtosis/
  *     - lifted from from AOCP, Donal Knuth. Vol 2, page 232, 3rd edition
  */
class StreamingVariance
{
public:
    StreamingVariance(): _count(0),_sum(0.0),_mean(0.0),_rawVariance(0.0)
    {
    }

    void clear()
    {
        _count = 0;
        _sum = _mean = _rawVariance = 0.0;
    }

    // See Knuth TAOCP vol 2, 3rd edition, page 232 for complete description of algorithm
    void push_back(double x)
    {
        double oldMean;

        // update _count and _sum
        _count++;
        _sum += x;

        // store the old value of _mean, since we'll need it for updating _s
        oldMean = _mean;

        // update _mean
        _mean = oldMean + (x - oldMean)/_count;

        // update _rawVariance
        _rawVariance = _rawVariance + (x - oldMean)*(x - _mean);

//        double delta, delta_n, term1;
//
//        samples.push_back(x);
//        _sum += x;
//
//        long long n1 = _count;
//        _count++;
//        delta = x - _mean;
//        delta_n = delta / _count;
//        term1 = delta * delta_n * n1;
//        _mean += delta_n;
//        _rawVariance += term1;
    }

    long long count() const
    {
        return _count;
    }

    double sum() const
    {
        return _sum;
    }

    double mean() const
    {
        return _mean;
    }

    double variance() const
    {
        return _rawVariance/(_count - 1.0);
    }

    /*
     * - an estimate of the variance of the sample variance
     *     - relies on samples having been originally drawn from a normal distribution
     *     - see https://math.stackexchange.com/q/72975/13983 for details
     */
    double varianceOfVariance()
    {
        return pow(variance(), 2)*2/(_count - 1.0);
    }

    /*
     * - a more conservative estimate of the sample variance
     *     - relies on samples having been originally drawn from a normal distribution
     */
    double varianceUpperBound(double confidencePercentile=.95) const
    {
        return variance()*(1 + normalZ(confidencePercentile)*sqrt(2/(_count - 1.0)));
    }

    double standardDeviation() const
    {
        return sqrt(variance());
    }

    static StreamingVariance* combine(StreamingVariance* combined, const StreamingVariance& a, const StreamingVariance& b)
    {
        // update _rawVariance first, since _mean is used in the expression
        combined->_rawVariance = a._rawVariance + b._rawVariance + pow((b._mean - a._mean), 2) * a._count * b._count / (a._count + b._count);
        combined->_mean = (a._count*a._mean + b._count*b._mean) / (a._count + b._count);

        combined->_sum = a._sum + b._sum;
        combined->_count = a._count + b._count;

        return combined;

//        long long combinedCount = (a._count + b._count);
//        double delta = b._mean - a._mean;
//        double delta2 = delta*delta;
//        double delta3 = delta*delta2;
//        double delta4 = delta2*delta2;
//
//        combined->_mean = (a._count*a._mean + b._count*b._mean) / combinedCount;
//
//        combined->_rawVariance = a._rawVariance + b._rawVariance +
//                      delta2 * a._count * b._count / combinedCount;
//
//        combined->_sum = a._sum + b._sum;
//        combined->_count = combinedCount;
//        
//        return combined;
    }

    friend StreamingVariance operator+(const StreamingVariance& a, const StreamingVariance& b);
    StreamingVariance& operator+=(const StreamingVariance &rhs);

//    std::vector<double> samples;

private:
    long long _count;
    double _sum, _mean, _rawVariance;
};

inline StreamingVariance operator+(const StreamingVariance& a, const StreamingVariance& b)
{
    StreamingVariance combined;
    StreamingVariance::combine(&combined, a, b);

    return combined;
}

inline StreamingVariance& StreamingVariance::operator+=(const StreamingVariance& rhs)
{
    return *combine(this, *this, rhs);
}

class StreamingStats
{
public:
    StreamingStats(): n(0),M1(0.0),M2(0.0),M3(0.0),M4(0.0)
    {
    }

    void clear()
    {
        n = 0;
        M1 = M2 = M3 = M4 = 0.0;
    }

    void push_back(double x)
    {
        double delta, delta_n, delta_n2, term1;

        long long n1 = n;
        n++;
        delta = x - M1;
        delta_n = delta / n;
        delta_n2 = delta_n * delta_n;
        term1 = delta * delta_n * n1;
        M1 += delta_n;
        M4 += term1 * delta_n2 * (n*n - 3*n + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
        M3 += term1 * delta_n * (n - 2) - 3 * delta_n * M2;
        M2 += term1;
    }

    long long count() const
    {
        return n;
    }

    double mean() const
    {
        return M1;
    }

    double variance() const
    {
        return M2/(n-1.0);
    }

    double standardDeviation() const
    {
        return sqrt(variance());
    }

    double skewness() const
    {
        return sqrt(double(n)) * M3/ pow(M2, 1.5);
    }

    double kurtosis() const
    {
        return double(n)*M4 / (M2*M2) - 3.0;
    }

    friend StreamingStats operator+(const StreamingStats a, const StreamingStats b);
    StreamingStats& operator+=(const StreamingStats &rhs);

private:
    long long n;
    double M1, M2, M3, M4;
};

inline StreamingStats operator+(const StreamingStats a, const StreamingStats b)
{
    StreamingStats combined;

    combined.n = a.n + b.n;

    double delta = b.M1 - a.M1;
    double delta2 = delta*delta;
    double delta3 = delta*delta2;
    double delta4 = delta2*delta2;

    combined.M1 = (a.n*a.M1 + b.n*b.M1) / combined.n;

    combined.M2 = a.M2 + b.M2 +
                  delta2 * a.n * b.n / combined.n;

    combined.M3 = a.M3 + b.M3 +
                  delta3 * a.n * b.n * (a.n - b.n)/(combined.n*combined.n);
    combined.M3 += 3.0*delta * (a.n*b.M2 - b.n*a.M2) / combined.n;

    combined.M4 = a.M4 + b.M4 + delta4*a.n*b.n * (a.n*a.n - a.n*b.n + b.n*b.n) /
                                (combined.n*combined.n*combined.n);
    combined.M4 += 6.0*delta2 * (a.n*a.n*b.M2 + b.n*b.n*a.M2)/(combined.n*combined.n) +
                   4.0*delta*(a.n*b.M3 - b.n*a.M3) / combined.n;

    return combined;
}

inline StreamingStats& StreamingStats::operator+=(const StreamingStats& rhs)
{
    StreamingStats combined = *this + rhs;
    *this = combined;
    return *this;
}

#endif // LM_STATS_H_
