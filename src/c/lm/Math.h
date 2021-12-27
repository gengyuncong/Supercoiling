/*
 * University of Illinois Open Source License
 * Copyright 2008-2012 Luthey-Schulten Group,
 * Copyright 2012-2016 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
 *
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
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
 * - Neither the names of the Luthey-Schulten Group, University of Illinois at
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
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
#ifndef LM_MATH_H_
#define LM_MATH_H_

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <numeric>

#include "lm/Types.h"

/*
 * constant definitions
 */
#define TWOPI       6.283185307179586476925287
#define PI          3.141592653589793238462643
#define PID2        1.570796326794896619231322
#define PID4        0.7853981634
#define PIOVER180   0.0174532925
#define NA          6.02214179e23
#define EPS         1e-9

/*
 * prefix definitions
 */
#define KIBI        1024
#define MEBI        1048576
#define GIBI        1073741824

/*
 * unary operations
 */
inline bool isPower2(unsigned int x)
{
    return !(x&(x - 1)) && x;
}

inline bool isPower2(unsigned long x)
{
    return !(x&(x - 1)) && x;
}

inline bool isPower2(unsigned long long x)
{
    return !(x&(x - 1)) && x;
}

inline unsigned int log2(unsigned int x)
{
    unsigned int r = 0;
    while (x>>=1)
        r++;
    return r;
}

inline unsigned int log2(unsigned long x)
{
    unsigned int r = 0;
    while (x>>=1)
        r++;
    return r;
}

inline unsigned int log2(unsigned long long x)
{
    unsigned int r = 0;
    while (x>>=1)
        r++;
    return r;
}

//#include <cmath>
//// simple rounding function
//// for interesting corner cases where this function won't work, see http://stackoverflow.com/a/4572677/425458
//double round(double d)
//{
//    return (d >= 0.0) ? floor(d + 0.5) : ceil(d - 0.5);
//}

/** binary operations
  */

/** simple almost equal function.
  * Tests equality based on an absolute tolerance value (EPS) used elsewhere in lm
  */
template <typename FloatType>
inline bool almostEqual(FloatType x0, FloatType x1)
{
    return std::abs(x0 - x1) < EPS;
}

/** tests whether one float is (nearly) evenly divisible by another
  */
template <typename FloatType>
inline bool almostDivisible(FloatType x, FloatType denominator)
{
    return fmod(x, denominator) < EPS;
}

template <typename IntType>
inline IntType ceilDiv(IntType val0, IntType val1)
{
    return (IntType)ceil(val0/(double)val1);
}

/** rounds x up to the nearest multiple of base
  */
template <typename FloatType>
inline FloatType roundNextMultiple(FloatType x, FloatType base)
{
    return (floor((x+EPS)/base) + 1)*base;
}

/** operations on single containers
  */

template <typename T, typename OutputIterator>
inline void cumprod(T* begin, T* end, OutputIterator outputIt)
{
    std::partial_sum(begin, end, outputIt, std::multiplies<T>());
}

template <typename InputIterator, typename OutputIterator>
inline void cumprod(InputIterator begin, InputIterator end, OutputIterator outputIt)
{
    std::partial_sum(begin, end, outputIt, std::multiplies<typename InputIterator::value_type>());
}

// product functor. ProductFunctor<T>::call(first, last) will return the product of an iterator range if T is a numeric type (returns 0 if the range is empty)
template <typename T, typename=typename EnableIf<IsNumeric<T>::value>::type>
struct ProductFunctor
{
    template <typename iterT> static inline T call(iterT first, iterT last)
    {
        if (first==last) return static_cast<T>(0);
        else             return std::accumulate(first, last, static_cast<T>(1), std::multiplies<T>());
    }
};

template <typename T>
inline T sum(T* begin, T* end)
{
    return std::accumulate(begin, end, static_cast<T>(0));
}

template <typename InputIterator>
inline typename InputIterator::value_type sum(InputIterator begin, InputIterator end)
{
    return std::accumulate(begin, end, static_cast<typename InputIterator::value_type>(0));
}

//template <typename T> struct _ProductFunctor<T, true> {template <typename iterT> static T call(iterT first, iterT last) {return std::accumulate(first, last, static_cast<T>(1), std::multiplies<T>());}};   //mul);}};
//template <typename T> struct _ProductFunctor<T, false> {template <typename iterT> static T call(iterT first, iterT last) {return *first;}};
//template <typename T> struct ProductFunctor {template <typename iterT> static T call(iterT first, iterT last) {return _ProductFunctor<T, IsNumeric<T>::value>::call(first, last);}};

#ifndef __cuda_cuda_h__
using std::min;
using std::max;
/*
template <class T> inline T min(T x, T y) {return x<y?x:y;}
template <class T> inline T min(T x, T y, T z) {return x<y?(x<z?x:z):(y<z?y:z);}
template <class T> inline T max(T x, T y) {return x>y?x:y;}
template <class T> inline T max(Tt x, T y, T z) {return x>y?(x>z?x:z):(y>z?y:z);}
*/
#endif

#endif
