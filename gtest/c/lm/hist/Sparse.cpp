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

#include <cstdio>
#include <map>
#include <vector>

#include "gtest/gtest.h"
#include "gmock/gmock.h"
#include "lm/hist/Sparse.h"
#include "robertslab/pbuf/Sparse.pb.h"

typedef int ObsElem;
typedef std::vector<ObsElem> ObsVec;
typedef std::vector<ObsVec> ObsVec2D;

typedef typename lm::hist::_DefaultHistMap<ObsElem>::type Counter;
typedef typename lm::hist::_DefaultHistMap<ObsVec>::type Counter2D;

typedef lm::hist::Sparse<ObsElem> Hist1D;
typedef lm::hist::Sparse<ObsVec> HistND;

template <class T, size_t I, size_t J>
ObsVec2D initObsVec2D(T (&arr)[I][J]) {
    ObsVec2D vec;
    ObsVec tmpVec;

    T* rowBegin;
    T* rowEnd = &arr[0][0];
    for (int i=0; i<I; i++) {
        rowBegin = rowEnd;
        rowEnd += J;
        tmpVec.assign(rowBegin , rowEnd);
        vec.push_back(tmpVec);}
    return vec;}

template <class Cnt, class T>
Cnt makeCount(const T& t) {
    Cnt cnt;
    for (typename T::const_iterator it = t.begin(); it != t.end(); it++) cnt[*it]++;
    return cnt;}

#define INITVEC(vec, arr) vec(arr, arr + sizeof(arr) / sizeof(arr[0]))
#define INITVEC2D(vec, arr) vec(initObsVec2D(arr))
#define INITCOUNTER(cnt, vec) cnt(makeCount<Counter>(vec))
#define INITCOUNTER2D(cnt, vec) cnt(makeCount<Counter2D>(vec))

int obsArr0[] = {1,0,1,3,2,2,0,4,3,1,1,2,4,4,1,1,1,3,2,2,0,4,0,4,1,4,4,2,3,3,2,2,0,2,0,3,2,1,4,2,0,0,3,3,4,4,2,0,0,0};
int obsArr1[] = {4,1,3,2,1,2,0,0,1,1,3,4,3,1,0,1,1,4,0,4,0,4,4,1,4,4,4,3,1,0,2,4,2,3,2,3,2,1,2,1,2,1,4,1,1,1,3,2,4,2};
int obsArr2[] = {2,0,3,1,1,3,1,4,3,3,2,3,1,4,0,0,3,4,0,3,4,4,2,2,4,3,1,1,0,4,0,2,1,4,2,0,0,1,1,2,2,3,4,3,1,0,1,3,1,4};
int obsArr3[] = {1,3,0,4,3,2,3,4,3,3,2,2,3,4,0,3,3,4,3,3,1,4,2,2,2,3,3,4,0,1,1,0,3,1,1,0,1,3,1,4,3,0,2,2,4,1,0,0,2,0};

int obsArr2D0[100][3] = {{0,2,0},{0,1,1},{0,2,1},{1,1,0},{1,2,1},{2,0,0},{0,0,1},{1,0,0},{1,1,1},{1,2,0},
                         {2,2,0},{0,1,0},{1,2,1},{1,1,0},{2,0,0},{0,2,2},{0,1,2},{1,1,0},{2,0,2},{1,2,1},
                         {0,0,1},{1,2,2},{0,0,1},{1,0,0},{0,2,1},{1,0,1},{2,1,2},{1,1,1},{0,2,0},{1,0,2},
                         {2,1,0},{2,1,1},{2,0,2},{0,1,2},{0,2,1},{0,1,2},{1,0,2},{2,2,0},{1,0,0},{2,0,0},
                         {1,1,1},{0,1,1},{1,2,2},{0,1,0},{1,1,2},{2,1,2},{0,2,0},{2,1,2},{0,0,1},{0,0,0},
                         {2,1,0},{0,2,0},{1,2,1},{1,0,0},{0,0,2},{2,0,1},{1,0,2},{2,2,1},{2,1,1},{1,2,2},
                         {2,2,1},{2,0,2},{0,0,1},{2,2,2},{0,0,2},{1,0,2},{2,2,0},{0,2,2},{2,1,0},{0,2,1},
                         {2,0,1},{1,1,2},{1,2,0},{2,1,2},{0,0,0},{0,2,0},{2,1,2},{1,2,0},{0,0,0},{0,0,0},
                         {2,1,0},{1,0,2},{0,2,1},{1,1,0},{0,1,2},{0,2,2},{0,1,1},{0,0,0},{2,1,1},{2,1,1},
                         {1,1,1},{0,2,2},{2,0,0},{2,1,1},{1,2,1},{0,1,0},{2,1,2},{0,2,1},{1,1,1},{0,0,1},};

class SparseFixture : public ::testing::Test
{
public:
    SparseFixture(): 
        INITVEC(obsVec0, obsArr0), INITVEC(obsVec1, obsArr1), INITVEC(obsVec2, obsArr2), INITVEC(obsVec3, obsArr3),
        INITVEC2D(obsVec2D0, obsArr2D0),
        INITCOUNTER(counter0, obsVec0), INITCOUNTER(counter1, obsVec1), INITCOUNTER(counter2, obsVec2), INITCOUNTER(counter3, obsVec3),
        INITCOUNTER2D(counter2D0, obsVec2D0)
    {}

    template <class Iter, class Cnt, class Hist>
    static void const _check(const Cnt& cnt, const Hist& hist)
    {
        for (std::pair<Iter, Iter> its(Iter(cnt).begin(), Iter(hist->map).begin());
             its.first != Iter(cnt).end(); its.first++, its.second++)
        {
            EXPECT_EQ(*its.first, *its.second);
        }
    }

    template <class Iter, class Cnt, class T>
    static void const _check(const Cnt& cnt, const T* t)
    {
        int i = 0;
        for (Iter it = Iter(cnt).begin(); it != Iter(cnt).end(); it++, i++)
        {
            EXPECT_EQ(*it, t[i]);
        }
    }

    template <class Cnt, class Obs, class Count>
    static void const _check(const Cnt& cnt, const Obs* obsArr, const Count* countArr)
    {
        _check<typename KeyFlatConstIterPolicy<Cnt>::type>(cnt, obsArr);
        _check<ValConstIter<Cnt> >(cnt, countArr);
    }

    template <class Cnt, class Hist>
    static void const check(const Cnt& cnt, const Hist& hist)
    {
        // check the shape of the histogram
        EXPECT_EQ(cnt.size(), hist->shape()[0]);

        // check the values of the observations
        _check<typename RemovePointer<Hist>::type::obs_const_iterator>(cnt, hist);

        // check the values of the counts
        _check<typename RemovePointer<Hist>::type::count_const_iterator>(cnt, hist);
    }

    template <class Cnt>
    static void const check(const Cnt& cnt, const robertslab::pbuf::Sparse& msg)
    {
        ndarray<ObsElem>* obsArr = robertslab::pbuf::NDArraySerializer::deserializeAllocate<ObsElem>(msg.obs());
        ndarray<uint64_t>* countArr = robertslab::pbuf::NDArraySerializer::deserializeAllocate<uint64_t>(msg.counts());

        // check the shape of the deserialized sparse histogram ndarray messages
        EXPECT_EQ(cnt.size(), obsArr->shape[0]);
        EXPECT_EQ(cnt.size(), countArr->shape[0]);

        // check the values of the observations and the counts
        _check(cnt, obsArr->values, countArr->values);

        delete obsArr; delete countArr;
    }

    template <class Obs, class Cnt, class Hist>
    void _addObs_test(const Obs& obsVec, const Cnt& cnt, Hist* hist)
    {
        hist->clear();
        for (typename Obs::const_iterator it = obsVec.begin(); it != obsVec.end(); it++) hist->addObs(*it);
        check(cnt, hist);
    }

    template <class Obs, class Cnt, class Hist>
    void _addObses_test(const Obs& obsVec, const Cnt& cnt, Hist* hist)
    {
        hist->clear();
        hist->addObses(obsVec);
        check(cnt, hist);
    }

    template <class Obs, class Cnt, class Hist>
    void _serializeTo_test(const Obs& obsVec, const Cnt& cnt, Hist* hist)
    {
        hist->clear();
        msg.Clear();

        hist->addObses(obsVec);
        hist->serializeTo(&msg);

        check(cnt, msg);
    }

    template <class Obs, class Cnt, class Hist>
    void _serializeTo_deserializeFrom_test(const Obs& obsVec, const Cnt& cnt, Hist* hist)
    {
        hist->clear();
        msg.Clear();

        hist->addObses(obsVec);
        hist->serializeTo(&msg);
        hist->deserializeFrom(msg);

        check(cnt, hist);
    }

public:
    Hist1D sparse1D;
    HistND sparseND;

    ObsVec obsVec0,obsVec1,obsVec2,obsVec3;
    ObsVec2D obsVec2D0;

    Counter counter0,counter1,counter2,counter3;
    Counter2D counter2D0;

    robertslab::pbuf::Sparse msg;
};

TEST_F(SparseFixture, addObs_1d_test)
{
    _addObs_test(obsVec0, counter0, &sparse1D);
    _addObs_test(obsVec1, counter1, &sparse1D);
    _addObs_test(obsVec2, counter2, &sparse1D);
    _addObs_test(obsVec3, counter3, &sparse1D);
}

TEST_F(SparseFixture, addObses_1d_test)
{
    _addObses_test(obsVec0, counter0, &sparse1D);
    _addObses_test(obsVec1, counter1, &sparse1D);
    _addObses_test(obsVec2, counter2, &sparse1D);
    _addObses_test(obsVec3, counter3, &sparse1D);
}

TEST_F(SparseFixture, serializeTo_1d_test)
{
    _serializeTo_test(obsVec0, counter0, &sparse1D);
    _serializeTo_test(obsVec1, counter1, &sparse1D);
    _serializeTo_test(obsVec2, counter2, &sparse1D);
    _serializeTo_test(obsVec3, counter3, &sparse1D);
}

TEST_F(SparseFixture, serializeTo_deserializeFrom_1d_test)
{
    _serializeTo_deserializeFrom_test(obsVec0, counter0, &sparse1D);
    _serializeTo_deserializeFrom_test(obsVec1, counter1, &sparse1D);
    _serializeTo_deserializeFrom_test(obsVec2, counter2, &sparse1D);
    _serializeTo_deserializeFrom_test(obsVec3, counter3, &sparse1D);
}

TEST_F(SparseFixture, addObs_nd_test)
{
    _addObs_test(obsVec2D0, counter2D0, &sparseND);
}

TEST_F(SparseFixture, addObses_nd_test)
{
    _addObses_test(obsVec2D0, counter2D0, &sparseND);
}

TEST_F(SparseFixture, serializeTo_nd_test)
{
    _serializeTo_test(obsVec2D0, counter2D0, &sparseND);
}

TEST_F(SparseFixture, serializeTo_deserializeFrom_nd_test)
{
    _serializeTo_deserializeFrom_test(obsVec2D0, counter2D0, &sparseND);
}
