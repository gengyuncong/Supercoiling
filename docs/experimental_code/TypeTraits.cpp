/*
 * University of Illinois Open Source License
 * Copyright 2012-2016 Roberts Group,
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
#include <iomanip>
#include <iostream>

#include <deque>
#include <list>
#include <vector>

#include "TypeTraits.h"

using std::cout;
using std::endl;
using std::deque;
using std::list;
using std::setw;
using std::vector;

#define INITDEQUE(type, vec, arr) deque<type> vec(arr, arr + sizeof(arr) / sizeof(arr[0]))
#define INITLIST(type, vec, arr) list<type> vec(arr, arr + sizeof(arr) / sizeof(arr[0]))
#define INITVEC(type, vec, arr) vector<type> vec(arr, arr + sizeof(arr) / sizeof(arr[0]))

double dubsArr[] = {46.855381265262764,85.9037082022656,69.58012107454921,62.65411941429959,8.381362232111467,71.74948490391098,59.834342857011144,68.58416960093038,28.022222667997177,51.02684642264451};

int main()
{
    // direct tests of the type traits themselves

    cout << setw(80) << "HasBegin<std::vector<double> >::value: " << HasBegin<vector<double> >::value << endl;
    cout << setw(80) << "HasBegin<double>::value: "               << HasBegin<double>::value          << endl;
    cout << setw(80) << "HasBegin<double*>::value: "              << HasBegin<double*>::value         << endl;

    cout << setw(80) << "HasBeginEnd<std::vector<double> >::value: " << HasBeginEnd<vector<double> >::value << endl;
    cout << setw(80) << "HasBeginEnd<double>::value: "               << HasBeginEnd<double>::value          << endl;
    cout << setw(80) << "HasBeginEnd<double*>::value: "              << HasBeginEnd<double*>::value         << endl;

    cout << setw(80) << "HasPushBack<std::vector<double> >::value: "                  << HasPushBack<vector<double> >::value                   << endl;
    cout << setw(80) << "HasPushBack<typename std::vector<double>::iteator>::value: " << HasPushBack<typename vector<double>::iterator>::value << endl;
    cout << setw(80) << "HasPushBack<double>::value: "                                << HasPushBack<double>::value                            << endl;
    cout << setw(80) << "HasPushBack<double*>::value: "                               << HasPushBack<double*>::value                           << endl;
    
    cout << endl;
    
    // test of functions that utilize the type traits
    INITDEQUE(double, dubsdeque, dubsArr);
    INITLIST(double, dubslist, dubsArr);
    INITVEC(double, dubsvec, dubsArr);

    cout << setw(80) << "enabled printing of std::deque<double> using HasBeginEnd: "  << dubsdeque << endl;
    cout << setw(80) << "enabled printing of std::list<double> using HasBeginEnd: "   << dubslist  << endl;
    cout << setw(80) << "enabled printing of std::vector<double> using HasBeginEnd: " << dubsvec   << endl;
}
