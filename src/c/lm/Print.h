/*
 * University of Illinois Open Source License
 * Copyright 2008-2010 Luthey-Schulten Group,
 * All rights reserved.
 * 
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
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
 * Urbana-Champaign, nor the names of its contributors may be used to endorse or
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
 * Author(s): Elijah Roberts, Max Klein
 */
#ifndef LM_PRINT_H_
#define LM_PRINT_H_

#include <cstdio>
#include <google/protobuf/message.h>
#include <google/protobuf/repeated_field.h>
#include <ostream>
#include <string>
#include <sstream>
#include <vector>

#include "lm/Traits.h"

namespace lm {

/**
 * Class for verbosity-configurable print function.
 */
class Print
{
public:
    static const int VERBOSE_DEBUG              = 10;
    static const int DEBUG                      =  9;
    static const int INFO                       =  4;
    static const int WARNING                    =  3;
    static const int ERROR                      =  2;
    static const int FATAL                      =  1;

    static std::string getDateTimeString();
    static void printDateTimeString();
    static void printf(int verbosity, const char * fmt, ...);
    static void printfStart(int verbosity, const char * fmt, ...);
    static void printfCont(int verbosity, const char * fmt, ...);
    static void printfEnd(int verbosity, const char * fmt, ...);
    static void printMsgDebug(int verbosity, const google::protobuf::Message& msg, size_t halfMaxSize=1048576);
    template <typename T> static const char* printf_format_string();
};

}

// forward declare printIterable to allow for printing of nested vectors (via recursive template resolution)
template <class T>
inline std::ostream& printIterable(std::ostream& stream, const T& iterable);

// put printNumeric in the top-level namespace
template <typename T> static void printNumeric(T num)
{
    std::printf(lm::Print::printf_format_string<T>(), num);
}

// EZ printing of std::pairs
template <class T0, class T1>
inline std::ostream& operator << (std::ostream& stream, const std::pair<T0, T1>& p)
{
    // no std::pair::const_iterator, so can't use printIterable(...)
    stream << "[";
    stream << p.first;
    stream << ", " << p.second;
    stream << "]";

    return stream;
}

// EZ printing of non-string iterables (anything with a .begin and .end)
template <class T> inline typename EnableIf<HasBeginEnd<T>::value and not IsSame<T, std::string>::value, std::ostream&>::type
operator << (std::ostream& stream, const T& iterable)
{
    return printIterable(stream, iterable);
}

// generic function for printing the contents of an iterable using a std::ostream
template <class T>
inline std::ostream& printIterable(std::ostream& stream, const T& iterable)
{
    typename T::const_iterator it = iterable.begin();

    stream << "[";
    if (it!=iterable.end())
    {
        stream << *it;
        it++;
    }
    for (;it!=iterable.end();++it)
    {
        stream << ", " << *it;
    }
    stream << "]";

    return stream;
}

#endif /* LM_PRINT_H_ */
