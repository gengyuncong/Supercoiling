/*
 * Copyright 2016-2017 Johns Hopkins University
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

#ifndef ROBERTSLAB_TYPES_H
#define ROBERTSLAB_TYPES_H

#include <algorithm>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <google/protobuf/repeated_field.h>
#include <iostream>
#include <list>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <vector>

typedef unsigned int uint;

/*
 *  Array types.
 */

template <typename T> const char* printf_format_string();

template <typename T> struct tuple
{
    tuple()
    :len(0),values(NULL)
    {
    }

    tuple(const tuple& t)
    :len(t.len),values((len!=0)?(new T[len]):(NULL))
    {
        if (len != 0 && values != NULL && t.values != NULL)
            memcpy(values, t.values, sizeof(T)*len);
    }

    tuple(const T v1)
    :len(1),values(new T[len])
    {
        values[0] = v1;
    }

    tuple(const T v1, const T v2)
    :len(2),values(new T[len])
    {
        values[0] = v1;
        values[1] = v2;
    }

    tuple(const T v1, const T v2, const T v3)
    :len(3),values(new T[len])
    {
        values[0] = v1;
        values[1] = v2;
        values[2] = v3;
    }

    tuple(const T v1, const T v2, const T v3, const T v4)
    :len(4),values(new T[len])
    {
        values[0] = v1;
        values[1] = v2;
        values[2] = v3;
        values[3] = v4;
    }

    tuple(uint len, const T* valuesArray)
    :len(len),values(new T[len])
    {
        memcpy(values, valuesArray, sizeof(T)*len);
    }

    tuple(const std::list<T>& valuesList)
    :len(valuesList.size()),values(new T[len])
    {
        int i=0;
        for (typename std::list<T>::iterator it = valuesList.begin(); it != valuesList.end(); it++)
            values[i++] = *it;
    }

    tuple(const std::vector<T>& valuesVector)
    :len(valuesVector.size()),values(new T[len])
    {
        for (uint i=0; i<valuesVector.size(); i++)
            values[i] = valuesVector[i];
    }

    tuple(const google::protobuf::RepeatedField<T>& repFieldRef)
    :len(repFieldRef.size()),values(new T[repFieldRef.size()])
    {
        memcpy(values, repFieldRef.data(), sizeof(T)*len);
    }

    virtual ~tuple()
    {
        if (values != NULL) delete[] values; values = NULL;
    }

    tuple& operator=(const tuple<T>& t)
    {
        // See if the length of this tuple is different than the argument tuple.
        if (len != t.len)
        {
            // Free the old memory.
            if (values != NULL) delete[] values; values = NULL;

            // Save the new length.
            len = t.len;

            // Allocate the new memory.
            if (len != 0) values = new T[len];
        }

        // Copy the values.
        if (len != 0 && values != NULL && t.values != NULL)
            memcpy(values, t.values, sizeof(T)*len);

        return *this;
    }

    bool operator==(const tuple<T>& t) const
    {
        if (len != t.len)
        {
            return false;
        }
        for (uint i=0; i<len; i++)
        {
            if (values[i] != t.values[i])
            {
                return false;
            }
        }
        return true;
    }

    bool operator!=(const tuple<T>& t) const
    {
        return !(*this== t);
    }

    const T operator[](const uint index) const
    {
        return get(index);
    }

    T* begin() {return values;}
    T* end() {return values + len;}
    const T* begin() const {return values;}
    const T* end() const {return values + len;}

    const T get(const uint index) const
    {
        if (index < len) return values[index];
        else throw std::invalid_argument("index: index exceeded length of tuple");
    }

    // check if the tuple contains no more than one value i!=1 and no values i==0
    bool is1d() const {return max() == product();}

    T max() const {return std::max(begin(), end());}
    T product() const {return std::accumulate(begin(), end(), static_cast<T>(1), std::multiplies<T>());}
    T sum() const {return std::accumulate(begin(), end(), static_cast<T>(1), std::plus<T>());}

    // print contents to stdout
    void print(const char* suffix="") const
    {
        std::cout << toString(suffix);
    }

    // print contents to a string
    #define BUFSIZE 1024
    std::string toString(const char* suffix="") const
    {
        // allocate a buffer for C-style formatting of numbers
        char buffer[BUFSIZE + 1];
        std::stringstream reprStream("(");

        for (uint i=0; i<len; i++)
        {
            if (i > 0) reprStream << ',';
            snprintf(buffer, BUFSIZE, printf_format_string<T>(), values[i]);
            reprStream << buffer;
        }
        reprStream << ")" << suffix;

        return reprStream.str();
    }
    #undef BUFSIZE


public:
    uint len;
    T* values;
};

typedef tuple<uint> utuple;

struct ndarray_ArrayOrder
{
public:
	enum Order {
    		ROW_MAJOR             = 0,    // Last dimension contiguous.
    		COLUMN_MAJOR          = 1,    // First dimension contiguous
    		IMPL_ORDER            = 2     // Ordering specific to the implementation.
	};
};

template <typename T> struct ndarray
{
public:

public:
    // construct empty 0D ndarray
    ndarray()
    :arrayOrder(ndarray_ArrayOrder::ROW_MAJOR),shape(tuple<uint>()),size(0),alignment(0),allocatedValues(false),values(NULL)
    {
    }

    // construct empty 1D ndarray based on size
    ndarray(const uint size, size_t alignment=0, ndarray_ArrayOrder::Order arrayOrder=ndarray_ArrayOrder::ROW_MAJOR)
    :arrayOrder(arrayOrder),shape(tuple<uint>(size)),size(size),alignment(alignment),allocatedValues(false),values(allocateMemory(size,alignment))
    {
        memset(values, 0, sizeof(T)*size);
    }

    // construct empty ND ndarray based on shape
    ndarray(const tuple<uint>& shape, size_t alignment=0, ndarray_ArrayOrder::Order arrayOrder=ndarray_ArrayOrder::ROW_MAJOR)
    :arrayOrder(arrayOrder),shape(shape),size(calculateSize(shape)),alignment(alignment),allocatedValues(false),values(allocateMemory(size,alignment))
    {
        memset(values, 0, sizeof(T)*size);
    }

    // construct ndarray from pointer valuesArray, if copyValues == false this is done by taking ownership of the pointer
    ndarray(const tuple<uint>& shape, const T* valuesArray, bool copyValues=true, size_t alignment=0, ndarray_ArrayOrder::Order arrayOrder=ndarray_ArrayOrder::ROW_MAJOR)
    :arrayOrder(arrayOrder),shape(shape),size(calculateSize(shape)),alignment(copyValues?alignment:0),allocatedValues(false),values(copyValues?allocateMemory(size,alignment):const_cast<T*>(valuesArray))
    {
        if (copyValues) memcpy(values, valuesArray, sizeof(T)*size);
    }

    // copy constructor
    ndarray(const ndarray& a)
    :arrayOrder(a.arrayOrder),shape(a.shape),size(a.size),alignment(a.alignment),allocatedValues(false),values(allocateMemory(size,alignment))
    {
        if (values != NULL && a.values != NULL)
            memcpy(values, a.values, sizeof(T)*size);
    }

    virtual ~ndarray()
    {
        deallocateMemory();
    }

    void reshape(const tuple<uint>& newShape)
    {
        // If the shape is different, reinitialize the array.
        if (shape != newShape)
        {
            // Free the existing memory.
            deallocateMemory();

            // Store the new shape and size.
            shape = newShape;
            size = calculateSize(shape);
            values = allocateMemory(size,alignment);

            // Zero out the new memory.
            memset(values, 0, sizeof(T)*size);

            // TODO: copy as much of the old memory as possible.
        }
    }

    ndarray& operator=(const ndarray& a)
    {
        // If the source and destination are not the same size, free the current memory and allocate new memory.
        if (shape != a.shape || size != a.size)
        {
            // Free the memory.
            deallocateMemory();

            // Store the new shape and size.
            arrayOrder = a.arrayOrder;
            shape = a.shape;
            size = a.size;
            alignment = a.alignment;

            // Allocate new memory.
            values = allocateMemory(size,alignment);
        }

        // Copy the data.
        memcpy(values, a.values, sizeof(T)*size);

        return *this;
    }

    ndarray& operator=(const T& v)
    {
        for (uint i=0; i<size; i++)
            values[i] = v;
        return *this;
    }

    const T& operator[](const uint index) const
    {
        return (const_cast<ndarray *>(this))->get(utuple(index));
    }

    T& operator[](const uint index)
    {
        return get(utuple(index));
    }

    const T& operator[](const tuple<uint>& index) const
    {
        return (const_cast<ndarray *>(this))->get(index);
    }

    T& operator[](const tuple<uint>& index)
    {
        return get(index);
    }

    T* begin() {return values;}
    T* end() {return values + size;}
    const T* begin() const {return values;}
    const T* end() const {return values + size;}

    T& get(const uint index)
    {
        return get(utuple(index));
    }

    T& get(const tuple<uint>& index)
    {
        // Validate the index.
        if (index.len != shape.len) throw std::invalid_argument("index: index tuple must have the same length as the shape of an ndarray");
        for (uint i=0; i<shape.len; i++)
            if (index[i] >= shape[i]) throw std::invalid_argument("index: value of index exceeded ndarry length for the dimension");

        // Calculate the position.
        uint position=0;
        switch (arrayOrder)
        {
        case ndarray_ArrayOrder::ROW_MAJOR:
            for (uint i=0; i<shape.len; i++)
            {
                uint offset=1;
                for (uint j=i+1; j<shape.len; j++)
                    offset *= shape[j];
                position += index[i]*offset;
            }
            break;
        case ndarray_ArrayOrder::COLUMN_MAJOR:
            for (uint i=0; i<shape.len; i++)
            {
                uint offset=1;
                for (uint j=0; j<i; j++)
                    offset *= shape[j];
                position += index[i]*offset;
            }
            break;
        case ndarray_ArrayOrder::IMPL_ORDER:
            throw std::runtime_error("an ndarray with ordering IMPL_ORDER cannot be accessed by index");
        }

        // Return a reference to the element.
        return values[position];
    }

    // equivalent to numpy function of same name
    static uint ravelMultiIndex(const utuple& multiIndex, const utuple& shape)
    {
        uint position=0;
        for (uint i=0; i<shape.len; i++)
        {
            uint offset=1;
            for (uint j=i+1; j<shape.len; j++)
                offset *= shape[j];
            position += multiIndex[i]*offset;
        }

        return position;
    }

    uint ravelMultiIndex(const utuple& multiIndex)
    {
        return ravelMultiIndex(multiIndex, shape);
    }

    void set(const uint index, const T value)
    {
        set(utuple(index), value);
    }

    void set(const tuple<uint>& index, const T value)
    {
        get(index) = value;
    }

    template <typename Iter>
    void setRange(Iter beginRange, Iter endRange)
    {
        for (int i = 0; (beginRange!=endRange) and (i < size); ++i, ++beginRange)
        {
            values[i] = *beginRange;
        }
    }

    // equivalent to numpy function of same name
    static utuple unravelIndex(uint index, const utuple& shape)
    {
        std::vector<uint> multiIndex(shape.len, 0);
        std::vector<uint> minorShapes(shape.values + 1, shape.values + shape.len + 1);
        // for shape->(x, y, z), the partial_sum will store (y*z, z, 0) in multiIndex
        std::partial_sum (minorShapes.rbegin(), minorShapes.rend(), multiIndex.rbegin() + 1, std::multiplies<uint>());

        uint rem;
        for (uint i=0; i<shape.len - 1; i++)
        {
            rem = index % multiIndex[i];

            multiIndex[i] = index/multiIndex[i];
            index = rem;
        }
        multiIndex.back() = index;

        return utuple(multiIndex);
    }

    utuple unravelIndex(uint index)
    {
        return unravelIndex(index, shape);
    }

    ndarray& equalsDifference(const ndarray& a1, const ndarray& a2)
    {
        if (shape != a1.shape) std::invalid_argument("a1: all ndarrays during equals difference must be of the same shape");
        if (shape != a2.shape) std::invalid_argument("a2: all ndarrays during equals difference must be of the same shape");
        if (size != a1.size) std::invalid_argument("a1: all ndarrays during equals difference must have the same size");
        if (size != a2.size) std::invalid_argument("a2: all ndarrays during equals difference must have the same size");
        for (uint i=0; i<size; i++)
            values[i] = a1.values[i]-a2.values[i];
        return *this;
    }

    T max() const {return std::max(begin(), end());}

    T prod() const
    {
        T total = T(1);
        for (size_t i=0; i<size; i++)
            total *= values[i];
        return total;
    }

    T sum() const
    {
        T total = T(0);
        for (size_t i=0; i<size; i++)
            total += values[i];
        return total;
    }

    double mean() const
    {
        return static_cast<double>(sum())/static_cast<double>(size);
    }

    double variance() const
    {
        double mu = mean();
        double sumSquare = 0.0;
        for (size_t i=0; i<size; i++)
            sumSquare += (values[i]-mu)*(values[i]-mu);
        return static_cast<double>(sumSquare)/static_cast<double>(size-1);
    }

    uint rank() const {return shape.len;}

    // print contents to stdout
    void print(const char* suffix="") const
    {
        std::cout << toString(NULL, suffix);
    }

    // print contents to a string
    #define BUFSIZE 1024
    std::string toString(const char* format=NULL, const char* suffix="") const
    {
        // Get the format string to use.
        const char* formatCode = (format!=NULL)?format:printf_format_string<T>();

        // Allocate a buffer for C-style formatting of numbers
        char buffer[BUFSIZE + 1];
        std::stringstream debugSS;

        if (shape.len == 1 && arrayOrder != ndarray_ArrayOrder::IMPL_ORDER)
        {
            debugSS << "[";
            for (uint i=0; i<shape[0]; i++)
            {
                if (i > 0) debugSS << (",");
                snprintf(buffer, BUFSIZE, formatCode, (*this)[tuple<uint>(i)]);
                debugSS << buffer;
            }
            debugSS << "]" << suffix;
        }
        else if (shape.len == 2 && arrayOrder != ndarray_ArrayOrder::IMPL_ORDER)
        {
            debugSS <<  "[[";
            for (uint i=0; i<shape[0]; i++)
            {
                if (i > 0) debugSS <<  " [";
                for (uint j=0; j<shape[1]; j++)
                {
                    if (j > 0) debugSS <<  ",";
                    snprintf(buffer, BUFSIZE, formatCode, (*this)[tuple<uint>(i,j)]);
                    debugSS << buffer;
                }
                debugSS << "]\n";
            }
            debugSS << "]" << suffix;
        }
        else if (shape.len == 3 && arrayOrder != ndarray_ArrayOrder::IMPL_ORDER)
        {
            debugSS << "[[[";
            for (uint k=0; k<shape[2]; k++)
            {
                if (k > 0) debugSS << " [[";
                for (uint i=0; i<shape[0]; i++)
                {
                    if (i > 0) debugSS << "  [";
                    for (uint j=0; j<shape[1]; j++)
                    {
                        if (j > 0) printf (",");
                        snprintf(buffer, BUFSIZE, formatCode, (*this)[tuple<uint>(i,j,k)]);
                        debugSS << buffer;
                    }
                    debugSS << "]\n";
                }
                debugSS << " ]\n";
            }
            debugSS << "]" << suffix;
        }
        else
        {
            debugSS << "[ordering: " << arrayOrder <<", rank: " << rank() << ", shape: " << shape.toString() << ", entries: " << int(size) << "]" << suffix;
        }

        return debugSS.str();
    }
    #undef BUFSIZE

private:
    uint calculateSize(tuple<uint> s)
    {
        if (s.len == 0) return 0;

        uint r = 1U;
        for (uint i=0; i<s.len; i++)
            r *= s[i];
        return r;
    }

    T* allocateMemory(uint size, size_t alignment)
    {
        if (size == 0)
        {
            allocatedValues = false;
            return NULL;
        }

        T* tmp;
        if (alignment > 0)
        {
            int _posix_ret_=posix_memalign(reinterpret_cast<void**>(&tmp), alignment*sizeof(T), size*sizeof(T));
            if (_posix_ret_ != 0) throw std::invalid_argument("alignment: could not allocate aligned memory");
        }
        else
        {
            tmp = static_cast<T*>(malloc(size*sizeof(T)));
        }

        // Mark that we allocated these values.
        allocatedValues = true;

        return tmp;
    }

    void deallocateMemory()
    {
        if (allocatedValues && values != NULL) free(values);
        values = NULL;
        allocatedValues = false;
    }


public:
    ndarray_ArrayOrder::Order arrayOrder;
    tuple<uint> shape;
    size_t size;
    size_t alignment;
    bool allocatedValues;
    T* values;
};

#endif // ROBERTSLAB_TYPES_H
