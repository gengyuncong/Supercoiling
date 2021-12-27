/*
 * Copyright 2016 Johns Hopkins University
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

#ifndef NDARRAYSERIALIZER_H
#define NDARRAYSERIALIZER_H

#include <cstring>
#include <zlib.h>

#ifdef OPT_SNAPPY
#include <snappy.h>
#endif

#include "lptf/Profile.h"
#include "lptf/ProfileCodes.h"
#include "robertslab/Exceptions.h"
#include "robertslab/Types.h"
#include "robertslab/pbuf/NDArray.pb.h"


namespace robertslab {
namespace pbuf {

template <typename T> robertslab::pbuf::NDArray_DataType NDArray_datatype_code();

// Configure the default compression algorithm.
#ifdef OPT_SNAPPY
#define DEFAULT_COMPRESSION robertslab::pbuf::NDArraySerializer::SNAPPY
#else
#define DEFAULT_COMPRESSION robertslab::pbuf::NDArraySerializer::DEFLATE
#endif

class NDArraySerializer
{
public:
    enum CompressionType {DEFLATE,
#ifdef OPT_SNAPPY
                          SNAPPY,
#endif
                          NONE};

public:
    template <typename T> static robertslab::pbuf::NDArray* serializeAllocate(const ndarray<T>& array, CompressionType compressionType=DEFAULT_COMPRESSION)
    {
        robertslab::pbuf::NDArray* msg = new robertslab::pbuf::NDArray();
        serializeInto(msg, array, compressionType);
        return msg;
    }

    template <typename T> static void serializeInto(robertslab::pbuf::NDArray* msg, const T* data, utuple shape, CompressionType compressionType=DEFAULT_COMPRESSION)
    {
        ndarray<T> array(shape, data, false);
        serializeInto(msg, array, compressionType);
    }

    template <typename T> static void serializeInto(robertslab::pbuf::NDArray* msg, const ndarray<T>& array, CompressionType compressionType=DEFAULT_COMPRESSION)
    {
        PROF_BEGIN(PROF_NDARRAY_SERIALIZE);

        // Clear the message.
        msg->Clear();

        // Set the array order.
        msg->set_array_order(convertArrayOrder(array.arrayOrder));

        // Set the data type.
        msg->set_data_type(NDArray_datatype_code<T>());

        // Set the shape.
        for (uint i=0; i<array.shape.len; i++)
            msg->add_shape(array.shape[i]);

        // See if we need to compress the data.
        if (compressionType == DEFLATE)
        {
            // Store deflated.
            msg->set_compressed_deflate(true);
            msg->set_compressed_snappy(false);
            std::string* data = msg->mutable_data();
            size_t dataSizeEstimate=compressBound(array.size*sizeof(T));
            data->resize(dataSizeEstimate);
            RL_ZLIB_EXCEPTION_CHECK(compress((unsigned char*)&((*data)[0]), &dataSizeEstimate, (const unsigned char*)array.values, array.size*sizeof(T)));
            data->resize(dataSizeEstimate);
        }
#ifdef OPT_SNAPPY
        else if (compressionType == SNAPPY)
        {
                // Store with snappy compression.
                msg->set_compressed_deflate(false);
                msg->set_compressed_snappy(true);
                std::string* data = msg->mutable_data();
                size_t dataSizeEstimate=snappy::MaxCompressedLength(array.size*sizeof(T));
                data->resize(dataSizeEstimate);
                snappy::RawCompress((const char*)array.values, array.size*sizeof(T), (char*)&((*data)[0]), &dataSizeEstimate);
                data->resize(dataSizeEstimate);
        }
#endif
        else
        {
            // Store without compression.
            msg->set_compressed_deflate(false);
            msg->set_compressed_snappy(false);
            std::string* data = msg->mutable_data();
            data->resize(array.size*sizeof(T));
            memcpy((unsigned char*)&((*data)[0]), (const unsigned char*)array.values, array.size*sizeof(T));
        }

        PROF_END(PROF_NDARRAY_SERIALIZE);
    }

    template <typename T> static void serializeIntoNoCopy(robertslab::pbuf::NDArray* msg, const T* data, utuple shape, CompressionType compressionType=DEFAULT_COMPRESSION)
    {
        serializeInto(msg, ndarray<T>(shape, data, false), compressionType);
    }

    template <typename T> static ndarray<T> deserialize(const robertslab::pbuf::NDArray& msg, size_t alignment=0)
    {
        // Get the shape of the ndarray.
        tuple<uint> shape(msg.shape().size(), (const uint*)msg.shape().data());

        // Create the ndarray.
        ndarray<T> array(shape, alignment, convertArrayOrder(msg.array_order()));

        // Deserialize the message.
        deserializeInto(&array, msg);

        return array;
    }

    template <typename T> static ndarray<T>* deserializeAllocate(const robertslab::pbuf::NDArray& msg, size_t alignment=0)
    {
        // Get the shape of the ndarray.
        tuple<uint> shape(msg.shape().size(), (const uint*)msg.shape().data());

        // Allocate the ndarray.
        ndarray<T>* array = new ndarray<T>(shape, alignment, convertArrayOrder(msg.array_order()));

        // Deserialize the message.
        deserializeInto(array, msg);

        return array;
    }

    template <typename T> static void deserializeInto(T* data, utuple shape, const robertslab::pbuf::NDArray& msg)
    {
        ndarray<T> array(shape, convertArrayOrder(msg.array_order()));
        deserializeInto(&array, msg);
        memcpy(data, array.values, array.size*sizeof(T));
    }

    template <typename T> static void deserializeInto(ndarray<T>* array, const robertslab::pbuf::NDArray& msg)
    {
        PROF_BEGIN(PROF_NDARRAY_DESERIALIZE);

        // Check that the datatype matches.
        if (msg.data_type() != NDArray_datatype_code<T>()) throw robertslab::InvalidArgException("In deserializeInto, type of NDArray message does not match that of destination.\n"
                                                                                                    "msg type: %s, destination type: %s", NDArray_DataType_Name(msg.data_type()).c_str(), NDArray_DataType_Name(NDArray_datatype_code<T>()).c_str());

        // Check that the shapes match.
        tuple<uint> shape(msg.shape().size(), (const uint*)msg.shape().data());
        if (shape != array->shape) throw robertslab::InvalidArgException("array", "the array must have the same size as the message", msg.data_type());

        // See if we need to decompress the data.
        if (msg.compressed_deflate())
        {
            size_t tmpBufferSize = array->size*sizeof(T);
            RL_ZLIB_EXCEPTION_CHECK(uncompress((unsigned char *)array->values, &tmpBufferSize, (unsigned char*)&(msg.data()[0]), msg.data().size()));
            if (tmpBufferSize != array->size*sizeof(T))
                throw robertslab::Exception("error during ndarray inflate deserialization, wrong number of bytes decompressed", tmpBufferSize, array->size*sizeof(T));

        }
        else if (msg.compressed_snappy())
        {
#ifdef OPT_SNAPPY
            size_t outputSize;
            if (!snappy::GetUncompressedLength((const char*)&(msg.data()[0]), msg.data().size(), &outputSize))
                throw robertslab::Exception("error during ndarray snappy deserialization, could not get uncompressed length");
            if (outputSize != array->size*sizeof(T))
                throw robertslab::Exception("error during ndarray snappy deserialization, wrong number of bytes to decompressed", outputSize, array->size*sizeof(T));
            if (!snappy::RawUncompress((const char*)&(msg.data()[0]), msg.data().size(), (char *)array->values))
                throw robertslab::Exception("error during ndarray snappy deserialization, could not get uncompressed data");
#else
            throw robertslab::InvalidArgException("msg", "support for snappy decompression is not available");
#endif
        }
        else
        {
            if (msg.data().size() != array->size*sizeof(double)) throw robertslab::InvalidArgException("msg", "inconsistent size during ndarray deserialization", msg.data().size(), array->size);
            memcpy(array->values, (const unsigned char*)&(msg.data()[0]), array->size*sizeof(T));
        }

        PROF_END(PROF_NDARRAY_DESERIALIZE);
    }

    static ndarray_ArrayOrder::Order convertArrayOrder(const robertslab::pbuf::NDArray::ArrayOrder arrayOrder)
    {
        switch (arrayOrder)
        {
        case robertslab::pbuf::NDArray::ROW_MAJOR:
            return ndarray_ArrayOrder::ROW_MAJOR;
        case robertslab::pbuf::NDArray::COLUMN_MAJOR:
            return ndarray_ArrayOrder::COLUMN_MAJOR;
        case robertslab::pbuf::NDArray::IMPL_ORDER:
            return ndarray_ArrayOrder::IMPL_ORDER;
        default:
            throw InvalidArgException("Invalid array order passed to convertArrayOrder. \n"
                                      "arrayOrder: %s", arrayOrder);
        }
    }

    static robertslab::pbuf::NDArray::ArrayOrder convertArrayOrder(ndarray_ArrayOrder::Order arrayOrder)
    {
        switch (arrayOrder)
        {
        case ndarray_ArrayOrder::ROW_MAJOR:
            return robertslab::pbuf::NDArray::ROW_MAJOR;
        case ndarray_ArrayOrder::COLUMN_MAJOR:
            return robertslab::pbuf::NDArray::COLUMN_MAJOR;
        case ndarray_ArrayOrder::IMPL_ORDER:
            return robertslab::pbuf::NDArray::IMPL_ORDER;
        default:
            throw InvalidArgException("Invalid array order passed to convertArrayOrder. \n"
                                      "arrayOrder: %s", arrayOrder);
        }
    }


};

}
}

#endif // NDARRAYSERIALIZER_H
