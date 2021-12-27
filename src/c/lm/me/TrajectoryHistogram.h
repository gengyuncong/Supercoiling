/*
 * Copyright 2019 Johns Hopkins University
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

#ifndef LM_ME_TRAJECTORY_HISTOGRAM_H
#define LM_ME_TRAJECTORY_HISTOGRAM_H

#include <limits>
#include <vector>

#include "lm/Exceptions.h"
#include "lm/Math.h"
#include "robertslab/pbuf/NDArray.pb.h"
#include "robertslab/pbuf/NDArraySerializer.h"

namespace lm {
namespace me {


template <typename T> class TrajectoryHistogram
{
public:
    TrajectoryHistogram(ndarray<int32_t> valueIndicesArg, ndarray<T> edges)
    :valueIndices(valueIndicesArg.shape),edges(edges),histogram(calculateHistogramShape(edges))
    {
        if (valueIndicesArg.shape.len != 1) THROW_EXCEPTION(lm::RuntimeException, "Invalid shape for value indices array: %d", valueIndices.shape.len);
        for (size_t i=0; i<valueIndicesArg.shape[0]; i++)
            this->valueIndices[i] = valueIndicesArg[i];
    }

    TrajectoryHistogram(ndarray<size_t> valueIndices, ndarray<T> edges)
    :valueIndices(valueIndices),edges(edges),histogram(calculateHistogramShape(edges))
    {
    }

    void clear()
    {
        histogram = 0.0;
    }

    void deserializeFrom(const robertslab::pbuf::NDArray& histogramMsg)
    {
        // Deserialize the values.
        robertslab::pbuf::NDArraySerializer::deserializeInto<double>(&histogram, histogramMsg);
    }

    void incrementCount(T* values, double count)
    {
    }

    void serializeInto(robertslab::pbuf::NDArray* histogramMsg)
    {
        // Serialize the values.
        robertslab::pbuf::NDArraySerializer::serializeInto<double>(histogramMsg, histogram);
    }

private:
    utuple calculateHistogramShape(const ndarray<T>& edges)
    {
        // The number of bins in each dimensions is the number of edges minus one.
        uint32_t shape[edges.shape.len];
        for (uint32_t i=0; i<edges.shape.len; i++)
        {
            shape[i] = edges.shape[i]-1;
        }
        return utuple(edges.shape.len, shape);
    }

private:
    ndarray<size_t> valueIndices;
    ndarray<T> edges;
    ndarray<double> histogram;
};

}
}

#endif // LM_ME_TRAJECTORY_HISTOGRAM_H
