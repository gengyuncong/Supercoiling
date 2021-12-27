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

#ifndef LM_ME_TIMESERIESLIST_H
#define LM_ME_TIMESERIESLIST_H

#include <limits>
#include <vector>

#include "lm/Exceptions.h"
#include "lm/Math.h"
#include "robertslab/pbuf/NDArray.pb.h"
#include "robertslab/pbuf/NDArraySerializer.h"

namespace lm {
namespace me {


template <typename T> class TimeSeriesList
{
public:
    TimeSeriesList()
    {
        clear();
    }

    void clear()
    {
        numberColumns = 0;
        nextUpdateTime = std::numeric_limits<double>::infinity();
        updateInterval = std::numeric_limits<double>::infinity();
        values.clear();
        times.clear();
    }

    size_t size()
    {
        return times.size();
    }

    void initialize(size_t columns, double interval, double time)
    {
        numberColumns = columns;
        updateInterval = interval;
        nextUpdateTime = roundNextMultiple(time, interval);
    }

    void initialize(size_t columns, double interval, double time, T value)
    {
        numberColumns = columns;
        updateInterval = interval;
        nextUpdateTime = roundNextMultiple(time, interval);

        if (1 != numberColumns) THROW_EXCEPTION(lm::RuntimeException, "Invalid number of values added to time series: %d %d", 1, numberColumns);
        values.push_back(value);
        times.push_back(time);
    }

    void initialize(size_t columns, double interval, double time, T* valuesArray, size_t numberValues)
    {
        numberColumns = columns;
        updateInterval = interval;
        nextUpdateTime = roundNextMultiple(time, interval);

        if (numberValues != numberColumns) THROW_EXCEPTION(lm::RuntimeException, "Invalid number of values added to time series: %d %d", numberValues, numberColumns);
        for (size_t i=0; i<numberValues; i++)
            values.push_back(valuesArray[i]);
        times.push_back(time);
    }

    void append(T value, double time)
    {
        if (1 != numberColumns) THROW_EXCEPTION(lm::RuntimeException, "Invalid number of values added to time series: %d %d", 1, numberColumns);

        while (nextUpdateTime <= (time+EPS))
        {
            values.push_back(value);
            times.push_back(nextUpdateTime);
            nextUpdateTime += updateInterval;
        }
    }

    void append(T* valuesArray, size_t numberValues, double time)
    {
        if (numberValues != numberColumns) THROW_EXCEPTION(lm::RuntimeException, "Invalid number of values added to time series: %d %d", numberValues, numberColumns);

        while (nextUpdateTime <= (time+EPS))
        {
            for (size_t i=0; i<numberValues; i++)
                values.push_back(valuesArray[i]);
            times.push_back(nextUpdateTime);
            nextUpdateTime += updateInterval;
        }
    }

    void appendFinal(T value, double time)
    {
        if (1 != numberColumns) THROW_EXCEPTION(lm::RuntimeException, "Invalid number of values added to time series: %d %d", 1, numberColumns);

        // Write out the last species count, or ensure that it already has been
        if ((times.empty() || times.back() + EPS < time))
        {
            values.push_back(value);
            times.push_back(nextUpdateTime);
        }
    }

    void appendFinal(T* valuesArray, size_t numberValues, double time)
    {
        if (numberValues != numberColumns) THROW_EXCEPTION(lm::RuntimeException, "Invalid number of values added to time series: %d %d", numberValues, numberColumns);

        // Write out the last time, or ensure that it already has been
        if ((times.empty() || times.back() + EPS < time))
        {
            for (size_t i=0; i<numberValues; i++)
                values.push_back(valuesArray[i]);
            times.push_back(time);
        }
    }

    void serializeInto(robertslab::pbuf::NDArray* valuesMsg, robertslab::pbuf::NDArray* timesMsg)
    {
        if (times.size()*numberColumns != values.size()) THROW_EXCEPTION(lm::RuntimeException, "Invalid number of times and values in time series: %d %d %d", times.size(), values.size(), numberColumns);

        // Serialize the values.
        robertslab::pbuf::NDArraySerializer::serializeInto<T>(valuesMsg, values.data(), utuple(uint(times.size()),numberColumns));

        // Serialize the times.
        robertslab::pbuf::NDArraySerializer::serializeInto<double>(timesMsg, times.data(), utuple(uint(times.size())));
    }

private:
    size_t numberColumns;
    double nextUpdateTime;
    double updateInterval;
    std::vector<T> values;
    std::vector<double> times;
};

}
}

#endif // LM_ME_FPTDEQUE_H
