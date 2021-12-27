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
 * Author(s): Elijah Roberts, Max Klein
 */
#ifndef LM_PROTWRAP_FFLUXPHASEOUTPUT_H_
#define LM_PROTWRAP_FFLUXPHASEOUTPUT_H_

#include <algorithm>
#include <limits>
#include <map>
#include <vector>

#include "lm/ffpilot/FFPilotPhaseZeroTrajectory.h"
#include "lm/ffpilot/io/FFPilotPhaseOutput.pb.h"
#include "lm/limit/LimitCheckFunctions.h"
#include "lm/limit/LimitTrackingWrap.h"
#include "lm/io/LimitTracking.pb.h"
#include "lm/protowrap/NDArray.h"
#include "lm/protowrap/RepeatedMap.h"
#include "lm/rng/RandomGenerator.h"
#include "lm/rng/XORShift.h"
#include "lm/Stats.h"
#include "lm/trajectory/Trajectory.h"
#include "lm/Types.h"

#ifdef OPT_CUDA
#include "lm/rng/XORWow.h"
#endif

namespace lm {
namespace protowrap {

typedef lm::ffpilot::io::EndPoint EndPointMsg;
typedef lm::ffpilot::io::FFPilotPhaseOutput FFPilotPhaseOutputMsg;
typedef lm::ffpilot::io::StartPoint StartPointMsg;

typedef std::vector<int32_t> PointKey;

// map key getter and setter for the EndPointMap typedef
template <typename PointMsg>
PointKey getPointKey(const PointMsg& pointMsgConst)
{
    return PointKey(pointMsgConst.species_coordinates().begin(), pointMsgConst.species_coordinates().end());
}
template <typename PointMsg>
void setPointKey(PointMsg* pointMsg, const PointKey& pointKey)
{
    pointMsg->mutable_species_coordinates()->Clear();
    for (PointKey::const_iterator it=pointKey.begin();it!=pointKey.end();it++)
    {
        pointMsg->add_species_coordinates(*it);
    }
}
typedef RepeatedMap<EndPointMsg, PointKey, &getPointKey, &setPointKey> EndPointMap;
// each entry in an EndPointVector is a Pair of a pointer to an endpoint and an index. The index allows you to lookup the time
typedef PairVector<lm::ffpilot::io::EndPoint*, int> EndPointVector;

typedef RepeatedMap<StartPointMsg, PointKey, &getPointKey, &setPointKey> StartPointMap;

//class EndPoint
//{
//    EndPointMsg* wrappedMsgPtr;
//};
//
//class StartPoint
//{
//public:
//    StartPointMsg* wrappedMsgPtr;
//    EndPointMap successfulEndPointMap;
//    EndPointMap failedEndPointMap;
//};
//typedef std::map<PointKey, StartPoint> StartPointMap;

class FFPilotPhaseOutputWrap
{
public:
    typedef FFPilotPhaseOutputMsg Msg;

    FFPilotPhaseOutputWrap(size_t randomCacheSize=10*KIBI)
    :msgPtr(NULL),rng(NULL),randomDoublesStart(NULL),randomDoubles(NULL),randomDoublesEnd(NULL),randomIndexesStart(NULL),
     randomIndexes(NULL),randomIndexesEnd(NULL),randomCacheSize(randomCacheSize),randomIndexesDirty(true)
    {
    }

    FFPilotPhaseOutputWrap(Msg* msgMutablePtr, size_t randomCacheSize=10*KIBI)
    :msgPtr(NULL),rng(NULL),randomDoublesStart(NULL),randomDoubles(NULL),randomDoublesEnd(NULL),randomIndexesStart(NULL),
     randomIndexes(NULL),randomIndexesEnd(NULL),randomCacheSize(randomCacheSize),randomIndexesDirty(true)
    {
        setWrappedMsg(msgMutablePtr);
    }

    ~FFPilotPhaseOutputWrap() {destructRng(); destructRandomIndexes();}

// mutators
    void addEndPointPhaseZero(const lm::io::TrajectoryState& trajectoryState, lm::trajectory::Trajectory* trajectory, int burnInCount)
    {
        // TODO: fix burn in count
        // TODO: include consistency check constraining (phaseZeroSamples > burnInCount) somewhere

        // cast the Trajectory reference to a PhaseZeroTrajectory reference (its true type)
        lm::ffpilot::FFPilotPhaseZeroTrajectory* phaseZeroTrajectory = static_cast<lm::ffpilot::FFPilotPhaseZeroTrajectory*>(trajectory);

        // set wrapper on the limit_trackings field
        limitTrackingsWrap.setWrappedField(trajectoryState.limit_tracking_list().limit_trackings());

        // consistency checks
        if (limitTrackingsWrap.size()!=3)
        {
            throw ConsistencyException("In FFPilotPhaseOutputWrap::addEndPointPhaseZero, finished Forward Flux phase zero trajectories should have 3 tracked limits in their outputs; trajectory id %llu has %d", trajectoryState.trajectory_id(), limitTrackingsWrap.size());
        }
        for (int i=0;i<3;i++)
        {
            if (limitTrackingsWrap.Get(i).limit_id()!=i)
            {
                throw ConsistencyException("In FFPilotPhaseOutputWrap::addEndPointPhaseZero, finished Forward Flux phase zero trajectories should have 3 tracked limits in their outputs with limit_ids {0, 1, 2}; trajectory id %llu has limit tracking index %d with limit_id %d", trajectoryState.trajectory_id(), i, limitTrackingsWrap.Get(i).limit_id());
            }
        }

        // get the start and end times for the entire work unit part
        double workUnitStartTime = phaseZeroTrajectory->getLastTime();
        double workUnitEndTime = trajectoryState.cme_state().species_counts().time(trajectoryState.cme_state().species_counts().time_size() - 1);

        // fetch forth some data from the limit trackings (ie forward flux, basin entry, and basin exit)
        speciesCountWrap.setWrappedMsg(limitTrackingsWrap.Get(0).species_counts());
        timeWrapForwardFlux.setWrappedMsg(limitTrackingsWrap.Get(0).times());

        // TODO: figure out what to do with burnInCount/burnInTime
        //double burnInTime = 0.0;
        // update the count of flux events only if
        if (speciesCountWrap.size() > 0)
        {
            int32_t* speciesCountDataForwardFlux = speciesCountWrap.get_data();
            double* timeDataForwardFlux = timeWrapForwardFlux.get_data();

            // get some metadata about the species counts
            uint rows = speciesCountWrap.shape(0);
            uint columns = speciesCountWrap.shape(1);

            // check that the shapes of times and species count ndarrays are consistent
            if (rows!=timeWrapForwardFlux.size())
            {
                throw ConsistencyException("In FFPilotPhaseOutputWrap::addEndPointPhaseZero, the number of rows in speciesCountDataForwardFlux: %d did not equal the number of rows in timeDataForwardFlux: %d.", rows, timeWrapForwardFlux.size());
            }

            // TODO: figure out what to do with burnInCount/burnInTime
            // if burnInCount is greater than zero, figure out the time until the first non-burned flux event
            // double burnInTime = burnInCount > 0 ? timeDataForwardFlux[burnInCount - 1] : 0.0;

            // load points from forward flux events into successful endpoints
            for (int i=0;i<rows;i++)
            {
                // get the endpoint key from the final species counts of the trajectory
                pointKey.assign(speciesCountDataForwardFlux + i*columns, speciesCountDataForwardFlux + (i + 1)*columns);
                EndPointMsg* endPointMsg = successfulEndPointMap[pointKey];
                // increment the counter of observations of this endpoint
                endPointMsg->set_count(endPointMsg->count() + 1);
                endPointMsg->add_times(timeDataForwardFlux[i]);

                // TODO: decide if the creation of endPointVector should be done one at a time (as below) or all at once
                endPointVector.push_back(std::make_pair(endPointMsg, endPointMsg->count() - 1));
                randomIndexesDirty = true;
            }

            if (speciesCountWrap.compressed_deflate()) delete[] speciesCountDataForwardFlux;
            if (timeWrapForwardFlux.compressed_deflate()) delete[] timeDataForwardFlux;

            // add to the summary metrics
            msgPtr->set_successful_trajectories_launched_count(msgPtr->successful_trajectories_launched_count() + rows);
        }

        // combine the latest streaming samples to get the variance
        phaseWeightSV += phaseZeroTrajectory->phaseWeightSV;

        // clear the streamingVariance to ensure against any double counting
        phaseZeroTrajectory->phaseWeightSV.clear();

        // correct workUnitEndTime for burn in and for time spent outside of the region of the starting basin (see Valeriani 2007, Dinner 2010)
        msgPtr->set_successful_trajectories_launched_total_time(phaseWeightSV.sum());   //msgPtr->successful_trajectories_launched_total_time() + (workUnitEndTime - workUnitStartTime));
        // total uncorrected work unit time
        msgPtr->set_failed_trajectories_launched_total_time(msgPtr->failed_trajectories_launched_total_time() + (workUnitEndTime - workUnitStartTime)); //msgPtr->failed_trajectories_launched_total_time() + phaseZeroTrajectory->timeInOtherBasinsLast);
        //msgPtr->set_failed_trajectories_launched_total_time(phaseZeroTrajectory->timeInOtherBasins);

        // TODO: replace this upper bound estimation (which is based on an assumption of sample normality, which is not true in this case) with a better one based on resampling
        msgPtr->set_variance(phaseWeightSV.varianceUpperBound(.99));
    }

    void addEndPoint(const lm::io::TrajectoryState& trajectoryState, const lm::trajectory::Trajectory& trajectory)
    {
        // set wrapper on the limit_trackings field
        limitTrackingsWrap.setWrappedField(trajectoryState.limit_tracking_list().limit_trackings());

        // consistency checks
        if (limitTrackingsWrap.size()!=2) throw ConsistencyException("Finished Forward Flux phase n>0 trajectories should have 2 tracked limits in their outputs; trajectory id %llu has %d", trajectoryState.trajectory_id(), limitTrackingsWrap.size());
        for (int i=0;i<2;i++) {if (limitTrackingsWrap.Get(i).limit_id()!=i) throw ConsistencyException("Finished Forward Flux phase n>0 trajectories should have 2 tracked limits in their outputs with limit_ids {0, 1}; trajectory id %llu has limit tracking index %d with limit_id %d", trajectoryState.trajectory_id(), i, limitTrackingsWrap.Get(i).limit_id());}

        // fetch forth some time data from limit 0 (ie backward flux) and limit 1 (ie forward flux) tracking
        timeWrapBackwardFlux.setWrappedMsg(limitTrackingsWrap.Get(0).times());
        double* timeDataBackwardFlux = timeWrapBackwardFlux.get_data(true);
        timeWrapForwardFlux.setWrappedMsg(limitTrackingsWrap.Get(1).times());
        double* timeDataForwardFlux = timeWrapForwardFlux.get_data(true);

        // check if this trajectory fluxed backwards or forwards (and make sure it didn't somehow do both)
        if (timeWrapBackwardFlux.size()==1 and timeWrapForwardFlux.size()==0)       // branch for "failed" trajectories (ie ones that fluxed backward)
        {
            msgPtr->set_failed_trajectories_launched_count(msgPtr->failed_trajectories_launched_count() + 1);
            msgPtr->set_failed_trajectories_launched_total_time(msgPtr->failed_trajectories_launched_total_time() + timeDataBackwardFlux[0] - trajectory.getLastTime());

            // TODO: decide whether startpoint tracking is any good/helpful here
            speciesCountWrap.setWrappedMsg(limitTrackingsWrap.Get(0).species_counts());
            int32_t* speciesCountData = speciesCountWrap.get_data(true);
            uint columns = speciesCountWrap.shape(1);

            processFailedStartPoint(speciesCountData, speciesCountData + columns, timeDataBackwardFlux[0], trajectory);

            phaseWeightSV.push_back(0);
        }
        else if (timeWrapBackwardFlux.size()==0 and timeWrapForwardFlux.size()==1)  // branch for "successful" trajectories (ie ones that fluxed forward)
        {
            msgPtr->set_successful_trajectories_launched_count(msgPtr->successful_trajectories_launched_count() + 1);
            msgPtr->set_successful_trajectories_launched_total_time(msgPtr->successful_trajectories_launched_total_time() + timeDataForwardFlux[0] - trajectory.getLastTime());

            // since this is data from a "successful" trajectory (ie one that fluxed forward), add its endpoint to the list used to initialize the next phase
            speciesCountWrap.setWrappedMsg(limitTrackingsWrap.Get(1).species_counts());
            int32_t* speciesCountData = speciesCountWrap.get_data(true);
            uint columns = speciesCountWrap.shape(1);

            EndPointMsg* endPointMsg = processSuccessfulEndPoint(speciesCountData, speciesCountData + columns, timeDataForwardFlux[0]);

//            pointKey.assign(speciesCountData, speciesCountData + columns);
//            EndPointMsg* endPointMsg = successfulEndPointMap[pointKey];
//            endPointMsg->set_count(endPointMsg->count() + 1);
//            endPointMsg->add_times(timeDataForwardFlux[0]);

            // TODO: decide whether startpoint tracking is any good/helpful here
            processSuccessfulStartPoint(speciesCountData, speciesCountData + columns, timeDataForwardFlux[0], trajectory);

            // TODO: decide if the creation of endPointVector should be done one at a time (as below) or all at once
            endPointVector.push_back(std::make_pair(endPointMsg, endPointMsg->count() - 1));
            randomIndexesDirty = true;

            phaseWeightSV.push_back(1);

            if (speciesCountWrap.compressed_deflate()) delete[] speciesCountData;
        }
        else throw ConsistencyException("Finished Forward Flux phase n>0 trajectory %llu has recorded %d backward flux events and %d forward flux pilot events; it should have either 1 forward or 1 backward flux event, and not both", trajectoryState.trajectory_id(), timeWrapBackwardFlux.size(), timeWrapForwardFlux.size());

        if (timeWrapBackwardFlux.compressed_deflate()) delete[] timeDataBackwardFlux;
        if (timeWrapForwardFlux.compressed_deflate()) delete[] timeDataForwardFlux;

        msgPtr->set_variance(phaseWeightSV.variance());
    }

    const EndPointVector::Pair& getEndPointCyclic(size_t index) const
    {
        // take the modulus of index to wrap it back around to somewhere within the bounds of endPointVector
        index %= endPointVector.size();
        return endPointVector[index];
    }

    const EndPointVector::Pair& getEndPointUniformRandom() const
    {
        uint32_t ri = getRandomIndex();
        return endPointVector[ri];
    }

    const EndPointVector::Pair& getEndPointCyclic_wellSampled(size_t index) const
    {
        int i = 0;
        while (true)
        {
            size_t ci = (index + i) % endPointVector.size();
            const EndPointVector::Pair& endPointPair = endPointVector[ci];

            if (endPointPair.first->count() > 2 or i > 100)
            {
                // try to get a well-sampled endpoint, but don't loop forever
                return endPointPair;
            }
            i++;
        }
    }

    const EndPointVector::Pair& getEndPointUniformRandom_wellSampled() const
    {
        int i = 0;
        while (true)
        {
            uint32_t ri = getRandomIndex();
            const EndPointVector::Pair& endPointPair = endPointVector[ri];

            if (endPointPair.first->count() > 2 or i > 100)
            {
                // try to get a well-sampled endpoint, but don't loop forever
                return endPointPair;
            }
            i++;
        }
    }

    void final_trajectory_id()
    {
        msgPtr->final_trajectory_id();
    }

    void first_trajectory_id()
    {
        msgPtr->first_trajectory_id();
    }

    EndPointMsg* processEndPoint(int32_t* endPointSpeciesCountsStart, int32_t* endPointSpeciesCountsEnd, double endPointTime, EndPointMap* endPointMap)
    {
        pointKey.assign(endPointSpeciesCountsStart, endPointSpeciesCountsEnd);

        EndPointMsg* endPointMsg = (*endPointMap)[pointKey];
        endPointMsg->set_count(endPointMsg->count() + 1);
        endPointMsg->add_times(endPointTime);

        return endPointMsg;
    }

    EndPointMsg* processSuccessfulEndPoint(int32_t* endPointSpeciesCountsStart, int32_t* endPointSpeciesCountsEnd, double endPointTime)
    {
        return processEndPoint(endPointSpeciesCountsStart, endPointSpeciesCountsEnd, endPointTime, &successfulEndPointMap);
    }

    StartPointMsg* processStartPoint(int32_t* endPointSpeciesCountsStart, int32_t* endPointSpeciesCountsEnd, double endPointTime, const lm::trajectory::Trajectory& trajectory, bool successful)
    {
        const lm::ffpilot::FFPilotTrajectory& ffpilotTraj = static_cast<const lm::ffpilot::FFPilotTrajectory&>(trajectory);

        StartPointMsg* startPointMsg = startPointMap[ffpilotTraj.getInitialSpeciesCounts()];
        startPointMsg->set_count(startPointMsg->count() + 1);
        //startPointMsg->add_times(endPointTime);

        if (not successful)
        {
            // trajectory fluxed backwards
            startPointMsg->set_failed_trajectories_launched_count(startPointMsg->failed_trajectories_launched_count() + 1);
            //startPointMsg->set_failed_trajectories_launched_total_time(startPointMsg->failed_trajectories_launched_total_time() + endPointTime - trajectory.getLastTime());

            startPointEndPointMap.setWrappedField(startPointMsg->mutable_failed_trajectory_end_points());
        }
        else
        {
            // trajectory fluxed forwards
            startPointMsg->set_successful_trajectories_launched_count(startPointMsg->successful_trajectories_launched_count() + 1);
            //startPointMsg->set_successful_trajectories_launched_total_time(startPointMsg->successful_trajectories_launched_total_time() + endPointTime - trajectory.getLastTime());

            startPointEndPointMap.setWrappedField(startPointMsg->mutable_successful_trajectory_end_points());
        }
        //processEndPoint(endPointSpeciesCountsStart, endPointSpeciesCountsEnd, endPointTime, &startPointEndPointMap);
        startPointEndPointMap.setWrappedFieldNull();

        return startPointMsg;
    }

    StartPointMsg* processFailedStartPoint(int32_t* endPointSpeciesCountsStart, int32_t* endPointSpeciesCountsEnd, double endPointTime, const lm::trajectory::Trajectory& trajectory)
    {
        return processStartPoint(endPointSpeciesCountsStart, endPointSpeciesCountsEnd, endPointTime, trajectory, false);
    }

    StartPointMsg* processSuccessfulStartPoint(int32_t* endPointSpeciesCountsStart, int32_t* endPointSpeciesCountsEnd, double endPointTime, const lm::trajectory::Trajectory& trajectory)
    {
        return processStartPoint(endPointSpeciesCountsStart, endPointSpeciesCountsEnd, endPointTime, trajectory, true);
    }

    void set_final_trajectory_id(uint64_t trajectory_id)
    {
        msgPtr->set_final_trajectory_id(trajectory_id);
    }

    void set_first_trajectory_id(uint64_t trajectory_id)
    {
        msgPtr->set_first_trajectory_id(trajectory_id);
    }

    void setWrappedMsg(Msg* newMsgMutablePtr)
    {
        msgPtr = newMsgMutablePtr;

        phaseWeightSV.clear();

        startPointMap.setWrappedField(wrappedMsg()->mutable_start_points());
        successfulEndPointMap.setWrappedField(wrappedMsg()->mutable_successful_trajectory_end_points());

        rebuildEndPointVector();
    }

    void setWrappedMsgNull()
    {
        msgPtr = NULL;

        phaseWeightSV.clear();

        startPointMap.setWrappedFieldNull();
        successfulEndPointMap.setWrappedFieldNull();

        endPointVector.clear();
    }

    Msg* wrappedMsg()
    {
        return msgPtr;
    }

protected:
    void initRng() const {destructRng(); rng = new lm::rng::XORShift(0, 0);}
    void initRng(int cudaDevice) const
    {
        destructRng();

        #ifdef OPT_CUDA
        // Create the cuda based rng.
        rng = new lm::rng::XORWow(cudaDevice, 0, 0, lm::rng::RandomGenerator::UNIFORM);
        #endif

        if (rng == NULL)
        {
            rng = new lm::rng::XORShift(0, 0);
        }
    }
    void initRandomIndexes(size_t size) const
    {
        initRng();
        destructRandomIndexes();

        randomCacheSize = size;
        randomIndexesStart = new uint32_t[randomCacheSize];
        randomDoublesStart = new double[randomCacheSize];
        
        // set range pointers to the ends of the arrays, to match the start pointers
        randomIndexesEnd = randomIndexesStart + randomCacheSize;
        randomDoublesEnd = randomDoublesStart + randomCacheSize;
    }

    uint32_t getRandomIndex() const
    {
        // refill the cache of random numbers, if needed
        if (randomIndexes==randomIndexesEnd) fillRandomIndex();

        // if endPointVector has changed in size since the last time we ran .fillRandomIndex(), rescale a randomDouble on the fly to get a randomIndex
        if (randomIndexesDirty)
        {
            randomIndexes++;
            return getIndexFromDouble(*randomDoubles++);
        }
            // otherwise, our cache of randomIndexes is still good, so just take from that
        else
        {
            randomDoubles++;
            return *randomIndexes++;
        }
    }

    inline uint32_t getIndexFromDouble(double d) const
    {
        uint32_t ri = static_cast<uint32_t>(floor(d*endPointVector.size()));
        if (ri >= endPointVector.size())
        {
            throw ConsistencyException("FFPilotPhaseOutputWrap generated a random index outside the bounds of its endPointVector. ri: %d, endPointVector.size(): %d", ri, endPointVector.size());
        }
        return ri;
    }

    void fillRandomIndex() const
    {
        // initialize the rng stuff for this instance of FFPilotPhaseOutputWrap, if needed
        if (randomIndexes==NULL) {initRandomIndexes(randomCacheSize);}

        // get a large quantity of random doubles
        rng->getRandomDoubles(randomDoublesStart, randomCacheSize);

        // set both randomIndexes and randomDoubles to the front of their arrays
        randomIndexes = randomIndexesStart;
        randomDoubles = randomDoublesStart;
        
        // convert the random doubles into random uints that can be used to randomly lookup values in our table of successful trajectory endpoints
        for (;randomIndexes!=randomIndexesEnd and randomDoubles!=randomDoublesEnd;randomIndexes++,randomDoubles++)
        {
            *randomIndexes = getIndexFromDouble(*randomDoubles);
        }

        // randomIndexes has been rebuilt according the current endPointVector.size(), so mark that they match
        randomIndexesDirty = false;

        // reset the ptrs
        randomIndexes = randomIndexesStart;
        randomDoubles = randomDoublesStart;
    }

    void rebuildEndPointVector()
    {
        size_t oldSize = endPointVector.size();
        endPointVector.clear();

        // iterate over all of the endpoints in the successful_trajectory_end_point repeated field
        for (EndPointMap::iterator epit=successfulEndPointMap.begin();epit!=successfulEndPointMap.end();epit++)
        {
            // add an entry to endPointVector for every time a particular endpoint was "seen"
            for (int i=0;i<epit->count();i++)
            {
                endPointVector.push_back(std::make_pair(&*epit, i));
            }
        }

        if (endPointVector.size()!=oldSize) randomIndexesDirty = true;
    }

    void destructRng() const {if (rng!=NULL) delete rng; rng=NULL;}
    void destructRandomIndexes() const
    {
        if (randomIndexesStart!=NULL) delete[] randomIndexesStart; randomIndexesStart = NULL;
        randomIndexes = NULL;
        randomIndexesEnd = NULL;

        if (randomDoublesStart!=NULL) delete[] randomDoublesStart; randomDoublesStart = NULL;
        randomDoubles = NULL;
        randomDoublesEnd = NULL;

        randomCacheSize = 0;
    }

public:
    EndPointMap successfulEndPointMap;

    StartPointMap startPointMap;
    EndPointMap startPointEndPointMap;

    StreamingVariance phaseWeightSV;

protected:
    Msg* msgPtr;
    PointKey pointKey;
    lm::protowrap::NDArray<int32_t> speciesCountWrap;
    lm::protowrap::NDArray<double> timeWrapForwardFlux, timeWrapBackwardFlux, timeWrapBasinExit;
    lm::protowrap::Repeated<lm::io::LimitTracking> limitTrackingsWrap;

    // list of pointers into the successful_trajectory_end_points field. Part of the system used to randomly choose some of them
    EndPointVector::T endPointVector;

    // rng used for randomly choosing points from one of the point lists. Caches large quantities of random numbers in an attempt to reduce the turnaround time of WorkUnitFinished messages on the supervisor
    mutable lm::rng::RandomGenerator * rng;
    mutable double *randomDoublesStart, *randomDoubles, *randomDoublesEnd;
    mutable uint32_t *randomIndexesStart, *randomIndexes, *randomIndexesEnd;
    mutable size_t randomCacheSize;

    // flag that indicates if the length of endPointVector has changed since we last refilled randomIndexes
    mutable bool randomIndexesDirty;
};

}
}


#endif /* LM_PROTOWRAP_FFLUXPHASEOUTPUT_H_ */
