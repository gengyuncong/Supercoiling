/*
 * Copyright 2012-2019 Johns Hopkins University
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

#ifndef LM_TRAJECTORY_TRAJECTORY_H_
#define LM_TRAJECTORY_TRAJECTORY_H_

#include <map>
#include <string>
#include <vector>

#include "lm/Types.h"
#include "lm/io/TrajectoryState.pb.h"
#include "lm/types/TrajectoryLimits.pb.h"

namespace lm {
namespace trajectory {

class Trajectory
{
public:
    // NB: any changes made to the enum Status *must* be made also to the array status_strings
    enum Status {ABORTED,
                 FINISHED,
                 NOT_STARTED,
                 RUNNING,
                 WAITING};
    static const std::string status_strings[];

    Trajectory(uint64_t id);
    virtual ~Trajectory();

    // accessors
    virtual uint64_t getID() const;
    virtual Status getStatus() const;
    virtual const lm::types::TrajectoryLimit& getLimitReached() const;
    virtual const lm::io::TrajectoryState& getState() const;
    virtual lm::io::TrajectoryState* getMutableState();
    virtual int64_t getWorkUnitsPerformed() const;

    // mutators
    virtual void setStatus(Status newStatus);
    virtual void setLimitReached(const lm::types::TrajectoryLimit& newLimitReached);
    virtual void setState(const lm::io::TrajectoryState& newState);
    virtual void incrementWorkUnitsPerformed();

protected:
    uint64_t id;
    int64_t workUnitsPerformed;
    lm::io::TrajectoryState state;
    lm::types::TrajectoryLimit limitReached;
    Status status;
};

}
}

#endif
