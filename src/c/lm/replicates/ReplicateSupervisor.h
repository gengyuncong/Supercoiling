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

#ifndef LM_REPLICATES_REPLICATESUPERVISOR_H_
#define LM_REPLICATES_REPLICATESUPERVISOR_H_

#include "lm/trajectory/TrajectoryListSupervisor.h"

namespace lm {
namespace replicates {

class ReplicateSupervisor : public lm::trajectory::TrajectoryListSupervisor
{
public:
    static bool registered;
    static bool registerClass();
    static void* allocateObject();

public:
    ReplicateSupervisor();
    virtual ~ReplicateSupervisor();
    virtual std::string getClassName();
    virtual void printPerformanceStatistics();

protected:
    virtual void startSimulation();
    virtual void buildTrajectoryList();
    virtual void buildTrajectoryOptions();
    virtual void finishSimulation();

private:
    uint64_t numberReplicates;
};

}
}

#endif /* LM_REPLICATES_REPLICATESUPERVISOR_H_ */
