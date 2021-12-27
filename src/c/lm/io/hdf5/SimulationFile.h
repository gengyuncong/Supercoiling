/*
 * University of Illinois Open Source License
 * Copyright 2008-2011 Luthey-Schulten Group,
 * Copyright 2012-2014 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Luthey-Schulten Group
 * 			     University of Illinois at Urbana-Champaign
 * 			     http://www.scs.uiuc.edu/~schulten
 *
 * Developed by: Roberts Group
 * 			     Johns Hopkins University
 * 			     http://biophysics.jhu.edu/roberts/
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
 * Urbana-Champaign, the Roberts Group, Johns Hopkins University, nor the names
 * of its contributors may be used to endorse or promote products derived from
 * this Software without specific prior written permission.
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

#ifndef LM_IO_HDF5_SIMULATIONFILE_H_
#define LM_IO_HDF5_SIMULATIONFILE_H_

#include <google/protobuf/repeated_field.h>
#include <map>
#include <string>
#include <vector>

#include "lm/Exceptions.h"
#include "lm/Types.h"
#include "lm/input/CMERestart.pb.h"
#include "lm/input/DiffusionModel.pb.h"
#include "lm/input/RDMERestart.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/input/SpatialModel.pb.h"
#include "lm/input/ffpilot/FFPilotOptions.pb.h"
#include "lm/io/DegreeAdvancementTimeSeries.pb.h"
#include "lm/io/SpeciesTimeSeries.pb.h"
#include "lm/io/hdf5/HDF5.h"
#include "lm/io/ffpilot/FFPilotOutput.pb.h"
#include "lm/types/OrderParameters.pb.h"
#include "lm/types/Tilings.pb.h"
#include "robertslab/pbuf/NDArray.pb.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using robertslab::pbuf::NDArraySerializer;

namespace lm {

namespace input {
class DiffusionModel;
class ReactionModel;
class SpatialModel;
class Tilings;
}

namespace rdme{
class Lattice;
}

namespace types {
class BoundaryConditions;
class Lattice;
}


namespace io {
class FirstPassageTimes;
class FFPilotOutput;
class LatticeTimeSeries;
class OrderParameterFirstPassageTimes;
class ParameterValues;
class SpeciesCounts;
class SpeciesTimeSeries;
class TilingHist;

namespace hdf5 {

using std::string;
using std::map;
using std::vector;

using lm::IOException;

typedef struct {
    lm::types::OrderParameters* orderParameters;
} CallbackDataOrderParameters;

typedef struct {
    lm::types::Tilings* tilings;
} CallbackDataTilings;

class SimulationFile
{
public:
    SimulationFile();
    virtual ~SimulationFile();
};

class Hdf5File : public SimulationFile
{
public:
    // struct that handles the set of arguments required to write out a dataset
    struct DatasetDescriptor
    {
    public:
        DatasetDescriptor(const std::string& groupPath, const std::string& datasetName, const utuple& shape, hid_t hdf5Type, void* data, hid_t rootGroup=-1);
        DatasetDescriptor(const std::string& groupPath, const std::string& datasetName, const robertslab::pbuf::NDArray& dataWrapper, hid_t rootGroup=-1);
        ~DatasetDescriptor();
        hid_t rootGroup;
        std::string groupPath;
        std::string datasetName;
        utuple shape;
        uint startingColumn;
        hid_t hdf5Type;
        byte* data;
        bool ownsData;
    private:
        template <typename dataT> void initalizeData(const robertslab::pbuf::NDArray& dataWrapper);
    };

    struct ReplicateHandles
    {
        hid_t group;
        hid_t speciesCountsDataset, speciesCountTimesDataset;
        ReplicateHandles():group(H5I_INVALID_HID),speciesCountsDataset(H5I_INVALID_HID),speciesCountTimesDataset(H5I_INVALID_HID) {}
    };
    typedef PairMap<string, uint64_t, ReplicateHandles *> ReplicateHandleMap;

public:
    static const uint MIN_VERSION;
    static const uint CURRENT_VERSION;
    static const uint MAX_REACTION_RATE_CONSTANTS;
    static const uint MAX_SHAPE_PARAMETERS;
    static const uint MAX_NAME_LENGTH = 256;

public:
    static bool isValidFile(const string filename);
    static bool isValidFile(const char * filename);
    static void create(const string filename);
    static void create(const char *  filename);
    static void create(const string filename, unsigned int numberSpecies);
    static void create(const char *  filename, unsigned int numberSpecies);
    static void create(const char * filename, bool initializeModel, unsigned int numberSpecies=0);

    Hdf5File(const string filename);
    Hdf5File(const char* filename);
	virtual ~Hdf5File();
    virtual void close();
    virtual string checkpoint();
    virtual void flush();
    virtual hid_t initGroup(const std::vector<std::string>& groupPathVector, hid_t rootGroup=-1);
    virtual hid_t initGroup(const std::string& groupPath, hid_t rootGroup=-1);

    // Methods for working with parameters.
    virtual bool hasParameter(string key) const;
    virtual string getParameter(string key, string defaultValue="") const;
    virtual void setParameter(string key, string value);

    // Methods for working with the model.
    virtual bool hasDiffusionModel() const;
    virtual void getDiffusionModel(lm::input::DiffusionModel* diffusionModel) const;
    virtual void setDiffusionModel(lm::input::DiffusionModel* diffusionModel);
    virtual bool hasBoundaryGradient() const;
    virtual void getBoundaryGradient(lm::types::BoundaryConditions* bc) const;
    virtual bool hasOrderParameters() const;
    virtual void getOrderParameters(lm::types::OrderParameters* orderParameters) const;
    virtual void setOrderParameters(lm::types::OrderParameters* orderParameters);
    virtual bool hasReactionModel() const;
    virtual void getReactionModel(lm::input::ReactionModel* reactionModel) const;
    virtual void setReactionModel(lm::input::ReactionModel* reactionModel);
    virtual void setSpatialModel(lm::input::SpatialModel* model);
    virtual void getSpatialModel(lm::input::SpatialModel* model) const;
    virtual bool hasTilings() const;
    virtual void getTilings(lm::types::Tilings* tilings) const;
    virtual void setTilings(lm::types::Tilings* tilings);
    virtual bool hasFFPilotOptions() const;
    virtual void getFFPilotOptions(lm::input::ffpilot::FFPilotOptions* options) const;
    virtual void setFFPilotOptions(lm::input::ffpilot::FFPilotOptions* options);

    // Methods for working with restart data.
    virtual bool hasReactionRestart() const;
    virtual void getReactionRestart(lm::input::CMERestart* restart) const;
    virtual bool hasDiffusionRestart() const;
    virtual void getDiffusionRestart(lm::input::RDMERestart* restart, uint numberTrajectories) const;

    // Methods for working with a replicate.
    virtual bool replicateExists(uint64_t replicate);
    virtual void openReplicate(uint64_t replicate);
    virtual void appendSpeciesCounts(uint64_t replicate, lm::io::SpeciesCounts* speciesCounts);
    virtual void appendSpeciesTimeSeries(uint64_t replicate, const lm::io::SpeciesTimeSeries& speciesCounts);
    virtual void appendSpeciesTimeSeries(uint64_t replicate, int numberEntries, int numberSpecies, const int32_t* counts, const double* times);
    virtual void appendLatticeTimeSeries(uint64_t replicate, const lm::io::LatticeTimeSeries& data);
    virtual void appendParameterValues(uint64_t replicate, lm::io::ParameterValues* parameterValues);
    virtual void setFirstPassageTimes(uint64_t replicate, const lm::io::FirstPassageTimes& speciesCounts);
    virtual vector<double> getLatticeTimes(uint64_t replicate);
    virtual void getLattice(uint64_t replicate, unsigned int latticeIndex, lm::rdme::Lattice* lattice);
    virtual void closeReplicate(uint64_t replicate);
    virtual void closeAllReplicates();

    // Methods for working with ffpilot output.
    virtual void setFFPilotOutput(const lm::io::ffpilot::FFPilotOutput& ffpilotOutput);

    virtual void setTilingHist(lm::io::TilingHist* tilingHist, std::string datasetName, hid_t superGroup);

    virtual void setRecordNamePrefix(const string& newRecordNamePrefix);

    //virtual void appendSpatialModelObjects(uint64_t replicate, lm::input::SpatialModel * model) throw(HDF5Exception,InvalidArgException);
    //virtual void getSpatialModelObjects(uint64_t replicate, lm::input::SpatialModel * model) throw(HDF5Exception);

    // Methods for working with NDArrays
    hsize_t setDatasetFromNDArray(const std::string& groupPath, const std::string& datasetName, const robertslab::pbuf::NDArray& ndarrayRef, hid_t rootGroup = -1);
    void setDatasetFromNDArrayReplicate(uint64_t replicate, const std::string& groupRelativePath, const std::string& datasetName, const robertslab::pbuf::NDArray& ndarray, bool condensed=false);
    void setDatasetFromNDArrayReplicateCondensed(uint64_t replicate, const std::string& groupRelativePath, const std::string& datasetName, const robertslab::pbuf::NDArray& ndarray);
    template <typename Container> hsize_t setDatasetFromContainer(const std::string& groupPath, const std::string& datasetName, const Container& container, hid_t rootGroup = -1)
    {
        utuple shape(container.size());
        hid_t hdf5Type = HDF5Type<typename Container::value_type>::T();

        return setDataset(DatasetDescriptor(groupPath, datasetName, shape, hdf5Type, (void*)container.data(), rootGroup));
    }

    // low(ish)-level methods for outputing abstract multi-dimensional array (ie a pointer plus a shape) as a dataset
    hsize_t setDataset(const DatasetDescriptor& dd);
//    void setDatasets(std::vector<DatasetDescriptor>* ddVector);

protected:
    static herr_t parseParameter(hid_t location_id, const char *attr_name, const H5A_info_t *ainfo, void *op_data);
    static herr_t getOrderParametersCallback (hid_t loc_id, const char *name, const H5L_info_t *info, void *callbackDataOrderParameters);
    static herr_t getTilingsCallback (hid_t loc_id, const char *name, const H5L_info_t *info, void *callbackDataTilings);
    static herr_t getDiffusionRestartCallback(hid_t loc_id, const char* name, const H5L_info_t* info, void* obj);

    virtual void open();
    virtual void openGroups();
    virtual void loadParameters();
    virtual void loadModel();
    virtual ReplicateHandles* openReplicateHandles(uint64_t replicate);
    virtual ReplicateHandles* createReplicateHandles(string replicateString);
    virtual void closeReplicateHandles(ReplicateHandles * handles);
	
protected:
    string          filename;
    hid_t           file;
    unsigned int    version;

    // Main group handles.
    hid_t           modelGroup, simulationsGroup;

    // for HDF5 output files, recordNamePrefix is used to
    string recordNamePrefix;

    // The parameters.
    map<string,string> parameterMap;

    // The model.
    bool            modelLoaded;
    unsigned int    numberSpecies;
    unsigned int    numberReactions;

    // Handles for each replicate that is open.
    ReplicateHandleMap::T openReplicates;
};

}
}
}

#endif /* LM_IO_HDF5_SIMULATIONFILE_H_ */
