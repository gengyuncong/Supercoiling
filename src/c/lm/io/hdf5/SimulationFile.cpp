/*
 * University of Illinois Open Source License
 * Copyright 2008-2011 Luthey-Schulten Group,
 * Copyright 2012-2019 Roberts Group,
 * All rights reserved.
 *
 * Developed by: Luthey-Schulten Group
 *                  University of Illinois at Urbana-Champaign
 *                  http://www.scs.uiuc.edu/~schulten
 *
 * Developed by: Roberts Group
 *                  Johns Hopkins University
 *                  http://biophysics.jhu.edu/roberts/
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
#include <cstdio>
#include <cstring>
#include <list>
#include <map>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>
#include <zlib.h>

#include <google/protobuf/repeated_field.h>

#include "lm/Exceptions.h"
#include "lm/Math.h"
#include "lm/Print.h"
#include "lm/String.h"
#include "lm/Tune.h"
#include "lm/Types.h"
#include "lm/input/DiffusionModel.pb.h"
#include "lm/input/ReactionModel.pb.h"
#include "lm/input/SpatialModel.pb.h"
#include "lm/io/FirstPassageTimes.pb.h"
#include "lm/io/LatticeTimeSeries.pb.h"
#include "lm/io/ParameterValues.pb.h"
#include "lm/io/SpeciesCounts.pb.h"
#include "lm/io/SpeciesTimeSeries.pb.h"
#include "lm/io/TilingHist.pb.h"
#include "lm/io/hdf5/HDF5.h"
#include "lm/io/hdf5/SimulationFile.h"
#include "lm/io/ffpilot/FFPilotOutput.pb.h"
#include "lm/rdme/Lattice.h"
#include "lm/types/Lattice.pb.h"
#include "lm/types/OrderParameters.pb.h"
#include "lm/types/Tilings.pb.h"
#include "robertslab/pbuf/NDArraySerializer.h"

using std::list;
using std::map;
using std::stringstream;
using std::string;
using std::vector;
using lm::IOException;
using robertslab::pbuf::NDArraySerializer;

namespace lm {
namespace io {
namespace hdf5 {

SimulationFile::SimulationFile()
{
}

SimulationFile::~SimulationFile()
{
}

Hdf5File::DatasetDescriptor::DatasetDescriptor(const std::string& groupPath, const std::string& datasetName, const utuple& shape, hid_t hdf5Type, void* data, hid_t rootGroup)
:rootGroup(rootGroup),groupPath(groupPath),datasetName(datasetName),shape(shape),startingColumn(0),hdf5Type(hdf5Type),data(static_cast<byte*>(data)),ownsData(false)
{
}

Hdf5File::DatasetDescriptor::DatasetDescriptor(const std::string& groupPath, const std::string& datasetName, const robertslab::pbuf::NDArray& dataWrapper, hid_t rootGroup)
:rootGroup(rootGroup),groupPath(groupPath),datasetName(datasetName),shape(dataWrapper.shape()),startingColumn(0),hdf5Type(-1),data(NULL),ownsData(false)
{
    if (dataWrapper.data_type() == robertslab::pbuf::NDArray::int32)
        initalizeData<int32_t>(dataWrapper);
    else if (dataWrapper.data_type() == robertslab::pbuf::NDArray::int64)
        initalizeData<int64_t>(dataWrapper);
    else if (dataWrapper.data_type() == robertslab::pbuf::NDArray::uint32)
        initalizeData<uint32_t>(dataWrapper);
    else if (dataWrapper.data_type() == robertslab::pbuf::NDArray::uint64)
        initalizeData<uint64_t>(dataWrapper);
    else if (dataWrapper.data_type() == robertslab::pbuf::NDArray::float32)
        initalizeData<float>(dataWrapper);
    else if (dataWrapper.data_type() == robertslab::pbuf::NDArray::float64)
        initalizeData<double>(dataWrapper);
}

Hdf5File::DatasetDescriptor::~DatasetDescriptor()
{
    if (ownsData && data == NULL)
    {
        delete[] data;
        data = NULL;
        ownsData = false;
    }
}

template <typename dataT> void Hdf5File::DatasetDescriptor::initalizeData(const robertslab::pbuf::NDArray& dataWrapper)
{
    hdf5Type = HDF5Type<dataT>::T();
    data = new byte[shape.product()*sizeof(dataT)];
    NDArraySerializer::deserializeInto<dataT>(reinterpret_cast<dataT*>(data), shape, dataWrapper);
    ownsData = true;
}

const uint Hdf5File::MIN_VERSION                   = 2;
const uint Hdf5File::CURRENT_VERSION               = 4;
const uint Hdf5File::MAX_REACTION_RATE_CONSTANTS   = 10;
const uint Hdf5File::MAX_SHAPE_PARAMETERS          = 10;

Hdf5File::Hdf5File(const string filename)
:filename(filename),file(H5I_INVALID_HID),version(0),modelGroup(H5I_INVALID_HID),
 simulationsGroup(H5I_INVALID_HID),recordNamePrefix(""),modelLoaded(false),numberSpecies(0)
{
    open();
}

Hdf5File::Hdf5File(const char* filename)
:filename(filename),file(H5I_INVALID_HID),version(0),modelGroup(H5I_INVALID_HID),
 simulationsGroup(H5I_INVALID_HID),recordNamePrefix(""),modelLoaded(false),numberSpecies(0)
{
    open();
}

Hdf5File::~Hdf5File()
{
    //Close the file, if it is still open.
    close();
}

void Hdf5File::open()
{
    // Make sure gzip is supported.
    unsigned int filter_info;
    if (!H5Zfilter_avail(H5Z_FILTER_DEFLATE)) throw lm::Exception("The HDF5 library does not support gzip compression.");
    HDF5_EXCEPTION_CHECK(H5Zget_filter_info (H5Z_FILTER_DEFLATE, &filter_info));
    if (!(filter_info & H5Z_FILTER_CONFIG_ENCODE_ENABLED) || !(filter_info & H5Z_FILTER_CONFIG_DECODE_ENABLED)) throw lm::Exception("The HDF5 library does not support gzip filtering for both encoding and decoding.");

    // Turn of error printing.
    HDF5_EXCEPTION_CHECK(H5Eset_auto2(H5E_DEFAULT, NULL, NULL));

    // Make sure the file exists.
    struct stat fileStats;
    if (stat(filename.c_str(), &fileStats) != 0) throw lm::IOException("the specified file did not exist", filename.c_str());

    // Make sure it is a regular file.
    if (!S_ISREG(fileStats.st_mode)) throw lm::IOException("the specified file was not a regular file", filename.c_str());

    // Make sure the file has the correct magic and hdf headers.
    if (!isValidFile(filename)) throw lm::IOException("the specified file was not of the correct format", filename.c_str());

    // Open the file.
    lm::Print::printfStart(Print::DEBUG, "Opening hdf5 file %s...", filename.c_str());
    HDF5_EXCEPTION_CALL(file,H5Fopen(filename.c_str(), H5F_ACC_RDWR, H5P_DEFAULT));
    lm::Print::printfEnd(Print::DEBUG, " %s opened, %d associated open hdf5 objects (negative means none).", filename.c_str(), H5Fget_obj_count(file, H5F_OBJ_ALL));

    // Get the size of the user block.
    hid_t creationProperties;
    hsize_t userblockSize;
    HDF5_EXCEPTION_CALL(creationProperties,H5Fget_create_plist(file));
    HDF5_EXCEPTION_CHECK(H5Pget_userblock(creationProperties, &userblockSize));
    if (userblockSize < 4) throw lm::IOException("the specified file did not have the correct user block size", userblockSize);

    // Make sure the version is supported.
    HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/", "formatVersion", &version));
    if (version < MIN_VERSION || version > CURRENT_VERSION) throw lm::IOException("the specified file format version is not supported", version);

    // Open the groups.
    openGroups();

    // Loaded the parameters.
    loadParameters();
}

hid_t Hdf5File::initGroup(const vector<string>& groupPathVector, hid_t rootGroup)
{
    hid_t currentGroup, nextGroup;
    currentGroup = rootGroup>=0 ? rootGroup : file;
    for (vector<string>::const_iterator it = groupPathVector.begin(); it!=groupPathVector.end(); it++)
    {
        if (H5Lexists(currentGroup, it->c_str(), H5P_DEFAULT) == 0)
        {
            HDF5_EXCEPTION_CALL(nextGroup, H5Gcreate2(currentGroup, it->c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        }
        else
        {
            HDF5_EXCEPTION_CALL(nextGroup, H5Gopen2(currentGroup, it->c_str(), H5P_DEFAULT));
        }
        currentGroup = nextGroup;
    }
    return currentGroup;
}

hid_t Hdf5File::initGroup(const string& groupPath, hid_t rootGroup)
{
    vector<string> groupPathVector;
    std::stringstream ss(groupPath);
    std::string item;
    while (std::getline(ss, item, '/'))
    {
        if (item.size() > 0)
        {
            groupPathVector.push_back(item);
        }
    }

    return initGroup(groupPathVector, rootGroup);
}

void Hdf5File::openGroups()
{
    HDF5_EXCEPTION_CALL(modelGroup,H5Gopen2(file, "/Model", H5P_DEFAULT));
    HDF5_EXCEPTION_CALL(simulationsGroup,H5Gopen2(file, "/Simulations", H5P_DEFAULT));
}

void Hdf5File::flush()
{
    if (file != H5I_INVALID_HID)
    {
        lm::Print::printf(Print::DEBUG, "Flushing file %s.", filename.c_str());
        HDF5_EXCEPTION_CHECK(H5Fflush(file, H5F_SCOPE_GLOBAL));
    }
}

string Hdf5File::checkpoint()
{
    // Close the file.
    close();

    // Copy the file to a checkpoint file.
    string checkpointFilename(filename+".chk");
    FILE * out=NULL, * in=NULL;
    in=fopen(filename.c_str(), "rb");
    if (in == NULL) throw lm::IOException("the checkpoint input file could not be opened", filename.c_str());
    out=fopen(checkpointFilename.c_str(), "wb");
    if (out == NULL) throw lm::IOException("the checkpoint output file could not be opened", checkpointFilename.c_str());
    size_t bufferSize = 1024*1024;
    char * buffer = new char[bufferSize];
    while (!feof(in))
    {
        size_t bytesRead = fread(buffer, sizeof(char), bufferSize, in);
        if (bytesRead != bufferSize && ferror(in)) throw lm::IOException("the checkpoint input file could not be read", filename.c_str());
        size_t bytesWritten = fwrite(buffer, sizeof(char), bytesRead, out);
        if (bytesWritten != bytesRead) throw lm::IOException("the checkpoint input file could not be read", checkpointFilename.c_str());
    }
    if (fclose(in) != 0)  throw lm::IOException("the checkpoint input file could not be closed", filename.c_str());
    if (fclose(out) != 0)  throw lm::IOException("the checkpoint output file could not be closed", checkpointFilename.c_str());
    delete[] buffer;

    // Open the file again.
    open();

    return checkpointFilename;
}

void Hdf5File::close()
{
    // Close any open replicate groups.
    closeAllReplicates();

    // Close any open groups.
    if (modelGroup != H5I_INVALID_HID)
    {
        HDF5_EXCEPTION_CHECK(H5Gclose(modelGroup));
        modelGroup = H5I_INVALID_HID;
    }
    if (simulationsGroup != H5I_INVALID_HID)
    {
        HDF5_EXCEPTION_CHECK(H5Gclose(simulationsGroup));
        simulationsGroup = H5I_INVALID_HID;
    }

    // Close the file.
    if (file != H5I_INVALID_HID)
    {
        lm::Print::printf(Print::DEBUG, "Closing hdf5 file %s.", filename.c_str(), int(file));
        HDF5_EXCEPTION_CHECK(H5Fclose(file));
        file = H5I_INVALID_HID;
    }
}

bool Hdf5File::isValidFile(const char * filename)
{
    return isValidFile(string(filename));
}

bool Hdf5File::isValidFile(const string filename)
{
    // Make sure the file has the right magic.
    FILE * fp = NULL;
    fp=fopen(filename.c_str(), "r");
    if (fp == NULL) throw lm::IOException("the specified file could not be opened", filename.c_str());
    char magic[5];
    magic[4] = '\0';
    for (int i=0; i<4; i++)
    {
        magic[i] = fgetc(fp);
        if (magic[i] == EOF)
        {
            magic[i] = '\0';
            break;
        }
    }
    if (fclose(fp) != 0)  throw lm::IOException("the specified file could not be closed", filename.c_str());
    if (strncmp(magic, "LMH5", 4) != 0) return false;

    // Make sure it is an hdf5 file.
    int isHdf5 = false;
    HDF5_EXCEPTION_CALL(isHdf5,H5Fis_hdf5(filename.c_str()));
    if (!isHdf5) return false;

    return true;

}

void Hdf5File::loadParameters()
{
    hsize_t n=0;
    hid_t parametersGroup;
    HDF5_EXCEPTION_CALL(parametersGroup,H5Gopen2(file, "/Parameters", H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Aiterate2(parametersGroup, H5_INDEX_CRT_ORDER, H5_ITER_NATIVE, &n, &Hdf5File::parseParameter, this));
    HDF5_EXCEPTION_CHECK(H5Gclose(parametersGroup));
}

herr_t Hdf5File::parseParameter(hid_t location_id, const char *attr_name, const H5A_info_t *ainfo, void *op_data)
{
    Hdf5File * file = reinterpret_cast<Hdf5File *>(op_data);

    // Open the attribute.
    hid_t attr, type;
    H5T_class_t typeClass;
    hsize_t size;
    HDF5_EXCEPTION_CALL(attr,H5Aopen(location_id, attr_name, H5P_DEFAULT));
    HDF5_EXCEPTION_CALL(type,H5Aget_type(attr));
    typeClass=H5Tget_class(type);
    HDF5_EXCEPTION_CALL(size,H5Aget_storage_size(attr));

    if (file->version == 2)
    {
        if (typeClass==H5T_FLOAT && size == sizeof(double))
        {
            double value;
            stringstream ss;
            HDF5_EXCEPTION_CHECK(H5Aread(attr, H5T_NATIVE_DOUBLE, &value));
            ss << value;
            file->parameterMap[attr_name] = ss.str();
        }
        else
        {
            // could not parse the parameter. "Warn" the user by raising an exception
            THROW_EXCEPTION(lm::InputException, "Unable to parse user defined value for parameter: %s\n", attr_name);
        }
    }

    if (file->version >= 3)
    {
        if (typeClass==H5T_INTEGER)
        {
            // read in any integer type as 64-bit int. Will cause problems with very large unsigned ints
            int64_t value;
            stringstream ss;
            HDF5_EXCEPTION_CHECK(H5Aread(attr, type, &value));
            ss << value;
            file->parameterMap[attr_name] = ss.str();
        }
        else if (typeClass==H5T_FLOAT)
        {
            // read in any float type as a double.
            double value;
            stringstream ss;
            HDF5_EXCEPTION_CHECK(H5Aread(attr, type, &value));
            ss << value;
            file->parameterMap[attr_name] = ss.str();
        }
        else if (typeClass==H5T_STRING)
        {
            // Get the dataspace.
            hid_t space;
            HDF5_EXCEPTION_CALL(space,H5Aget_space(attr));

            // hdf5 requires different handling for fixed and variable length string attributes
            htri_t is_variable_len = H5Tis_variable_str(type);

            hid_t memtype;
            char * value;
            if (is_variable_len)
            {
                // Create the memory datatype.
                HDF5_EXCEPTION_CALL(memtype, H5Tget_native_type(type, H5T_DIR_DEFAULT));

                // set the size to variable. Commented since it's unnecessary/causes problems
                //HDF5_EXCEPTION_CHECK(H5Tset_size(memtype, H5T_VARIABLE));

                // Read the data. For variable length strings, pass the output char array as a char**. Apparently H5Aread will also take care of allocation
                HDF5_EXCEPTION_CHECK(H5Aread(attr, memtype, &value));

                // Add the parameter to the map.
                file->parameterMap[attr_name] = value;

                // Reclaim the memory.
                #ifdef OLD_HDFREE
                    free(value);
                #else
                    H5free_memory(value);
                #endif
            }
            else
            {
                // Create the memory datatype.
                HDF5_EXCEPTION_CALL(memtype, H5Tcopy(H5T_C_S1));

                // set the fixed size
                HDF5_EXCEPTION_CHECK(H5Tset_size(memtype, size));

                // allocate space for the fixed string
                value = new char[size];

                // Read the data.
                HDF5_EXCEPTION_CHECK(H5Aread(attr, memtype, value));

                // Add the parameter to the map.
                file->parameterMap[attr_name] = value;

                // Reclaim the memory.
                delete [] value;
            }
            HDF5_EXCEPTION_CHECK(H5Sclose(space));
            HDF5_EXCEPTION_CHECK(H5Tclose(memtype));
        }
        else
        {
            // could not parse the parameter. "Warn" the user by raising an exception
            THROW_EXCEPTION(lm::InputException, "Unable to parse user defined value for parameter: %s\n", attr_name);
        }
    }
    HDF5_EXCEPTION_CHECK(H5Tclose(type));
    HDF5_EXCEPTION_CHECK(H5Aclose(attr));
    return 0;
}

bool Hdf5File::hasParameter(string key) const
{
    return parameterMap.count(key)>0;
}

string Hdf5File::getParameter(string key, string defaultValue) const
{
    // If we didn't find the key, return the default value.
    map<string,string>::const_iterator it = parameterMap.find(key);
    if (it == parameterMap.end()) return defaultValue;

    return it->second;
}

void Hdf5File::setParameter(string key, string value)
{
    // Set the parameter in the map.
    parameterMap[key] = value;

    // Save the parameter in the file.
    if (version <= 2)
    {
        double dvalue = atof(value.c_str());
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_double(file, "/Parameters", key.c_str(), &dvalue, 1));
    }
    else
    {
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_string(file, "/Parameters", key.c_str(), value.c_str()));
    }
}

void Hdf5File::loadModel()
{
    if (!modelLoaded)
    {
        if (version <= 3)
        {
            HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model", "numberSpecies", &numberSpecies));
            modelLoaded = true;
        }
        else
        {
            if (H5Lexists(file, "/Model/Reaction", H5P_DEFAULT) > 0)
            {
                HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Reaction", "numberSpecies", &numberSpecies));
                modelLoaded = true;
            }
            else
            {
                throw Exception("No model has been defined.");
            }
        }
    }
}

bool Hdf5File::hasDiffusionModel() const
{
    return (H5Lexists(file, "/Model/Diffusion", H5P_DEFAULT) > 0);
}

void Hdf5File::getDiffusionModel(lm::input::DiffusionModel* diffusionModel) const
{
    // Make sure the model is not null and then clear it.
    if (diffusionModel == NULL) throw InvalidArgException("diffusionModel", "cannot be null");
    diffusionModel->Clear();

    if (H5Lexists(file, "/Model/Diffusion", H5P_DEFAULT) > 0)
    {
        // Read the diffusion model attributes.
        uint numberSpecies, numberReactions, numberSiteTypes, latticeXSize, latticeYSize, latticeZSize, particlesPerSite;
        double latticeSpacing;
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "numberSpecies", &numberSpecies));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "numberReactions", &numberReactions));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "numberSiteTypes", &numberSiteTypes));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_double(file, "/Model/Diffusion", "latticeSpacing", &latticeSpacing));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "latticeXSize", &latticeXSize));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "latticeYSize", &latticeYSize));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "latticeZSize", &latticeZSize));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Diffusion", "particlesPerSite", &particlesPerSite));

        // Fill in the model.
        diffusionModel->set_number_species(numberSpecies);
        diffusionModel->set_number_reactions(numberReactions);
        diffusionModel->set_number_site_types(numberSiteTypes);
        diffusionModel->set_lattice_spacing(latticeSpacing);

        int ndims;
        hsize_t dims[4];
        H5T_class_t type;
        size_t size;

        // Read the diffusion matrix.
        H5LTget_dataset_info(file, "/Model/Diffusion/DiffusionMatrix", dims, &type, &size);
        if (dims[0] != numberSiteTypes || dims[1] != numberSiteTypes || dims[2] != numberSpecies || size != sizeof(double)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/DiffusionMatrix");
        double * D = new double[numberSiteTypes*numberSiteTypes*numberSpecies];
        H5LTread_dataset_double(file, "/Model/Diffusion/DiffusionMatrix", D);
        for (uint i=0; i<numberSiteTypes*numberSiteTypes*numberSpecies; i++) diffusionModel->add_diffusion_matrix(D[i]);
        delete [] D;

        // If we have reactions, read the reaction location matrix.
        if (numberReactions > 0)
        {
            H5LTget_dataset_info(file, "/Model/Diffusion/ReactionLocationMatrix", dims, &type, &size);
            if (dims[0] != numberReactions || dims[1] != numberSiteTypes || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/ReactionLocationMatrix");
            uint * RL = new uint[numberReactions*numberSiteTypes];
            H5LTread_dataset(file, "/Model/Diffusion/ReactionLocationMatrix", H5T_STD_U32LE, RL);
            for (uint i=0; i<numberReactions*numberSiteTypes; i++) diffusionModel->add_reaction_location_matrix(RL[i]);
            delete [] RL;
        }

        // Read the initial lattice.
        {
            HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(file, "/Model/Diffusion/Lattice", &ndims));
            if (ndims != 4) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/Lattice");
            HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(file, "/Model/Diffusion/Lattice",dims, &type, &size));
            if (size != sizeof(uint8_t)) throw Exception("Invalid data type", filename.c_str(), "/Model/Diffusion/Lattice");
            ndarray<uint8_t>* particles = new ndarray<uint8_t>(utuple(dims[0],dims[1],dims[2],dims[3]));
            HDF5_EXCEPTION_CHECK(H5LTread_dataset(file, "/Model/Diffusion/Lattice", H5T_NATIVE_UINT8, particles->values));
            NDArraySerializer::serializeInto(diffusionModel->mutable_initial_lattice()->mutable_particles(), *particles);
            delete particles;
        }

        // Read the initial lattice sites.
        {
            HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(file, "/Model/Diffusion/LatticeSites", &ndims));
            if (ndims != 3) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/LatticeSites");
            HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(file, "/Model/Diffusion/Lattice",dims, &type, &size));
            if (size != sizeof(uint8_t)) throw Exception("Invalid data type", filename.c_str(), "/Model/Diffusion/LatticeSites");
            ndarray<uint8_t>* sites = new ndarray<uint8_t>(utuple(dims[0],dims[1],dims[2]));
            HDF5_EXCEPTION_CHECK(H5LTread_dataset(file, "/Model/Diffusion/LatticeSites", H5T_NATIVE_UINT8, sites->values));
            NDArraySerializer::serializeInto(diffusionModel->mutable_initial_lattice()->mutable_sites(), *sites);
            delete sites;
        }
    }
}

void Hdf5File::setDiffusionModel(lm::input::DiffusionModel* diffusionModel)
{
    // Validate that the model is consistent.
    if (diffusionModel == NULL) throw InvalidArgException("diffusionModel", "cannot be NULL");
    if (diffusionModel->number_species() == 0) throw InvalidArgException("diffusionModel.number_species", "cannot be zero");
    if (diffusionModel->number_site_types() == 0) throw InvalidArgException("diffusionModel.number_site_types", "cannot be zero");
    if (diffusionModel->diffusion_matrix_size() != (int)(diffusionModel->number_site_types()*diffusionModel->number_site_types()*diffusionModel->number_species())) throw InvalidArgException("diffusion.diffusion_matrix", "inconsistent size");
    if (diffusionModel->reaction_location_matrix_size() != (int)(diffusionModel->number_reactions()*diffusionModel->number_site_types())) throw InvalidArgException("diffusion.reaction_location_matrix", "inconsistent size");
    if (diffusionModel->lattice_spacing() <= 0.0) throw InvalidArgException("diffusionModel.lattice_spacing", "must be greater than zero");
    if (diffusionModel->initial_lattice().particles().shape_size() != 4) throw InvalidArgException("diffusionModel.initial_lattice.particles", "must have four dimensions");
    if (diffusionModel->initial_lattice().particles().shape(0) == 0) throw InvalidArgException("diffusionModel.initial_lattice.particles x size", "cannot be zero");
    if (diffusionModel->initial_lattice().particles().shape(1) == 0) throw InvalidArgException("diffusionModel.initial_lattice.particles y size", "cannot be zero");
    if (diffusionModel->initial_lattice().particles().shape(2) == 0) throw InvalidArgException("diffusionModel.initial_lattice.particles z size", "cannot be zero");
    if (diffusionModel->initial_lattice().particles().shape(3) == 0) throw InvalidArgException("diffusionModel.initial_lattice.particles particle per site", "cannot be zero");
    if (diffusionModel->initial_lattice().sites().shape_size() != 3) throw InvalidArgException("diffusionModel.initial_lattice.sites", "must have three dimensions");
    if (diffusionModel->initial_lattice().sites().shape(0) != diffusionModel->initial_lattice().particles().shape(0)) throw InvalidArgException("diffusionModel.initial_lattice.sites x size", "must be the same as diffusionModel.initial_lattice.particles x size");
    if (diffusionModel->initial_lattice().sites().shape(1) != diffusionModel->initial_lattice().particles().shape(1)) throw InvalidArgException("diffusionModel.initial_lattice.sites y size", "must be the same as diffusionModel.initial_lattice.particles y size");
    if (diffusionModel->initial_lattice().sites().shape(2) != diffusionModel->initial_lattice().particles().shape(2)) throw InvalidArgException("diffusionModel.initial_lattice.sites z size", "must be the same as diffusionModel.initial_lattice.particles z size");

    // If a diffusion model already exists, delete it.
    if (H5Lexists(file, "/Model/Diffusion", H5P_DEFAULT) > 0)
    {
        HDF5_EXCEPTION_CHECK(H5Ldelete(file, "/Model/Diffusion", H5P_DEFAULT));
    }

    // Create the group for the reaction model.
    hid_t group;
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Diffusion", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));

    // Write the attributes.
    numberSpecies = diffusionModel->number_species();
    uint numberReactions = diffusionModel->number_reactions();
    uint numberSiteTypes = diffusionModel->number_site_types();
    double latticeSpacing = diffusionModel->lattice_spacing();
    uint latticeXSize = diffusionModel->initial_lattice().particles().shape(0);
    uint latticeYSize = diffusionModel->initial_lattice().particles().shape(1);
    uint latticeZSize = diffusionModel->initial_lattice().particles().shape(2);
    uint particlesPerSite = diffusionModel->initial_lattice().particles().shape(3);
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "numberSpecies", &numberSpecies, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "numberReactions", &numberReactions, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "numberSiteTypes", &numberSiteTypes, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_double(file, "/Model/Diffusion", "latticeSpacing", &latticeSpacing, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "latticeXSize", &latticeXSize, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "latticeYSize", &latticeYSize, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "latticeZSize", &latticeZSize, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Diffusion", "particlesPerSite", &particlesPerSite, 1));

    // Write the diffusion matrix.
    {
        const unsigned int RANK=3;
        hsize_t dims[RANK];
        dims[0] = numberSiteTypes;
        dims[1] = numberSiteTypes;
        dims[2] = numberSpecies;
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Diffusion/DiffusionMatrix", RANK, dims, H5T_IEEE_F64LE, diffusionModel->diffusion_matrix().data()));
    }

    // Write the reaction location matrix.
    if (numberReactions > 0)
    {
        const unsigned int RANK=2;
        hsize_t dims[RANK];
        dims[0] = numberReactions;
        dims[1] = numberSiteTypes;
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Diffusion/ReactionLocationMatrix", RANK, dims, H5T_STD_U32LE, diffusionModel->reaction_location_matrix().data()));
        printf("Set rlm %d %d %d %d\n", numberReactions, numberSiteTypes, diffusionModel->reaction_location_matrix().data()[0], diffusionModel->reaction_location_matrix().data()[1]);
    }

    // Write the lattice data set.
    {
        ndarray<uint8_t>* particles = NDArraySerializer::deserializeAllocate<uint8_t>(diffusionModel->initial_lattice().particles());
        unsigned int RANK=4;
        hsize_t dims[RANK], chunk[RANK];
        dims[0] = particles->shape[0];
        dims[1] = particles->shape[1];
        dims[2] = particles->shape[2];
        dims[3] = particles->shape[3];
        chunk[0] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,particles->shape[0]);
        chunk[1] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,particles->shape[1]);
        chunk[2] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,particles->shape[2]);
        chunk[3] = particles->shape[3];
        hid_t dataspaceHandle, dcplHandle, datasetHandle;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, NULL));
        HDF5_EXCEPTION_CALL(dcplHandle,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_deflate (dcplHandle, TUNE_LATTICE_GZIP_COMPRESSION_LEVEL));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcplHandle, RANK, chunk));
        HDF5_EXCEPTION_CALL(datasetHandle,H5Dcreate2(file, "/Model/Diffusion/Lattice", H5T_STD_U8LE, dataspaceHandle, H5P_DEFAULT, dcplHandle, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Dwrite(datasetHandle, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, particles->values));
        HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcplHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
        delete particles;
    }

    // Write the lattice sites data set.
    {
        ndarray<uint8_t>* sites = NDArraySerializer::deserializeAllocate<uint8_t>(diffusionModel->initial_lattice().sites());
        unsigned int RANK=3;
        hsize_t dims[RANK], chunk[RANK];
        dims[0] = sites->shape[0];
        dims[1] = sites->shape[1];
        dims[2] = sites->shape[2];
        chunk[0] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,sites->shape[0]);
        chunk[1] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,sites->shape[1]);
        chunk[2] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,sites->shape[2]);
        hid_t dataspaceHandle, dcplHandle, datasetHandle;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, NULL));
        HDF5_EXCEPTION_CALL(dcplHandle,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_deflate (dcplHandle, TUNE_LATTICE_GZIP_COMPRESSION_LEVEL));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcplHandle, RANK, chunk));
        HDF5_EXCEPTION_CALL(datasetHandle,H5Dcreate2(file, "/Model/Diffusion/LatticeSites", H5T_STD_U8LE, dataspaceHandle, H5P_DEFAULT, dcplHandle, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Dwrite(datasetHandle, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, sites->values));
        HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcplHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
        delete sites;
    }
}

bool Hdf5File::hasBoundaryGradient() const
{
    return (H5Lexists(file, "/Model/Diffusion/Gradient", H5P_DEFAULT) > 0);
}

void Hdf5File::getBoundaryGradient(lm::types::BoundaryConditions* bc) const
{
    // Make sure the model is not null and then clear it.
    if (bc == NULL) throw InvalidArgException("bc", "cannot be null");

    if (H5Lexists(file, "/Model/Diffusion/Gradient", H5P_DEFAULT) > 0)
    {
        // Read the lattice size.
        int latticeXSize,latticeYSize,latticeZSize;
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_int(file, "/Model/Diffusion", "latticeXSize", &latticeXSize));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_int(file, "/Model/Diffusion", "latticeYSize", &latticeYSize));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_int(file, "/Model/Diffusion", "latticeZSize", &latticeZSize));

        // Read the gradient.
        int ndims;
        hsize_t dims[3];
        H5T_class_t type;
        size_t size;
        HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(file, "/Model/Diffusion/Gradient", &ndims));
        if (ndims != 3) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/Lattice");
        HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(file, "/Model/Diffusion/Gradient",dims, &type, &size));
        if (latticeXSize+2 != (int)dims[0] || latticeYSize+2 != (int)dims[1] || latticeZSize+2 != (int)dims[2]) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/Gradient");
        if (size != sizeof(double)) throw Exception("Invalid dataset type", filename.c_str(), "/Model/Diffusion/Lattice");
        ndarray<double>* gradient = new ndarray<double>(utuple(dims[0],dims[1],dims[2]));
        HDF5_EXCEPTION_CHECK(H5LTread_dataset(file, "/Model/Diffusion/Gradient", H5T_NATIVE_DOUBLE, gradient->values));
        NDArraySerializer::serializeInto(bc->mutable_boundary_gradient(), *gradient);
        delete gradient;
    }
}

void Hdf5File::setFFPilotOutput(const lm::io::ffpilot::FFPilotOutput& output)
{
    // Make sure that the root group exists and is open.
    bool rootIsFile = true;
    hid_t rootGroup = file;
    if (recordNamePrefix != "")
    {
        if (H5Lexists(file, recordNamePrefix.c_str(), H5P_DEFAULT) == 0)
            HDF5_EXCEPTION_CHECK(H5Gcreate2(file, recordNamePrefix.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        HDF5_EXCEPTION_CALL(rootGroup, H5Gopen2(file, recordNamePrefix.c_str(), H5P_DEFAULT));
        rootIsFile = false;
    }

    // Make sure we have the ffpilot top level group.
    if (H5Lexists(rootGroup, "FFPilot", H5P_DEFAULT) == 0)
        HDF5_EXCEPTION_CHECK(H5Gcreate2(rootGroup, "FFPilot", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    // Make sure we have the output group and open it.
    hid_t outputGroup;
    if (H5Lexists(rootGroup, "FFPilot/Output", H5P_DEFAULT) == 0)
        HDF5_EXCEPTION_CHECK(H5Gcreate2(rootGroup, "FFPilot/Output", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CALL(outputGroup, H5Gopen2(rootGroup, "FFPilot/Output", H5P_DEFAULT));

    // Make sure we have a stage folder for the outputs.
    if (H5Lexists(outputGroup, "Stage", H5P_DEFAULT) == 0)
        HDF5_EXCEPTION_CHECK(H5Gcreate2(outputGroup, "Stage", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    // Write each stage output.
    for (int i=0; i<output.stage_output_size(); i++)
    {
        // Get the stage output.
        const lm::io::ffpilot::FFPilotStageOutput& stageOutput = output.stage_output(i);

        // Make the stage group, deleting any previous version of it.
        hid_t stageGroup;
        char stageName[64];
        snprintf(stageName, 64, "Stage/%d", stageOutput.id());
        if (H5Lexists(outputGroup, stageName, H5P_DEFAULT) > 0)
            HDF5_EXCEPTION_CHECK(H5Ldelete(outputGroup, stageName, H5P_DEFAULT));
        HDF5_EXCEPTION_CALL(stageGroup, H5Gcreate2(outputGroup, stageName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

        // Declare some variables.
        int32_t ivalue;
        const unsigned int RANK=1;
        hsize_t dims[RANK];

        // Set the direction.
        ivalue = stageOutput.direction();
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_int(outputGroup, stageName, "direction", &ivalue, 1));

        // Set the type.
        ivalue = stageOutput.type();
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_int(outputGroup, stageName, "type", &ivalue, 1));

        // Write the edges.
        ndarray<double> edges = NDArraySerializer::deserialize<double>(stageOutput.edges());
        if (edges.shape.len != 1) THROW_EXCEPTION(lm::RuntimeException, "the edges array was not 1D: %d", edges.shape.len);
        dims[0] = edges.shape[0];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(stageGroup, "Edges", RANK, dims, H5T_IEEE_F64LE, edges.values));

        // Write the counts.
        ndarray<uint64_t> counts = NDArraySerializer::deserialize<uint64_t>(stageOutput.trajectory_counts());
        if (counts.shape.len != 1) THROW_EXCEPTION(lm::RuntimeException, "the counts array was not 1D: %d", counts.shape.len);
        if (counts.shape[0] != edges.shape[0]) THROW_EXCEPTION(lm::RuntimeException, "the counts array did not have the correct number of rows: %d,%d", counts.shape[0],edges.shape[0]);
        dims[0] = counts.shape[0];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(stageGroup, "TrajectoryCounts", RANK, dims, H5T_STD_U64LE, counts.values));

        // Write the costs.
        ndarray<double> costs = NDArraySerializer::deserialize<double>(stageOutput.costs());
        if (costs.shape.len != 1) THROW_EXCEPTION(lm::RuntimeException, "the costs array was not 1D: %d", costs.shape.len);
        if (costs.shape[0] != edges.shape[0]) THROW_EXCEPTION(lm::RuntimeException, "the costs array did not have the correct number of rows: %d,%d", costs.shape[0],edges.shape[0]);
        dims[0] = costs.shape[0];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(stageGroup, "Costs", RANK, dims, H5T_IEEE_F64LE, costs.values));

        // Write the weights.
        ndarray<double> weights = NDArraySerializer::deserialize<double>(stageOutput.weights());
        if (weights.shape.len != 1) THROW_EXCEPTION(lm::RuntimeException, "the weights array was not 1D: %d", weights.shape.len);
        if (weights.shape[0] != edges.shape[0]) THROW_EXCEPTION(lm::RuntimeException, "the weights array did not have the correct number of rows: %d,%d", weights.shape[0],edges.shape[0]);
        dims[0] = weights.shape[0];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(stageGroup, "Weights", RANK, dims, H5T_IEEE_F64LE, weights.values));

        // Write the fpts.
        ndarray<double> fpts = NDArraySerializer::deserialize<double>(stageOutput.first_passage_times());
        if (fpts.shape.len != 1) THROW_EXCEPTION(lm::RuntimeException, "the fpts array was not 1D: %d", fpts.shape.len);
        if (fpts.shape[0] != edges.shape[0]) THROW_EXCEPTION(lm::RuntimeException, "the fpts array did not have the correct number of rows: %d,%d", fpts.shape[0],edges.shape[0]);
        dims[0] = fpts.shape[0];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(stageGroup, "FirstPassageTimes", RANK, dims, H5T_IEEE_F64LE, fpts.values));

        // Write the weight variances.
        ndarray<double> weight_variances = NDArraySerializer::deserialize<double>(stageOutput.weight_variances());
        if (weight_variances.shape.len != 1) THROW_EXCEPTION(lm::RuntimeException, "the weight_variances array was not 1D: %d", weight_variances.shape.len);
        if (weight_variances.shape[0] != edges.shape[0]) THROW_EXCEPTION(lm::RuntimeException, "the weight_variances array did not have the correct number of rows: %d,%d", weight_variances.shape[0],edges.shape[0]);
        dims[0] = weight_variances.shape[0];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(stageGroup, "WeightVariances", RANK, dims, H5T_IEEE_F64LE, weight_variances.values));

        // Close the stage group.
        HDF5_EXCEPTION_CHECK(H5Gclose(stageGroup));
    }

    // Close the output group.
    HDF5_EXCEPTION_CHECK(H5Gclose(outputGroup));

    // Close the root group, if it is not the file.
    if (!rootIsFile) HDF5_EXCEPTION_CHECK(H5Gclose(rootGroup));
}

void Hdf5File::setTilingHist(lm::io::TilingHist* tilingHist, std::string datasetName, hid_t superGroup)
{
    if (tilingHist->tile_vals_size() > 0)
    {
        hid_t thGroup;
        hsize_t dims[1];
        uint number_tiles, tiling_id;

        dims[0] = tilingHist->tile_vals_size();
        HDF5_EXCEPTION_CALL(thGroup, H5Gcreate2(superGroup, datasetName.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(thGroup, "TileIndices", 1, dims, H5T_STD_U32LE, tilingHist->tile_indices().data()));
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(thGroup, "TileVals", 1, dims, H5T_IEEE_F64LE, tilingHist->tile_vals().data()));
        number_tiles = tilingHist->number_tiles();
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(superGroup, datasetName.c_str(), "NumberTiles", &number_tiles, 1));
        tiling_id = tilingHist->tiling_id();
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(superGroup, datasetName.c_str(), "TilingID", &tiling_id, 1));
        HDF5_EXCEPTION_CHECK(H5Gclose(thGroup));
    }
}

bool Hdf5File::hasOrderParameters() const
{
    return (H5Lexists(file, "/OrderParameters", H5P_DEFAULT) > 0);
}

herr_t Hdf5File::getOrderParametersCallback(hid_t loc_id, const char * name, const H5L_info_t * info, void * callbackDataOrderParameters)
{
    // Declare and initialize handle for the order parameter group
    hid_t opGroup;
    HDF5_EXCEPTION_CALL(opGroup, H5Gopen(loc_id, name, H5P_DEFAULT));

    // recast the callbackData structure away from void *
    CallbackDataOrderParameters* cdOP = (CallbackDataOrderParameters*)callbackDataOrderParameters;

    // create a new order parameter in the ffpilotParameters protobuf
    lm::types::OrderParameter* newOP = cdOP->orderParameters->add_order_parameter();

    // get the order parameter type and ID
    uint type;
    HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(loc_id, name, "Type", &type));
    newOP->set_type(type);

    // read in the ids of the species that this order parameter uses, as well as their associated coefficients
    hsize_t dims[2];
    H5T_class_t hdf5Type;
    size_t size;

    herr_t speciesIDsExists;
    HDF5_EXCEPTION_CALL(speciesIDsExists, H5LTfind_dataset(opGroup, "SpeciesIDs"))
    if (speciesIDsExists)
    {
        H5LTget_dataset_info(opGroup, "SpeciesIDs", dims, &hdf5Type, &size);
        int* speciesBuffer = new int[dims[0]];
        H5LTread_dataset_int(opGroup, "SpeciesIDs", speciesBuffer);
        for (hsize_t i=0; i < dims[0]; i++)
        {
            newOP->add_species_index(uint32_t(speciesBuffer[i]));
        }
        // free the buffer
        delete[] speciesBuffer;
    }

    herr_t speciesCoefficientsExists;
    HDF5_EXCEPTION_CALL(speciesCoefficientsExists, H5LTfind_dataset(opGroup, "SpeciesCoefficients"))
    if (speciesCoefficientsExists)
    {
        H5LTget_dataset_info(opGroup, "SpeciesCoefficients", dims, &hdf5Type, &size);
        double* coefficientBuffer = new double[dims[0]];
        H5LTread_dataset_double(opGroup, "SpeciesCoefficients", coefficientBuffer);
        for (hsize_t i=0; i < dims[0]; i++)
        {
            newOP->add_species_coefficient(coefficientBuffer[i]);
        }

        // free the buffer
        delete[] coefficientBuffer;
    }

    herr_t speciesExponentsExists;
    HDF5_EXCEPTION_CALL(speciesExponentsExists, H5LTfind_dataset(opGroup, "SpeciesExponents"))
    if (speciesExponentsExists)
    {
        H5LTget_dataset_info(opGroup, "SpeciesExponents", dims, &hdf5Type, &size);
        double* exponentsBuffer = new double[dims[0]];
        H5LTread_dataset_double(opGroup, "SpeciesExponents", exponentsBuffer);
        for (hsize_t i=0; i < dims[0]; i++)
        {
            newOP->add_species_exponent(exponentsBuffer[i]);
        }

        // free the buffer
        delete[] exponentsBuffer;
    }

//    // read in the values of the tiling's basins
//    herr_t basinsExists;
//    HDF5_EXCEPTION_CALL(basinsExists, H5LTfind_dataset(opGroup, "Basins"))
//    if (basinsExists)
//    {
//        // sanity check the rank of the basins dataset
//        int rank;
//        HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(opGroup, "Basins", &rank));
//        if (rank != 2) THROW_EXCEPTION(IOException, "Rank of Basins dataset invalid.\n"
//            "Please ensure that all Basins datasets are 2D in your input .lm file. rank: %d", rank);

//        HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(opGroup, "Basins", dims, &hdf5Type, &size));
//        double* basinsBuffer = new double[dims[0]*dims[1]];
//        HDF5_EXCEPTION_CHECK(H5LTread_dataset_double(opGroup, "Basins", basinsBuffer));

//        for (uint i=0; i<dims[0]; i++)
//        {
//            lm::input::Basin* newBasin = newOP->add_basins();
//            for (uint j=0; j<dims[1]; j++)
//            {
//                newBasin->add_species_count(basinsBuffer[i*dims[1]+j]);
//            }
//        }
//        // free the buffers
//        delete[] basinsBuffer;
//    }

    // free the group handle
    HDF5_EXCEPTION_CHECK(H5Gclose(opGroup));

//    // sanity check the contents of the order parameter we just read in
//    if (not basinsExists and not (speciesIDsExists and speciesCoefficientsExists))
//    {
//        THROW_EXCEPTION(IOException, "Incomplete order parameter in input.\n"
//                "Should either have (SpeciesIDs, SpeciesCoefficients, [SpeciesExponents]) or (Basins).\n"
//                "id: %d, speciesIDsExists: %d, speciesCoefficientsExists: %d, speciesExponentsExists: %d, basinsExists: %d", id, speciesIDsExists, speciesCoefficientsExists, speciesExponentsExists, basinsExists);
//    }

    return 0;
}

void Hdf5File::getOrderParameters(lm::types::OrderParameters* orderParameters) const
{
    // Make sure the orderParameters protobuf is not null and then clear it
    if (orderParameters == NULL) throw InvalidArgException("orderParameters", "cannot be null");
    orderParameters->Clear();

    // Declare and initialize the data structure for the callbacks in the order parameters iterator
    CallbackDataOrderParameters* opCD = new CallbackDataOrderParameters;
    opCD->orderParameters = orderParameters;

    if (H5Lexists(file, "/OrderParameters", H5P_DEFAULT) > 0)
    {
        H5Literate_by_name(file, "/OrderParameters", H5_INDEX_NAME, H5_ITER_INC, NULL, getOrderParametersCallback, (void *)opCD, H5P_DEFAULT);
    }
}

void Hdf5File::setOrderParameters(lm::types::OrderParameters * orderParameters)
{
    // Validate the set of order parameters
    if (orderParameters==NULL) throw InvalidArgException("orderParameters", "cannot be NULL");
    // TODO_LOW: add more checks

    // If a set of order parameters exist, delete it
    if (H5Lexists(file, "/OrderParameters", H5P_DEFAULT) > 0)
    {
        HDF5_EXCEPTION_CHECK(H5Ldelete(file, "/OrderParameters", H5P_DEFAULT));
    }

    hid_t opsGroup;
    HDF5_EXCEPTION_CALL(opsGroup, H5Gcreate2(file, "/OrderParameters", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    // If there are any order parameters, write out the relevant tables
    for (int i=0; i<orderParameters->order_parameter().size(); i++)
    {
        // Great the group for the order parameter.
        hid_t opGroup;
        char opName[64];
        snprintf(opName, 64, "%d", i);
        HDF5_EXCEPTION_CALL(opGroup, H5Gcreate2(opsGroup, opName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

        // Set the type.
        uint type = orderParameters->order_parameter(i).type();
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(opsGroup, opName, "Type", &type, 1));

        // write the order parameter's datasets
        hsize_t opDims[1];
        opDims[0] = hsize_t(orderParameters->order_parameter(i).species_index().size());
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(opGroup, "SpeciesIDs", 1, opDims, H5T_STD_U32LE, orderParameters->order_parameter(i).species_index().data()));
        opDims[0] = hsize_t(orderParameters->order_parameter(i).species_coefficient().size());
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(opGroup, "SpeciesCoefficients", 1, opDims, H5T_IEEE_F64LE, orderParameters->order_parameter(i).species_coefficient().data()));
        HDF5_EXCEPTION_CHECK(H5Gclose(opGroup));
    }
    HDF5_EXCEPTION_CHECK(H5Gclose(opsGroup));
}

bool Hdf5File::hasReactionModel() const
{
    return (H5Lexists(file, "/Model/Reaction", H5P_DEFAULT) > 0);
}

void Hdf5File::getReactionModel(lm::input::ReactionModel * reactionModel) const
{
    // Make sure the model is not null and then clear it.
    if (reactionModel == NULL) throw InvalidArgException("reactionModel", "cannot be null");
    reactionModel->Clear();

    if (H5Lexists(file, "/Model/Reaction", H5P_DEFAULT) > 0)
    {
        // Read at least the numbers of species.
        unsigned int numberSpecies;
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Reaction", "numberSpecies", &numberSpecies));
        reactionModel->set_number_species(numberSpecies);
        reactionModel->set_number_reactions(0);

        // If we have the number of reactions, we must have a full model so read it.
        if (H5Aexists_by_name(file, "/Model/Reaction", "numberReactions", H5P_DEFAULT) > 0)
        {
            hsize_t dims[3];
            H5T_class_t type;
            size_t size;

            // Read the initial species counts.
            H5LTget_dataset_info(file, "/Model/Reaction/InitialSpeciesCounts", dims, &type, &size);
            if (dims[0] != numberSpecies || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Reaction/InitialSpeciesCounts");
            int* intialSpeciesCountsBuffer = new int[numberSpecies];
            H5LTread_dataset_int(file, "/Model/Reaction/InitialSpeciesCounts", intialSpeciesCountsBuffer);
            for (uint i=0; i<numberSpecies; i++) reactionModel->add_initial_species_count(static_cast<uint>(intialSpeciesCountsBuffer[i]));
            delete[] intialSpeciesCountsBuffer;

            // Read the number of reactions.
            uint numberReactions;
            HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Reaction", "numberReactions", &numberReactions));
            reactionModel->set_number_reactions(numberReactions);

            // Read the reaction tables.
            if (numberReactions > 0)
            {
                // Make sure all of the data sets are the correct size.
                H5LTget_dataset_info(file, "/Model/Reaction/ReactionTypes", dims, &type, &size);
                if (dims[0] != numberReactions || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Reaction/ReactionTypes");
                H5LTget_dataset_info(file, "/Model/Reaction/ReactionRateConstants", dims, &type, &size);
                if (dims[0] != numberReactions || dims[1] != MAX_REACTION_RATE_CONSTANTS || size != sizeof(double)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Reaction/ReactionRateConstants");
                H5LTget_dataset_info(file, "/Model/Reaction/StoichiometricMatrix", dims, &type, &size);
                if (dims[0] != numberSpecies || dims[1] != numberReactions || size != sizeof(int)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Reaction/StoichiometricMatrix");
                H5LTget_dataset_info(file, "/Model/Reaction/DependencyMatrix", dims, &type, &size);
                if (dims[0] != numberSpecies || dims[1] != numberReactions || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Reaction/DependencyMatrix");

                // Allocate some buffers for reading the data.
                int * intBuffer = new int[numberSpecies*numberReactions];
                double * doubleBuffer = new double[numberReactions*MAX_REACTION_RATE_CONSTANTS];

                // Read the reaction info.
                H5LTread_dataset_int(file, "/Model/Reaction/ReactionTypes", intBuffer);
                H5LTread_dataset_double(file, "/Model/Reaction/ReactionRateConstants", doubleBuffer);
                for (uint i=0; i<numberReactions; i++)
                {
                    reactionModel->add_reaction();
                    reactionModel->mutable_reaction(i)->set_type((uint)intBuffer[i]);
                    for (uint j=0; j<MAX_REACTION_RATE_CONSTANTS; j++)
                    {
                        double k = doubleBuffer[i*MAX_REACTION_RATE_CONSTANTS+j];
                        if (!std::isnan(k))
                            reactionModel->mutable_reaction(i)->add_rate_constant(k);
                        else
                            break;
                    }
                }

                // Read the matrices.
                H5LTread_dataset_int(file, "/Model/Reaction/StoichiometricMatrix", intBuffer);
                for (uint i=0; i<numberSpecies*numberReactions; i++) reactionModel->add_stoichiometric_matrix(intBuffer[i]);
                H5LTread_dataset_int(file, "/Model/Reaction/DependencyMatrix", intBuffer);
                for (uint i=0; i<numberSpecies*numberReactions; i++) reactionModel->add_dependency_matrix((uint)intBuffer[i]);

                // Free the buffers.
                delete [] doubleBuffer;
                delete [] intBuffer;
            }

            // If we have a noise model, read it.
            if (H5Lexists(file, "/Model/Reaction/Noise", H5P_DEFAULT) > 0)
            {
                // Read the number of processes and update interval.
                uint32_t numberProcesses;
                double updateInterval;
                HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Reaction/Noise", "numberProcesses", &numberProcesses));
                HDF5_EXCEPTION_CHECK(H5LTget_attribute_double(file, "/Model/Reaction/Noise", "updateInterval", &updateInterval));
                reactionModel->mutable_noise_model()->set_number_processes(numberProcesses);
                reactionModel->mutable_noise_model()->set_process_update_interval(updateInterval);

                // Read the noise types matrx.
                H5LTget_dataset_info(file, "/Model/Reaction/Noise/NoiseTypes", dims, &type, &size);
                if (dims[0] != numberProcesses || size != sizeof(uint32_t)) THROW_EXCEPTION(lm::InvalidArgException, "invalid shape for noise_model.process_types");
                ndarray<uint32_t> T(utuple(static_cast<uint>(dims[0])));
                H5LTread_dataset_int(file, "/Model/Reaction/Noise/NoiseTypes", reinterpret_cast<int*>(T.values));
                NDArraySerializer::serializeInto(reactionModel->mutable_noise_model()->mutable_process_types(), T);

                // Write the noise parameters matrx.
                H5LTget_dataset_info(file, "/Model/Reaction/Noise/NoiseParameters", dims, &type, &size);
                if (dims[0] != numberProcesses || dims[1] == 0 || size != sizeof(double)) THROW_EXCEPTION(lm::InvalidArgException, "invalid shape for noise_model.process_parameters");
                ndarray<double> K(utuple(static_cast<uint>(dims[0]),static_cast<uint>(dims[1])));
                H5LTread_dataset_double(file, "/Model/Reaction/Noise/NoiseParameters", K.values);
                NDArraySerializer::serializeInto(reactionModel->mutable_noise_model()->mutable_process_parameters(), K);

                // Read the noise dependencies matrx.
                H5LTget_dataset_info(file, "/Model/Reaction/Noise/NoiseDependencies", dims, &type, &size);
                if (dims[0] != numberReactions || dims[1] == 0 || dims[2] != numberProcesses || size != sizeof(uint32_t)) THROW_EXCEPTION(lm::InvalidArgException, "invalid shape for noise_model.reaction_dependencies");
                ndarray<uint32_t> D(utuple(static_cast<uint>(dims[0]),static_cast<uint>(dims[1]),static_cast<uint>(dims[2])));
                H5LTread_dataset_int(file, "/Model/Reaction/Noise/NoiseDependencies", reinterpret_cast<int*>(D.values));
                NDArraySerializer::serializeInto(reactionModel->mutable_noise_model()->mutable_reaction_dependencies(), D);
            }
        }
    }
}

void Hdf5File::setReactionModel(lm::input::ReactionModel * reactionModel)
{
    // Validate that the model is consistent.
    if (reactionModel == NULL) throw InvalidArgException("reactionModel", "cannot be NULL");
    if (reactionModel->number_species() == 0) throw InvalidArgException("reactionModel.number_species", "cannot be zero");
    if (reactionModel->initial_species_count_size() != (int)reactionModel->number_species()) throw InvalidArgException("reactionModel.initial_species_count", "inconsistent size");
    if (reactionModel->reaction_size() != (int)reactionModel->number_reactions()) throw InvalidArgException("reactionModel.reaction", "inconsistent size");
    if (reactionModel->stoichiometric_matrix_size() != (int)(reactionModel->number_species()*reactionModel->number_reactions())) throw InvalidArgException("reactionModel.stoichiometric_matrix", "inconsistent size");
    if (reactionModel->dependency_matrix_size() != (int)(reactionModel->number_species()*reactionModel->number_reactions())) throw InvalidArgException("reactionModel.dependency_matrix", "inconsistent size");

    // If a reaction model already exists, delete it.
    if (H5Lexists(file, "/Model/Reaction", H5P_DEFAULT) > 0)
    {
        HDF5_EXCEPTION_CHECK(H5Ldelete(file, "/Model/Reaction", H5P_DEFAULT));
    }

    // Create the group for the reaction model.
    {
        hid_t group;
        HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Reaction", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Gclose(group));
    }

    // Write the numbers of species and reactions.
    numberSpecies = reactionModel->number_species();
    uint numberReactions = reactionModel->number_reactions();
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Reaction", "numberSpecies", &numberSpecies, 1));
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Reaction", "numberReactions", &numberReactions, 1));

    hsize_t dims[3];

    // Write the initial species counts.
    dims[0] = numberSpecies;
    HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/InitialSpeciesCounts", 1, dims, H5T_STD_U32LE, reactionModel->initial_species_count().data()));

    // Write the species names, if we have any.
    if (reactionModel->species_name_size() > 0)
    {
        // Make sure we have the correct number of names.
        if (reactionModel->species_name_size() != numberSpecies) throw InvalidArgException("reactionModel.species_name", "inconsistent size");

        // Copy the names into a fixed-size string buffer.
        char *nameData = new char[numberSpecies*MAX_NAME_LENGTH];
        for (int i=0; i<numberSpecies; i++)
            strncpy(&(nameData[i*MAX_NAME_LENGTH]), reactionModel->species_name(i).c_str(), MAX_NAME_LENGTH);
        dims[0] = numberSpecies;
        dims[1] = MAX_NAME_LENGTH;

        // Make the dataset.
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/SpeciesNames", 2, dims, H5T_C_S1, nameData));

        // Free the buffer.
        delete[] nameData;
    }

    // If we have any reactions, write out the reaction tables.
    if (reactionModel->number_reactions())
    {
        // Write the reaction tables.
        uint * types = new uint[numberReactions];
        double * constants = new double[numberReactions*MAX_REACTION_RATE_CONSTANTS];
        bool hasNameTable = false;
        char *nameData = new char[numberReactions*MAX_NAME_LENGTH];
        for (uint i=0; i<numberReactions; i++)
        {
            types[i] = reactionModel->reaction(i).type();
            uint j=0;
            for (; j<(uint)(reactionModel->reaction(i).rate_constant_size()) && j<MAX_REACTION_RATE_CONSTANTS; j++)
                constants[i*MAX_REACTION_RATE_CONSTANTS+j] = reactionModel->reaction(i).rate_constant(j);
            for (; j<MAX_REACTION_RATE_CONSTANTS; j++)
                constants[i*MAX_REACTION_RATE_CONSTANTS+j] = NAN;

            // If we have a name, fill in the row of the table.
            if (reactionModel->reaction(i).has_name())
            {
                hasNameTable = true;
                strncpy(&(nameData[i*MAX_NAME_LENGTH]), reactionModel->reaction(i).name().c_str(), MAX_NAME_LENGTH);
            }
        }
        dims[0] = numberReactions;
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/ReactionTypes", 1, dims, H5T_STD_U32LE, types));
        dims[0] = numberReactions;
        dims[1] = MAX_REACTION_RATE_CONSTANTS;
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/ReactionRateConstants", 2, dims, H5T_IEEE_F64LE, constants));

        // Write the reaction names, if we have any.
        if (hasNameTable)
        {
            dims[0] = numberReactions;
            dims[1] = MAX_NAME_LENGTH;
            HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/ReactionNames", 2, dims, H5T_C_S1, nameData));
        }

        delete[] nameData;
        delete[] constants;
        delete[] types;

        // Write the matrices.
        dims[0] = numberSpecies;
        dims[1] = numberReactions;
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/StoichiometricMatrix", 2, dims, H5T_STD_I32LE, reactionModel->stoichiometric_matrix().data()));
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/DependencyMatrix", 2, dims, H5T_STD_U32LE, reactionModel->dependency_matrix().data()));
    }

    // If there is a noise model, write it out.
    if (reactionModel->has_noise_model())
    {
        // Create the group for the noise model.
        hid_t group;
        HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Reaction/Noise", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Gclose(group));

        // Write the numbers of processes and update time.
        uint32_t numberProcesses = reactionModel->noise_model().number_processes();
        double updateInterval = reactionModel->noise_model().process_update_interval();
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Reaction/Noise", "numberProcesses", &numberProcesses, 1));
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_double(file, "/Model/Reaction/Noise", "updateInterval", &updateInterval, 1));

        // Write the noise types matrx.
        ndarray<uint32_t> T = NDArraySerializer::deserialize<uint32_t>(reactionModel->noise_model().process_types());
        if (T.shape.len != 1 || T.shape[0] != numberProcesses) THROW_EXCEPTION(lm::InvalidArgException, "invalid shape for noise_model.process_types");
        dims[0] = T.shape[0];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/Noise/NoiseTypes", 1, dims, H5T_STD_U32LE, T.values));

        // Write the noise types matrx.
        ndarray<double> K = NDArraySerializer::deserialize<double>(reactionModel->noise_model().process_parameters());
        if (K.shape.len != 2 || K.shape[0] != numberProcesses || K.shape[1] == 0) THROW_EXCEPTION(lm::InvalidArgException, "invalid shape for noise_model.process_parameters");
        dims[0] = K.shape[0];
        dims[1] = K.shape[1];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/Noise/NoiseParameters", 2, dims, H5T_IEEE_F64LE, K.values));

        // Write the noise types matrx.
        ndarray<uint32_t> D = NDArraySerializer::deserialize<uint32_t>(reactionModel->noise_model().reaction_dependencies());
        if (D.shape.len != 3 || D.shape[0] != numberReactions || D.shape[1] == 0 || D.shape[2] != numberProcesses) THROW_EXCEPTION(lm::InvalidArgException, "invalid shape for noise_model.reaction_dependencies");
        dims[0] = D.shape[0];
        dims[1] = D.shape[1];
        dims[2] = D.shape[2];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Reaction/Noise/NoiseDependencies", 3, dims, H5T_STD_U32LE, D.values));
    }
}

void Hdf5File::getSpatialModel(lm::input::SpatialModel * spatialModel) const
{
    // Make sure the model is not null and then clear it.
    if (spatialModel == NULL) throw InvalidArgException("spatialModel", "cannot be null");
    spatialModel->Clear();

    if (H5Lexists(file, "/Model/Spatial", H5P_DEFAULT) > 0)
    {
        // Read the diffusion model attributes.
        uint numberRegions, numberObstacles;
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Spatial/Regions", "numberRegions", &numberRegions));
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_uint(file, "/Model/Spatial/Obstacles", "numberObstacles", &numberObstacles));

        hsize_t dims[2];
        H5T_class_t type;
        size_t size;

        if (numberRegions > 0)
        {
            // Read the region types table.
            {
                H5LTget_dataset_info(file, "/Model/Spatial/Regions/Types", dims, &type, &size);
                if (dims[0] != numberRegions || dims[1] != 2 || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Spatial/Regions/Types");
                uint * types = new uint[dims[0]*dims[1]];
                H5LTread_dataset(file, "/Model/Spatial/Regions/Types", H5T_STD_U32LE, types);
                for (uint i=0; i<dims[0]; i++)
                {
                    spatialModel->add_region();
                    spatialModel->mutable_region(i)->set_shape(types[i*dims[1]]);
                    spatialModel->mutable_region(i)->set_site_type(types[i*dims[1]+1]);
                }
                delete [] types;
            }

            // Read the shape parameters table.
            {
                H5LTget_dataset_info(file, "/Model/Spatial/Regions/ShapeParameters", dims, &type, &size);
                if (dims[0] != numberRegions || dims[1] != MAX_SHAPE_PARAMETERS || size != sizeof(double)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Spatial/Regions/ShapeParameters");
                double * shapeParams = new double[dims[0]*dims[1]];
                H5LTread_dataset_double(file, "/Model/Spatial/Regions/ShapeParameters", shapeParams);
                for (uint i=0; i<dims[0]; i++)
                {
                    for (uint j=0; j<dims[1]; j++)
                    {
                        double shapeParam = shapeParams[i*dims[1]+j];
                        if (!std::isnan(shapeParam))
                            spatialModel->mutable_region(i)->add_shape_parameter(shapeParam);
                        else
                            break;
                    }
                }
                delete [] shapeParams;
            }
        }

        if (numberObstacles > 0)
        {
            // Read the obstacle types table.
            {
                H5LTget_dataset_info(file, "/Model/Spatial/Obstacles/Types", dims, &type, &size);
                if (dims[0] != numberObstacles || dims[1] != 2 || size != sizeof(uint)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Spatial/Obstacles/Types");
                uint * types = new uint[dims[0]*dims[1]];
                H5LTread_dataset(file, "/Model/Spatial/Obstacles/Types", H5T_STD_U32LE, types);
                for (uint i=0; i<dims[0]; i++)
                {
                    spatialModel->add_obstacle();
                    spatialModel->mutable_obstacle(i)->set_shape(types[i*dims[1]]);
                    spatialModel->mutable_obstacle(i)->set_site_type(types[i*dims[1]+1]);
                }
                delete [] types;
            }

            // Read the obstacle parameters table.
            {
                H5LTget_dataset_info(file, "/Model/Spatial/Obstacles/ShapeParameters", dims, &type, &size);
                if (dims[0] != numberObstacles || dims[1] != MAX_SHAPE_PARAMETERS || size != sizeof(double)) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Spatial/Obstacles/ShapeParameters");
                double * shapeParams = new double[dims[0]*dims[1]];
                H5LTread_dataset_double(file, "/Model/Spatial/Obstacles/ShapeParameters", shapeParams);
                for (uint i=0; i<dims[0]; i++)
                {
                    for (uint j=0; j<dims[1]; j++)
                    {
                        double shapeParam = shapeParams[i*dims[1]+j];
                        if (!std::isnan(shapeParam))
                            spatialModel->mutable_obstacle(i)->add_shape_parameter(shapeParam);
                        else
                            break;
                    }
                }
                delete [] shapeParams;
            }
        }
    }
}

void Hdf5File::setSpatialModel(lm::input::SpatialModel * spatialModel)
{
    // Validate that the model is consistent.
    if (spatialModel == NULL) throw InvalidArgException("spatialModel", "cannot be NULL");

    // If a diffusion model already exists, delete it.
    if (H5Lexists(file, "/Model/Spatial", H5P_DEFAULT) > 0)
    {
        HDF5_EXCEPTION_CHECK(H5Ldelete(file, "/Model/Spatial", H5P_DEFAULT));
    }

    // Create the group for the spatial model.
    hid_t group;
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Spatial", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Spatial/Regions", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));
    HDF5_EXCEPTION_CALL(group,H5Gcreate2(file, "/Model/Spatial/Obstacles", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    HDF5_EXCEPTION_CHECK(H5Gclose(group));

    // Write the attributes.
    uint numberRegions = spatialModel->region_size();
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Spatial/Regions", "numberRegions", &numberRegions, 1));
    uint numberObstacles = spatialModel->obstacle_size();
    HDF5_EXCEPTION_CHECK(H5LTset_attribute_uint(file, "/Model/Spatial/Obstacles", "numberObstacles", &numberObstacles, 1));

    if (numberRegions > 0)
    {
        // Write the region types table.
        {
            const unsigned int RANK=2;
            hsize_t dims[RANK];
            dims[0] = numberRegions;
            dims[1] = 2;
            uint * types = new uint[dims[0]*dims[1]];
            for (uint i=0; i<dims[0]; i++)
            {
                types[i*dims[1]] = spatialModel->region(i).shape();
                types[i*dims[1]+1] = spatialModel->region(i).site_type();
            }
            HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Spatial/Regions/Types", RANK, dims, H5T_STD_U32LE, types));
            delete [] types;
        }

        // Write the shape parameters table.
        {
            const unsigned int RANK=2;
            hsize_t dims[RANK];
            dims[0] = numberRegions;
            dims[1] = MAX_SHAPE_PARAMETERS;
            double * shapeParams = new double[dims[0]*dims[1]];
            for (uint i=0; i<dims[0]; i++)
            {
                for (uint j=0; j<dims[1]; j++)
                {
                    if (j < (uint)spatialModel->region(i).shape_parameter_size())
                        shapeParams[i*dims[1]+j] = spatialModel->region(i).shape_parameter(j);
                    else
                        shapeParams[i*dims[1]+j] = NAN;
                }
            }
            HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Spatial/Regions/ShapeParameters", RANK, dims, H5T_IEEE_F64LE, shapeParams));
            delete [] shapeParams;
        }
    }

    if (numberObstacles > 0)
    {
        // Write the obstacle types table.
        {
            const unsigned int RANK=2;
            hsize_t dims[RANK];
            dims[0] = numberObstacles;
            dims[1] = 2;
            uint * types = new uint[dims[0]*dims[1]];
            for (uint i=0; i<dims[0]; i++)
            {
                types[i*dims[1]] = spatialModel->obstacle(i).shape();
                types[i*dims[1]+1] = spatialModel->obstacle(i).site_type();
            }
            HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Spatial/Obstacles/Types", RANK, dims, H5T_STD_U32LE, types));
            delete [] types;
        }

        // Write the shape parameters table.
        {
            const unsigned int RANK=2;
            hsize_t dims[RANK];
            dims[0] = numberObstacles;
            dims[1] = MAX_SHAPE_PARAMETERS;
            double * shapeParams = new double[dims[0]*dims[1]];
            for (uint i=0; i<dims[0]; i++)
            {
                for (uint j=0; j<dims[1]; j++)
                {
                    if (j < (uint)spatialModel->obstacle(i).shape_parameter_size())
                        shapeParams[i*dims[1]+j] = spatialModel->obstacle(i).shape_parameter(j);
                    else
                        shapeParams[i*dims[1]+j] = NAN;
                }
            }
            HDF5_EXCEPTION_CHECK(H5LTmake_dataset(file, "/Model/Spatial/Obstacles/ShapeParameters", RANK, dims, H5T_IEEE_F64LE, shapeParams));
            delete [] shapeParams;
        }
    }
}

bool Hdf5File::hasTilings() const
{
    return (H5Lexists(file, "/Tilings", H5P_DEFAULT) > 0);
}

herr_t Hdf5File::getTilingsCallback(hid_t loc_id, const char* name, const H5L_info_t* info, void* callbackDataTilings)
{
    // Declare and initialize handle for the tiling group
    hid_t tilingGroup;
    HDF5_EXCEPTION_CALL(tilingGroup, H5Gopen(loc_id, name, H5P_DEFAULT));

    // recast the callbackData structure away from void *
    CallbackDataTilings* cdT = reinterpret_cast<CallbackDataTilings *>(callbackDataTilings);

    // create a new message in the tilings protobuf.
    lm::types::Tiling* newTiling = cdT->tilings->add_tiling();

    // Make sure the data we need exists.
    if (H5Aexists(tilingGroup, "OrderParameterIndex") > 0 && H5Lexists(tilingGroup, "Edges", H5P_DEFAULT))
    {
        // Get the order parameter associated with this tiling.
        int orderParametersIndex;
        HDF5_EXCEPTION_CHECK(H5LTget_attribute_int(loc_id, name, "OrderParameterIndex", &orderParametersIndex));
        ndarray<int32_t> orderParametersIndices(1);
        orderParametersIndices[0] = orderParametersIndex;
        NDArraySerializer::serializeInto(newTiling->mutable_order_parameter_indices(), orderParametersIndices);

        // Read in the values of the tiling's edges.
        hsize_t dims[1];
        H5T_class_t hdf5Type;
        size_t size;
        HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(tilingGroup, "Edges", dims, &hdf5Type, &size));
        ndarray<double> edges(static_cast<uint>(dims[0]));
        HDF5_EXCEPTION_CHECK(H5LTread_dataset_double(tilingGroup, "Edges", edges.values));
        NDArraySerializer::serializeInto(newTiling->mutable_edges(), edges);
    }

    // Free the group handle
    HDF5_EXCEPTION_CHECK(H5Gclose(tilingGroup));

    return 0;
}

void Hdf5File::getTilings(lm::types::Tilings* tilings) const
{
    // Make sure the tilings protobuf is not null and then clear it
    if (tilings == NULL) throw InvalidArgException("tilings", "cannot be null");
    tilings->Clear();

    // Declare and initialize the data structure for the callbacks in the tilings iterator
    CallbackDataTilings* cdT = new CallbackDataTilings;
    cdT->tilings = tilings;

    if (H5Lexists(file, "/Tilings", H5P_DEFAULT) > 0)
    {
        H5Literate_by_name(file, "/Tilings", H5_INDEX_NAME, H5_ITER_INC, NULL, getTilingsCallback, (void *)cdT, H5P_DEFAULT);
    }
}

void Hdf5File::setTilings(lm::types::Tilings* tilings)
{
    // Validate the set of tilings
    if (tilings == NULL) throw InvalidArgException("tilings", "cannot be NULL");
    // TODO_LOW: add more checks

    // If a set of tilings exist, delete it
    if (H5Lexists(file, "/Tilings", H5P_DEFAULT) > 0)
    {
        HDF5_EXCEPTION_CHECK(H5Ldelete(file, "/Tilings", H5P_DEFAULT));
    }

    hid_t tilingsGroup;
    HDF5_EXCEPTION_CALL(tilingsGroup, H5Gcreate2(file, "/Tilings", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    // If there are any sets of tilings, write out the relevant tables
    for (int i=0; i<tilings->tiling_size(); i++)
    {
        hid_t tilingGroup;
        char tilingName[64];
        snprintf(tilingName, 64, "%d", i);
        HDF5_EXCEPTION_CALL(tilingGroup, H5Gcreate2(tilingsGroup, tilingName, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

        if (tilings->tiling(i).has_order_parameter_indices() && tilings->tiling(i).has_edges())
        {
            // Write the tiling's order parameter.
            ndarray<int32_t> orderParameterIndices = NDArraySerializer::deserialize<int32_t>(tilings->tiling(i).order_parameter_indices());
            if (orderParameterIndices.shape.len != 1) THROW_EXCEPTION(lm::RuntimeException, "the tiling.order_parameter_indices array was not 1D: %d", orderParameterIndices.shape.len);
            if (orderParameterIndices.shape[0] != 1) THROW_EXCEPTION(lm::RuntimeException, "only single dimensional tilings are currently supported");
            int32_t orderParameterIndex = orderParameterIndices[0];
            HDF5_EXCEPTION_CHECK(H5LTset_attribute_int(tilingsGroup, tilingName, "OrderParameterIndex", &orderParameterIndex, 1));

            // Write the tiling's edges.
            ndarray<double> edges = NDArraySerializer::deserialize<double>(tilings->tiling(i).edges());
            if (edges.shape.len != orderParameterIndices.shape[0]) THROW_EXCEPTION(lm::RuntimeException, "inconsistent tiling data: %d %d", edges.shape.len, orderParameterIndices.shape[0]);
            if (edges.shape.len != 1) THROW_EXCEPTION(lm::RuntimeException, "only single dimensional tilings are currently supported");
            hsize_t edgesDims[1];
            edgesDims[0] = edges.shape[0];
            HDF5_EXCEPTION_CHECK(H5LTmake_dataset(tilingGroup, "Edges", 1, edgesDims, H5T_IEEE_F64LE, edges.values));
        }

        // Close the group.
        HDF5_EXCEPTION_CHECK(H5Gclose(tilingGroup));
    }
    HDF5_EXCEPTION_CHECK(H5Gclose(tilingsGroup));
}

bool Hdf5File::hasFFPilotOptions() const
{
    return (H5Lexists(file, "/FFPilot/Input", H5P_DEFAULT) > 0);
}

void Hdf5File::getFFPilotOptions(lm::input::ffpilot::FFPilotOptions* options) const
{
    // Make sure the tilings protobuf is not null and then clear it
    if (options == NULL) throw InvalidArgException("options", "cannot be null");
    options->Clear();

    if (H5Lexists(file, "/FFPilot/Input", H5P_DEFAULT) > 0)
    {
        // Open the group.
        hid_t optionsGroup;
        HDF5_EXCEPTION_CALL(optionsGroup,H5Gopen2(file, "/FFPilot/Input", H5P_DEFAULT));

        // Read in the options.
        if (H5Lexists(optionsGroup, "Basins", H5P_DEFAULT) > 0)
        {
            int ndims;
            HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(optionsGroup, "Basins", &ndims));
            hsize_t* dims = new hsize_t[ndims];
            H5T_class_t hdf5Type;
            size_t size;
            HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(optionsGroup, "Basins", dims, &hdf5Type, &size));
            if (ndims != 2) THROW_EXCEPTION(lm::RuntimeException, "the Basins array was not 2D: %d", ndims);
            if (dims[0] != 2) THROW_EXCEPTION(lm::RuntimeException, "the Basins array did not have two rows: %d", dims[0]);
            ndarray<int32_t> counts(utuple(static_cast<uint>(dims[0]),static_cast<uint>(dims[1])));
            HDF5_EXCEPTION_CHECK(H5LTread_dataset_int(optionsGroup, "Basins", counts.values));
            NDArraySerializer::serializeInto(options->mutable_phase_zero_basins(), counts);
            delete[] dims;
        }
        if (H5Aexists_by_name(file, "/FFPilot/Input", "errorGoal", H5P_DEFAULT) > 0)
        {
            double value;
            HDF5_EXCEPTION_CHECK(H5LTget_attribute_double(file, "/FFPilot/Input", "errorGoal", &value));
            options->set_error_goal(value);
        }
        if (H5Aexists_by_name(file, "/FFPilot/Input", "fallbackMethod", H5P_DEFAULT) > 0)
        {
            int value;
            HDF5_EXCEPTION_CHECK(H5LTget_attribute_int(file, "/FFPilot/Input", "fallbackMethod", &value));
            options->set_fallback_method(static_cast<lm::input::ffpilot::FFPilotOptions_FallbackMethod>(value));
        }
        if (H5Aexists_by_name(file, "/FFPilot/Input", "pilotSkip", H5P_DEFAULT) > 0)
        {
            int value;
            HDF5_EXCEPTION_CHECK(H5LTget_attribute_int(file, "/FFPilot/Input", "pilotSkip", &value));
            options->set_pilot_skip(value!=0);
        }
        if (H5Aexists_by_name(file, "/FFPilot/Input", "pilotStageCrossings", H5P_DEFAULT) > 0)
        {
            int value;
            HDF5_EXCEPTION_CHECK(H5LTget_attribute_int(file, "/FFPilot/Input", "pilotStageCrossings", &value));
            options->set_pilot_stage_crossings(value);
        }
        if (H5Aexists_by_name(file, "/FFPilot/Input", "prodSkip", H5P_DEFAULT) > 0)
        {
            int value;
            HDF5_EXCEPTION_CHECK(H5LTget_attribute_int(file, "/FFPilot/Input", "prodSkip", &value));
            options->set_prod_skip(value!=0);
        }
        if (H5Lexists(optionsGroup, "ProdTrajectoryCounts", H5P_DEFAULT) > 0)
        {
            int ndims;
            HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(optionsGroup, "ProdTrajectoryCounts", &ndims));
            hsize_t* dims = new hsize_t[ndims];
            H5T_class_t hdf5Type;
            size_t size;
            HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(optionsGroup, "ProdTrajectoryCounts", dims, &hdf5Type, &size));
            if (ndims != 2) THROW_EXCEPTION(lm::RuntimeException, "the ProdTrajectoryCounts array was not 2D: %d", ndims);
            if (dims[0] != 2) THROW_EXCEPTION(lm::RuntimeException, "the ProdTrajectoryCounts array did not have two rows: %d", dims[0]);
            ndarray<uint64_t> counts(utuple(static_cast<uint>(dims[0]),static_cast<uint>(dims[1])));
            HDF5_EXCEPTION_CHECK(H5LTread_dataset_long(optionsGroup, "ProdTrajectoryCounts", reinterpret_cast<long*>(counts.values)));
            NDArraySerializer::serializeInto(options->mutable_prod_trajectory_counts(), counts);
            delete[] dims;
        }
        if (H5Lexists(optionsGroup, "SamplingMultipliers", H5P_DEFAULT) > 0)
        {
            int ndims;
            HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(optionsGroup, "SamplingMultipliers", &ndims));
            hsize_t* dims = new hsize_t[ndims];
            H5T_class_t hdf5Type;
            size_t size;
            HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(optionsGroup, "SamplingMultipliers", dims, &hdf5Type, &size));
            if (ndims != 2) THROW_EXCEPTION(lm::RuntimeException, "the SamplingMultipliers array was not 2D: %d", ndims);
            if (dims[0] != 2) THROW_EXCEPTION(lm::RuntimeException, "the SamplingMultipliers array did not have two rows: %d", dims[0]);
            ndarray<double> multipliers(utuple(static_cast<uint>(dims[0]),static_cast<uint>(dims[1])));
            HDF5_EXCEPTION_CHECK(H5LTread_dataset_double(optionsGroup, "SamplingMultipliers", reinterpret_cast<double*>(multipliers.values)));
            NDArraySerializer::serializeInto(options->mutable_optimize_sampling_multipliers(), multipliers);
            delete[] dims;
        }
        if (H5Aexists_by_name(file, "/FFPilot/Input", "trajectoryOutput", H5P_DEFAULT) > 0)
        {
            int value;
            HDF5_EXCEPTION_CHECK(H5LTget_attribute_int(file, "/FFPilot/Input", "trajectoryOutput", &value));
            options->set_trajectory_output(value!=0);
        }

        // Close the group.
        HDF5_EXCEPTION_CHECK(H5Gclose(optionsGroup));
    }
}

void Hdf5File::setFFPilotOptions(lm::input::ffpilot::FFPilotOptions* options)
{
    // Validate the set of tilings
    if (options == NULL) throw InvalidArgException("options", "cannot be NULL");

    // Make sure we have the ffpilot top level group.
    if (H5Lexists(file, "/FFPilot", H5P_DEFAULT) == 0) HDF5_EXCEPTION_CHECK(H5Gcreate2(file, "/FFPilot", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT))

    // If a set of options exist, delete it
    if (H5Lexists(file, "/FFPilot/Input", H5P_DEFAULT) > 0)
    {
        HDF5_EXCEPTION_CHECK(H5Ldelete(file, "/FFPilot/Input", H5P_DEFAULT));
    }

    // Create the group.
    hid_t optionsGroup;
    HDF5_EXCEPTION_CALL(optionsGroup, H5Gcreate2(file, "/FFPilot/Input", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    // Add the options.
    if (options->has_phase_zero_basins())
    {
        // Write the tiling's order parameter.
        ndarray<int32_t> counts = NDArraySerializer::deserialize<int32_t>(options->phase_zero_basins());
        if (counts.shape.len != 2) THROW_EXCEPTION(lm::RuntimeException, "the phase_zero_basins array was not 2D: %d", counts.shape.len);
        if (counts.shape[0] != 2) THROW_EXCEPTION(lm::RuntimeException, "the phase_zero_basins array did not have two rows: %d", counts.shape[0]);
        const unsigned int RANK=2;
        hsize_t dims[RANK];
        dims[0] = counts.shape[0];
        dims[1] = counts.shape[1];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(optionsGroup, "Basins", RANK, dims, H5T_STD_I32LE, counts.values));
    }
    if (options->has_error_goal())
    {
        double value = options->error_goal();
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_double(file, "/FFPilot/Input", "errorGoal", &value, 1));
    }
    if (options->has_fallback_method())
    {
        int value = options->fallback_method();
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_int(file, "/FFPilot/Input", "fallbackMethod", &value, 1));
    }
    if (options->has_pilot_skip())
    {
        int value = options->pilot_skip()?1:0;
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_int(file, "/FFPilot/Input", "pilotSkip", &value, 1));
    }
    if (options->has_pilot_stage_crossings())
    {
        int value = options->pilot_stage_crossings();
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_int(file, "/FFPilot/Input", "pilotStageCrossings", &value, 1));
    }
    if (options->has_prod_skip())
    {
        int value = options->prod_skip()?1:0;
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_int(file, "/FFPilot/Input", "prodSkip", &value, 1));
    }
    if (options->has_prod_trajectory_counts())
    {
        // Write the tiling's order parameter.
        ndarray<uint64_t> counts = NDArraySerializer::deserialize<uint64_t>(options->prod_trajectory_counts());
        if (counts.shape.len != 2) THROW_EXCEPTION(lm::RuntimeException, "the prod_trajectory_counts array was not 2D: %d", counts.shape.len);
        if (counts.shape[0] != 2) THROW_EXCEPTION(lm::RuntimeException, "the prod_trajectory_counts array did not have two rows: %d", counts.shape[0]);
        const unsigned int RANK=2;
        hsize_t dims[RANK];
        dims[0] = counts.shape[0];
        dims[1] = counts.shape[1];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(optionsGroup, "ProdTrajectoryCounts", RANK, dims, H5T_STD_U64LE, counts.values));
    }
    if (options->has_optimize_sampling_multipliers())
    {
        // Write the samplnig multipliers.
        ndarray<double> multipliers = NDArraySerializer::deserialize<double>(options->optimize_sampling_multipliers());
        if (multipliers.shape.len != 2) THROW_EXCEPTION(lm::RuntimeException, "the optimize_sampling_multipliers array was not 2D: %d", multipliers.shape.len);
        if (multipliers.shape[0] != 2) THROW_EXCEPTION(lm::RuntimeException, "the optimize_sampling_multipliers array did not have two rows: %d", multipliers.shape[0]);
        const unsigned int RANK=2;
        hsize_t dims[RANK];
        dims[0] = multipliers.shape[0];
        dims[1] = multipliers.shape[1];
        HDF5_EXCEPTION_CHECK(H5LTmake_dataset(optionsGroup, "SamplingMultipliers", RANK, dims, H5T_IEEE_F64LE, multipliers.values));
    }

    if (options->has_trajectory_output())
    {
        int value = options->trajectory_output()?1:0;
        HDF5_EXCEPTION_CHECK(H5LTset_attribute_int(file, "/FFPilot/Input", "trajectoryOutput", &value, 1));
    }

    // Close the group.
    HDF5_EXCEPTION_CHECK(H5Gclose(optionsGroup));
}


bool Hdf5File::hasReactionRestart() const
{
    Print::printf(Print::INFO, "Restart reaction check %d",H5Lexists(file, "/Restart/Reaction", H5P_DEFAULT));
    return (H5Lexists(file, "/Restart/Reaction", H5P_DEFAULT) > 0);
}

void Hdf5File::getReactionRestart(lm::input::CMERestart* restart) const
{
    // Make sure the tilings protobuf is not null and then clear it
    if (restart == NULL) throw InvalidArgException("restart", "cannot be null");
    restart->Clear();

    Print::printf(Print::INFO, "Restart reaction %d %d",H5Lexists(file, "/Restart/Reaction", H5P_DEFAULT),H5Lexists(file, "/Restart/Reaction/Trajectories", H5P_DEFAULT));

    if (H5Lexists(file, "/Restart/Reaction/Trajectories", H5P_DEFAULT) > 0)
    {
        int ndims;
        hsize_t* dims;
        H5T_class_t hdf5Type;
        size_t size;
        HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(file, "/Restart/Reaction/Trajectories", &ndims));
        dims = new hsize_t[ndims];
        HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(file, "/Restart/Reaction/Trajectories", dims, &hdf5Type, &size));
        if (ndims != 1 || hdf5Type == H5T_STD_U64LE) THROW_EXCEPTION(lm::InputException, "Invalid properties for array /Restart/Reaction/Trajectories: %d %d", ndims, hdf5Type);
        ndarray<uint64_t> trajectories(utuple(static_cast<uint>(dims[0])));
//        HDF5_EXCEPTION_CHECK(H5LTread_dataset_long(file, "/Restart/Reaction/Trajectories", reinterpret_cast<int64_t*>(trajectories.values)));
		HDF5_EXCEPTION_CHECK(H5LTread_dataset_long(file, "/Restart/Reaction/Trajectories", reinterpret_cast<long*>(trajectories.values)));
        NDArraySerializer::serializeInto(restart->mutable_restart_trajectory_ids(), trajectories);
        delete[] dims;
    }
    if (H5Lexists(file, "/Restart/Reaction/Times", H5P_DEFAULT) > 0)
    {
        int ndims;
        hsize_t* dims;
        H5T_class_t hdf5Type;
        size_t size;
        HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(file, "/Restart/Reaction/Times", &ndims));
        dims = new hsize_t[ndims];
        HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(file, "/Restart/Reaction/Times", dims, &hdf5Type, &size));
        if (ndims != 1 || hdf5Type == H5T_IEEE_F64LE) THROW_EXCEPTION(lm::InputException, "Invalid properties for array /Restart/Reaction/Times: %d %d", ndims, hdf5Type);
        ndarray<double> times(utuple(static_cast<uint>(dims[0])));
        HDF5_EXCEPTION_CHECK(H5LTread_dataset_double(file, "/Restart/Reaction/Times", reinterpret_cast<double*>(times.values)));
        NDArraySerializer::serializeInto(restart->mutable_restart_times(), times);
        delete[] dims;
    }
    if (H5Lexists(file, "/Restart/Reaction/SpeciesCounts", H5P_DEFAULT) > 0)
    {
        int ndims;
        hsize_t* dims;
        H5T_class_t hdf5Type;
        size_t size;
        HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(file, "/Restart/Reaction/SpeciesCounts", &ndims));
        dims = new hsize_t[ndims];
        HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(file, "/Restart/Reaction/SpeciesCounts", dims, &hdf5Type, &size));
        if (ndims != 2 || hdf5Type == H5T_STD_I32LE) THROW_EXCEPTION(lm::InputException, "Invalid properties for array /Restart/Reaction/SpeciesCounts: %d %d", ndims, hdf5Type);
        ndarray<int32_t> counts(utuple(static_cast<uint>(dims[0]),static_cast<uint>(dims[1])));
        HDF5_EXCEPTION_CHECK(H5LTread_dataset_int(file, "/Restart/Reaction/SpeciesCounts", reinterpret_cast<int32_t*>(counts.values)));
        NDArraySerializer::serializeInto(restart->mutable_restart_species_counts(), counts);
        delete[] dims;
    }
}

bool Hdf5File::hasDiffusionRestart() const
{
    Print::printf(Print::INFO, "Restart diffusion check %d",H5Lexists(file, "/Restart/Diffusion", H5P_DEFAULT));
    return (H5Lexists(file, "/Restart/Diffusion", H5P_DEFAULT) > 0);
}

herr_t Hdf5File::getDiffusionRestartCallback(hid_t loc, const char* name, const H5L_info_t* info, void* obj)
{
    Print::printf(Print::INFO, "Restart lattice :%s:", name);

    // Recast the restart object.
    lm::input::RDMERestart* restart = static_cast<lm::input::RDMERestart *>(obj);

    // Get the lattice index from the name.
    int index = atoi(name);

    // Make sure the index is valid.
    if (index < 0 || index >= restart->restart_lattice_size()) THROW_EXCEPTION(lm::InputException, "Invalid restart lattice index: %d (%s), must be between 0 and %d", index, name, index >= restart->restart_lattice_size());

    // Read the lattice.
    int ndims;
    hsize_t* dims;
    H5T_class_t hdf5Type;
    size_t size;
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(loc, name, &ndims));
    dims = new hsize_t[ndims];
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(loc, name, dims, &hdf5Type, &size));
    if (ndims != 4 || hdf5Type == H5T_STD_U8LE) THROW_EXCEPTION(lm::InputException, "Invalid properties for array /Restart/Diffusion/Lattice/%s: %d %d", name, ndims, hdf5Type);
    ndarray<uint8_t> lattice(utuple(static_cast<uint>(dims[0]),static_cast<uint>(dims[1]),static_cast<uint>(dims[2]),static_cast<uint>(dims[3])));
    HDF5_EXCEPTION_CHECK(H5LTread_dataset_uint8(loc, name, lattice.values));
    delete[] dims;

    // Add the lattice data to the lattice.
    NDArraySerializer::serializeInto(restart->mutable_restart_lattice(index)->mutable_particles(), lattice);

    return 0;
}

void Hdf5File::getDiffusionRestart(lm::input::RDMERestart* restart, uint numberTrajectories) const
{
    // Make sure the tilings protobuf is not null and then clear it
    if (restart == NULL) throw InvalidArgException("tilings", "cannot be null");
    restart->Clear();

    if (H5Lexists(file, "/Restart/Diffusion/Lattice", H5P_DEFAULT) > 0)
    {
        // Create an empty lattice for each trajectory.
        for (uint i=0; i<numberTrajectories; i++)
            restart->add_restart_lattice();

        // Iterate through the datasets and create the lattices.
        H5Literate_by_name(file, "/Restart/Diffusion/Lattice", H5_INDEX_NAME, H5_ITER_INC, NULL, getDiffusionRestartCallback, (void *)restart, H5P_DEFAULT);
    }
}


bool Hdf5File::replicateExists(uint64_t replicate)
{
    char replicateName[8];
    snprintf(replicateName, sizeof(replicateName), "%07d", (int)replicate);
    if (H5Lexists(simulationsGroup, replicateName, H5P_DEFAULT) > 0) return true;
    return false;
}

void Hdf5File::openReplicate(uint64_t replicate)
{
    openReplicateHandles(replicate);
}

void Hdf5File::appendSpeciesCounts(uint64_t replicate, lm::io::SpeciesCounts * speciesCounts)
{
    appendSpeciesTimeSeries(replicate, speciesCounts->number_entries(), speciesCounts->number_species(), speciesCounts->species_count().data(), speciesCounts->time().data());
}

void Hdf5File::appendSpeciesTimeSeries(uint64_t replicate, const lm::io::SpeciesTimeSeries& speciesTimeSeries)
{
    int numberEntries = speciesTimeSeries.counts().shape(0);
    int numberSpecies = speciesTimeSeries.counts().shape(1);

    if (speciesTimeSeries.times().shape(0) != numberEntries)
        InvalidArgException("speciesTimeSeries.times.shape", "Numebr of rows in time array incocnsistent with counts array.");

    // Extract the data.
    ndarray<int32_t> *counts = NDArraySerializer::deserializeAllocate<int32_t>(speciesTimeSeries.counts());
    ndarray<double> *times = NDArraySerializer::deserializeAllocate<double>(speciesTimeSeries.times());

    // Append  the data.
    appendSpeciesTimeSeries(replicate, numberEntries, numberSpecies, counts->values, times->values);

    // Free any allocated memory.
    delete counts;
    delete times;
}

void Hdf5File::appendSpeciesTimeSeries(uint64_t replicate, int numberEntries, int numberSpecies, const int32_t* counts, const double* times)
{
    ReplicateHandles * handles = openReplicateHandles(replicate);

    // Update the species counts dataset.
    {
        // Get the current size of the dataset.
        unsigned int RANK=2;
        hsize_t dims[RANK];
        hid_t dataspace_id;
        int result;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(handles->speciesCountsDataset));
        HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspace_id, dims, NULL));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));

        // Extend the dataset by the number of rows in the data set.
        dims[0] += numberEntries;
        HDF5_EXCEPTION_CHECK(H5Dset_extent(handles->speciesCountsDataset, dims));

        // Create the memory dataset.
        hid_t memspace_id;
        hsize_t memDims[RANK];
        memDims[0] = numberEntries;
        memDims[1] = numberSpecies;
        HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));

        // Write the new data.
        HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(handles->speciesCountsDataset));
        hsize_t start[RANK], count[RANK];
        start[0] = dims[0]-numberEntries;
        start[1] = 0;
        count[0] = memDims[0];
        count[1] = memDims[1];
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
        HDF5_EXCEPTION_CHECK(H5Dwrite(handles->speciesCountsDataset, H5T_NATIVE_INT32, memspace_id, dataspace_id, H5P_DEFAULT, counts));

        // Cleanup some resources.
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
    }

    // Update the species count times dataset.
    {
        // Get the current size of the dataset.
        unsigned int RANK=1;
        hsize_t dims[RANK];
        hid_t dataspace_id;
        int result;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(handles->speciesCountTimesDataset));
        HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspace_id, dims, NULL));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));

        // Extend the dataset by the number of rows in the data set.
        dims[0] += numberEntries;
        HDF5_EXCEPTION_CHECK(H5Dset_extent(handles->speciesCountTimesDataset, dims));

        // Create the memory dataset.
        hid_t memspace_id;
        hsize_t memDims[RANK];
        memDims[0] = numberEntries;
        HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));

        // Write the new data.
        HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(handles->speciesCountTimesDataset));
        hsize_t start[RANK], count[RANK];
        start[0] = dims[0]-numberEntries;
        count[0] = memDims[0];
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
        HDF5_EXCEPTION_CHECK(H5Dwrite(handles->speciesCountTimesDataset, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, times));

        // Cleanup some resources.
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
    }
}

void Hdf5File::appendLatticeTimeSeries(uint64_t replicate, const lm::io::LatticeTimeSeries& latticeTimeSeries)
{
    int numberEntries = latticeTimeSeries.lattices_size();
    if (latticeTimeSeries.times().shape(0) != numberEntries) InvalidArgException("latticeTimeSeries.times.shape", "Number of rows in time array incocnsistent with lattices array.");

    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open or create the lattice group.
    hid_t latticeGroupHandle, latticeTimesDatasetHandle;
    if (H5Lexists(replicateHandles->group, "Lattice", H5P_DEFAULT) > 0)
    {
        HDF5_EXCEPTION_CALL(latticeGroupHandle,H5Gopen2(replicateHandles->group, "Lattice", H5P_DEFAULT));
        HDF5_EXCEPTION_CALL(latticeTimesDatasetHandle,H5Dopen2(replicateHandles->group, "LatticeTimes", H5P_DEFAULT));
    }
    else
    {
        HDF5_EXCEPTION_CALL(latticeGroupHandle,H5Gcreate2(replicateHandles->group, "Lattice", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

        // Create the lattice times dataset.
        const unsigned int RANK=1;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        maxDims[0] = H5S_UNLIMITED;
        chunkDims[0] = 100;
        hid_t dataspaceHandle, propsHandle;
        HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(propsHandle,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(propsHandle, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(latticeTimesDatasetHandle,H5Dcreate2(replicateHandles->group, "LatticeTimes", H5T_IEEE_F64LE, dataspaceHandle, H5P_DEFAULT, propsHandle, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(propsHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
    }

    // Track the first lattice index.
    uint firstLatticeIndex;

    // Update the times dataset.
    {
        // Get the current size of the dataset.
        unsigned int RANK=1;
        hsize_t dims[RANK];
        hid_t dataspace_id;
        int result;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(latticeTimesDatasetHandle));
        HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspace_id, dims, NULL));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));

        // Save the first lattice index.
        firstLatticeIndex = dims[0];

        // Extend the dataset by the number of rows in the data set.
        dims[0] += numberEntries;
        HDF5_EXCEPTION_CHECK(H5Dset_extent(latticeTimesDatasetHandle, dims));

        // Create the memory dataset.
        hid_t memspace_id;
        hsize_t memDims[RANK];
        memDims[0] = numberEntries;
        HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));

        // Write the new data.
        ndarray<double> *times = NDArraySerializer::deserializeAllocate<double>(latticeTimeSeries.times());
        HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(latticeTimesDatasetHandle));
        hsize_t start[RANK], count[RANK];
        start[0] = dims[0]-numberEntries;
        firstLatticeIndex = start[0];
        count[0] = memDims[0];
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
        HDF5_EXCEPTION_CHECK(H5Dwrite(latticeTimesDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, times->values));

        // Cleanup some resources.
        delete times;
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
    }

    // Append the time to the times data set.
    for (int i=0; i<numberEntries; i++)
    {
        // Create the lattice data set.
        {
            ndarray<uint8_t>* particles = NDArraySerializer::deserializeAllocate<uint8_t>(latticeTimeSeries.lattices(i).particles());
            const unsigned int RANK=4;
            hsize_t dims[RANK], chunk[RANK];
            dims[0] = particles->shape[0];
            dims[1] = particles->shape[1];
            dims[2] = particles->shape[2];
            dims[3] = particles->shape[3];
            chunk[0] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,particles->shape[0]);
            chunk[1] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,particles->shape[1]);
            chunk[2] = min((uint)TUNE_LATTICE_GZIP_CHUNK_SIZE,particles->shape[2]);
            chunk[3] = particles->shape[3];
            hid_t dataspaceHandle, dcplHandle, datasetHandle;
            HDF5_EXCEPTION_CALL(dataspaceHandle,H5Screate_simple(RANK, dims, NULL));
            HDF5_EXCEPTION_CALL(dcplHandle,H5Pcreate(H5P_DATASET_CREATE));
            HDF5_EXCEPTION_CHECK(H5Pset_deflate (dcplHandle, TUNE_LATTICE_GZIP_COMPRESSION_LEVEL));
            HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcplHandle, RANK, chunk));
            char latticeDatasetName[11];
            snprintf(latticeDatasetName, sizeof(latticeDatasetName), "%010d", firstLatticeIndex+i);
            HDF5_EXCEPTION_CALL(datasetHandle,H5Dcreate2(latticeGroupHandle, latticeDatasetName, H5T_STD_U8LE, dataspaceHandle, H5P_DEFAULT, dcplHandle, H5P_DEFAULT));
            HDF5_EXCEPTION_CHECK(H5Dwrite(datasetHandle, H5T_NATIVE_UINT8, H5S_ALL, H5S_ALL, H5P_DEFAULT, particles->values));
            HDF5_EXCEPTION_CHECK(H5Dclose(datasetHandle));
            HDF5_EXCEPTION_CHECK(H5Pclose(dcplHandle));
            HDF5_EXCEPTION_CHECK(H5Sclose(dataspaceHandle));
            delete particles;
        }
    }

    // Cleanup some resources.
    HDF5_EXCEPTION_CHECK(H5Dclose(latticeTimesDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(latticeGroupHandle));

}

void Hdf5File::appendParameterValues(uint64_t replicate, lm::io::ParameterValues * parameterValues)
{
    if (parameterValues->value_size() != parameterValues->time_size()) throw InvalidArgException("parameterValues", "inconsistent number of entries");

    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open or create the parameter values group.
    hid_t pvGroupHandle;
    if ((pvGroupHandle=H5Gopen2(replicateHandles->group, "ParameterValues", H5P_DEFAULT)) < 0)
    {
        HDF5_EXCEPTION_CALL(pvGroupHandle,H5Gcreate2(replicateHandles->group, "ParameterValues", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }

    // Open or create the parameter values dataset.
    hid_t pvDatasetHandle;
    if ((pvDatasetHandle=H5Dopen2(pvGroupHandle, parameterValues->parameter().c_str(), H5P_DEFAULT)) < 0)
    {
        // Create the datasets.
        uint RANK=2;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        dims[1] = 2;
        maxDims[0] = H5S_UNLIMITED;
        maxDims[1] = 2;
        chunkDims[0] = 100;
        chunkDims[1] = 2;
        hid_t dataspace_id, dcpl_id;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(dcpl_id,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcpl_id, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(pvDatasetHandle,H5Dcreate2(pvGroupHandle, parameterValues->parameter().c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcpl_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    }

    // Get the current size of the dataset.
    unsigned int RANK=2;
    hsize_t dims[RANK];
    hid_t dataspace_id;
    int result;
    HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(pvDatasetHandle));
    HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspace_id, dims, NULL));
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));

    // Extend the dataset by the number of rows in the data set.
    dims[0] += parameterValues->value_size();
    HDF5_EXCEPTION_CHECK(H5Dset_extent(pvDatasetHandle, dims));

    // Create the memory dataset.
    hid_t memspace_id;
    hsize_t memDims[RANK];
    memDims[0] = parameterValues->value_size();
    memDims[1] = 1;
    HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));

    // Write the new times.
    HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(pvDatasetHandle));
    hsize_t start[RANK], count[RANK];
    start[0] = dims[0]-parameterValues->value_size();
    count[0] = parameterValues->value_size();
    start[1] = 0;
    count[1] = 1;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(pvDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, parameterValues->time().data()));
    start[1] = 1;
    count[1] = 1;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(pvDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, parameterValues->value().data()));

    // Cleanup some resources.
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
    HDF5_EXCEPTION_CHECK(H5Dclose(pvDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(pvGroupHandle));
}

void Hdf5File::setFirstPassageTimes(uint64_t replicate, const lm::io::FirstPassageTimes& fpt)
{
    // Get the data.
    uint32_t species = fpt.species();
    ndarray<int32_t>* newCounts = NDArraySerializer::deserializeAllocate<int32_t>(fpt.counts());
    ndarray<double>* newTimes = NDArraySerializer::deserializeAllocate<double>(fpt.first_passage_times());

    // Make sure the data is consistent.
    if (newCounts->shape.len != 1 || newCounts->shape != newTimes->shape) throw InvalidArgException("firstPassageTimes", "inconsistent number of first passage time entries");
    if (newCounts->shape == 0) throw InvalidArgException("firstPassageTimes", "no entries to save");

    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open the first passage time group.
    hid_t fptGroupHandle;
    if ((fptGroupHandle=H5Gopen2(replicateHandles->group, "FirstPassageTimes", H5P_DEFAULT)) < 0)
    {
        HDF5_EXCEPTION_CALL(fptGroupHandle,H5Gcreate2(replicateHandles->group, "FirstPassageTimes", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }

    // Construct a string representation of the species name.
    std::stringstream ss;
    ss.fill('0');
    ss.width(2);
    ss << species;
    string speciesString = ss.str();

    // Open the species group.
    hid_t speciesGroupHandle, countsDatasetHandle, timesDatasetHandle;
    if ((speciesGroupHandle=H5Gopen2(fptGroupHandle, speciesString.c_str(), H5P_DEFAULT)) >= 0)
    {
        // Open the data sets.
        HDF5_EXCEPTION_CALL(countsDatasetHandle,H5Dopen2(speciesGroupHandle, "Counts", H5P_DEFAULT));
        HDF5_EXCEPTION_CALL(timesDatasetHandle,H5Dopen2(speciesGroupHandle, "Times", H5P_DEFAULT));
    }
    else
    {
        // Create the group.
        HDF5_EXCEPTION_CALL(speciesGroupHandle,H5Gcreate2(fptGroupHandle, speciesString.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

        // Create the datasets.
        const uint RANK=1;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        maxDims[0] = H5S_UNLIMITED;
        chunkDims[0] = 100;
        hid_t dataspace_id, dcpl_id;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(dcpl_id,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcpl_id, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(countsDatasetHandle,H5Dcreate2(speciesGroupHandle, "Counts", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CALL(timesDatasetHandle,H5Dcreate2(speciesGroupHandle, "Times", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcpl_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    }

    // Read the lowest and highest first passage times in the file.
    const uint RANK=1;
    int result;
    hid_t countsDataspace;
    hsize_t countsDims[RANK];
    HDF5_EXCEPTION_CALL(countsDataspace,H5Dget_space(countsDatasetHandle));
    HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(countsDataspace, countsDims, NULL));
    int minCount=0, maxCount=0;
    if (countsDims[0] == 1)
    {
        HDF5_EXCEPTION_CHECK(H5Dread(countsDatasetHandle, H5T_NATIVE_INT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, &minCount));
        maxCount = minCount;
    }
    else if (countsDims[0] > 1)
    {
        // Create the memory dataspace.
        hid_t memDataspaceHandle;
        hsize_t memDims[RANK];
        memDims[0] = 1;
        HDF5_EXCEPTION_CALL(memDataspaceHandle,H5Screate_simple(RANK, memDims, NULL));

        hsize_t start[RANK], count[RANK];
        count[0] = 1;
        start[0] = 0;
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(countsDataspace, H5S_SELECT_SET, start, NULL, count, NULL));
        HDF5_EXCEPTION_CHECK(H5Dread(countsDatasetHandle, H5T_NATIVE_INT32, memDataspaceHandle, countsDataspace, H5P_DEFAULT, &minCount));
        start[0] = countsDims[0]-1;
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(countsDataspace, H5S_SELECT_SET, start, NULL, count, NULL));
        HDF5_EXCEPTION_CHECK(H5Dread(countsDatasetHandle, H5T_NATIVE_INT32, memDataspaceHandle, countsDataspace, H5P_DEFAULT, &maxCount));
        HDF5_EXCEPTION_CHECK(H5Sclose(memDataspaceHandle));
    }
    HDF5_EXCEPTION_CHECK(H5Sclose(countsDataspace));

    // Figure out the min and the max from the new counts.
    int32_t minNewCount=newCounts->get(0);
    int32_t maxNewCount=newCounts->get(newCounts->shape[0]-1);

    // If there are no existing records, just insert one large block.
    if (countsDims[0] == 0)
    {
        // Allocate the block.
        countsDims[0] = maxNewCount-minNewCount+1;
        uint32_t * counts = new uint32_t[countsDims[0]];
        double * times = new double[countsDims[0]];
        for (uint i=0; i<countsDims[0]; i++)
        {
            counts[i] = 0;
            times[i] = -1.0;
        }

        // Fill in the block.
        for (uint i=0; i<newCounts->shape[0]; i++)
        {
            uint32_t count = newCounts->get(i);
            counts[count-minNewCount] = count;
            times[count-minNewCount] = newTimes->get(i);
        }

        // Extend the datasets and write the block.
        HDF5_EXCEPTION_CHECK(H5Dset_extent(countsDatasetHandle, countsDims));
        HDF5_EXCEPTION_CHECK(H5Dwrite(countsDatasetHandle, H5T_NATIVE_UINT32, H5S_ALL, H5S_ALL, H5P_DEFAULT, counts));
        HDF5_EXCEPTION_CHECK(H5Dset_extent(timesDatasetHandle, countsDims));
        HDF5_EXCEPTION_CHECK(H5Dwrite(timesDatasetHandle, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, times));

        // Free the buffers.
        delete [] counts;
        delete [] times;
    }

    // Otherwise if there are only new records to append to the end, just add them
    else if (maxNewCount > maxCount && minNewCount > maxCount)
    {
        // Allocate the new buffer.
        hsize_t endingStart[RANK], endingCount[RANK];
        endingCount[0] = maxNewCount-maxCount;
        uint32_t * newEndingCounts = new uint32_t[endingCount[0]];
        double * newEndingTimes = new double[endingCount[0]];
        for (uint i=0; i<endingCount[0]; i++)
        {
            newEndingCounts[i] = 0;
            newEndingTimes[i] = -1.0;
        }

        // Fill in the buffer.
        for (int i=0; i<newCounts->shape[0]; i++)
        {
            int32_t count = newCounts->get(i);
            newEndingCounts[count-maxCount-1] = count;
            newEndingTimes[count-maxCount-1] = newTimes->get(i);
        }

        // Create the memory dataspace.
        hid_t memDataspaceHandle;
        hsize_t memDims[RANK];
        memDims[0] = endingCount[0];
        HDF5_EXCEPTION_CALL(memDataspaceHandle,H5Screate_simple(RANK, memDims, NULL));

        // Append the records.
        hid_t countsDataspaceHandle, timesDataspaceHandle;
        endingStart[0] = countsDims[0];
        countsDims[0] += endingCount[0];
        HDF5_EXCEPTION_CHECK(H5Dset_extent(countsDatasetHandle, countsDims));
        HDF5_EXCEPTION_CALL(countsDataspaceHandle,H5Dget_space(countsDatasetHandle));
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(countsDataspaceHandle, H5S_SELECT_SET, endingStart, NULL, endingCount, NULL));
        HDF5_EXCEPTION_CHECK(H5Dwrite(countsDatasetHandle, H5T_NATIVE_UINT32, memDataspaceHandle, countsDataspaceHandle, H5P_DEFAULT, newEndingCounts));
        HDF5_EXCEPTION_CHECK(H5Sclose(countsDataspaceHandle));
        HDF5_EXCEPTION_CHECK(H5Dset_extent(timesDatasetHandle, countsDims));
        HDF5_EXCEPTION_CALL(timesDataspaceHandle,H5Dget_space(timesDatasetHandle));
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(timesDataspaceHandle, H5S_SELECT_SET, endingStart, NULL, endingCount, NULL));
        HDF5_EXCEPTION_CHECK(H5Dwrite(timesDatasetHandle, H5T_NATIVE_DOUBLE, memDataspaceHandle, timesDataspaceHandle, H5P_DEFAULT, newEndingTimes));
        HDF5_EXCEPTION_CHECK(H5Sclose(timesDataspaceHandle));
        HDF5_EXCEPTION_CHECK(H5Sclose(memDataspaceHandle));

        // Free the buffers.
        if (newEndingCounts != NULL) delete [] newEndingCounts;
        if (newEndingTimes != NULL) delete [] newEndingTimes;
    }

    // Otherwise, we need to insert some records before and/or after the existing block so we have to copy the old data and rebuild.
    else
    {
        // Declare the new buffer.
        uint32_t combinedMinCount = min(minCount,minNewCount);
        uint32_t combinedMaxCount = max(maxCount,maxNewCount);
        size_t newSize = combinedMaxCount-combinedMinCount+1;
        uint32_t * counts = new uint32_t[newSize];
        double * times = new double[newSize];
        for (uint i=0; i<newSize; i++)
        {
            counts[i] = 0;
            times[i] = -1.0;
        }

        // Copy the old data into the buffer.
        hid_t memDataspaceHandle;
        hsize_t memDims[RANK];
        memDims[0] = maxCount-minCount+1;
        HDF5_EXCEPTION_CALL(memDataspaceHandle,H5Screate_simple(RANK, memDims, NULL));
        HDF5_EXCEPTION_CHECK(H5Dread(countsDatasetHandle, H5T_NATIVE_INT32, memDataspaceHandle, H5S_ALL, H5P_DEFAULT, &counts[minCount-combinedMinCount]));
        HDF5_EXCEPTION_CHECK(H5Dread(timesDatasetHandle, H5T_NATIVE_DOUBLE, memDataspaceHandle, H5S_ALL, H5P_DEFAULT, &times[minCount-combinedMinCount]));
        HDF5_EXCEPTION_CHECK(H5Sclose(memDataspaceHandle));

        // Fill in the buffer with the new data.
        for (uint i=0; i<newCounts->shape[0]; i++)
        {
            int32_t count = newCounts->get(i);
            if (count < minCount || count > maxCount)
            {
                newCounts[count-combinedMinCount] = count;
                newTimes[count-combinedMinCount] = newTimes->get(i);
            }
            else
            {
                throw InvalidArgException("firstPassageTimes", "contained duplicates of existing counts", replicate, minCount, minNewCount, maxCount, maxNewCount);
            }
        }

        // Update the dataset with the combined data.
        memDims[0] = newSize;
        countsDims[0] = newSize;
        HDF5_EXCEPTION_CALL(memDataspaceHandle,H5Screate_simple(RANK, memDims, NULL));
        HDF5_EXCEPTION_CHECK(H5Dset_extent(countsDatasetHandle, countsDims));
        HDF5_EXCEPTION_CHECK(H5Dwrite(countsDatasetHandle, H5T_NATIVE_UINT32, memDataspaceHandle, H5S_ALL, H5P_DEFAULT, newCounts));
        HDF5_EXCEPTION_CHECK(H5Dset_extent(timesDatasetHandle, countsDims));
        HDF5_EXCEPTION_CHECK(H5Dwrite(timesDatasetHandle, H5T_NATIVE_DOUBLE, memDataspaceHandle, H5S_ALL, H5P_DEFAULT, newTimes));
        HDF5_EXCEPTION_CHECK(H5Sclose(memDataspaceHandle));

        // Free the buffers.
        if (counts != NULL) delete [] newCounts;
        if (times != NULL) delete [] newTimes;
    }

    // Close any resources.
    HDF5_EXCEPTION_CHECK(H5Dclose(countsDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Dclose(timesDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(speciesGroupHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(fptGroupHandle));

    // Free the ndarrays.
    if (newCounts != NULL) delete newCounts;
    if (newTimes != NULL) delete newTimes;
}

/*void SimulationFile::appendSpatialModelObjects(unsigned int replicate, lm::input::SpatialModel * model) throw(HDF5Exception,InvalidArgException)
{
    if (model->sphere_xc_size() != model->sphere_type_size()) throw InvalidArgException("model", "inconsistent number of sphere entries");
    if (model->sphere_yc_size() != model->sphere_type_size()) throw InvalidArgException("model", "inconsistent number of sphere entries");
    if (model->sphere_zc_size() != model->sphere_type_size()) throw InvalidArgException("model", "inconsistent number of sphere entries");
    if (model->sphere_radius_size() != model->sphere_type_size()) throw InvalidArgException("model", "inconsistent number of sphere entries");

    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open or create the group.
    hid_t modelHandle, spatialHandle;
    if ((modelHandle=H5Gopen2(replicateHandles->group, "Model", H5P_DEFAULT)) < 0)
    {
        HDF5_EXCEPTION_CALL(modelHandle,H5Gcreate2(replicateHandles->group, "Model", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }
    if ((spatialHandle=H5Gopen2(modelHandle, "Spatial", H5P_DEFAULT)) < 0)
    {
        HDF5_EXCEPTION_CALL(spatialHandle,H5Gcreate2(modelHandle, "Spatial", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
    }

    // Open or create the datasets.
    hid_t spheresDatasetHandle;
    if ((spheresDatasetHandle=H5Dopen2(spatialHandle, "Spheres", H5P_DEFAULT)) < 0)
    {
        // Create the datasets.
        uint RANK=2;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        dims[1] = 5;
        maxDims[0] = H5S_UNLIMITED;
        maxDims[1] = 5;
        chunkDims[0] = 100;
        chunkDims[1] = 5;
        hid_t dataspace_id, dcpl_id;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(dcpl_id,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcpl_id, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(spheresDatasetHandle,H5Dcreate2(spatialHandle, "Spheres", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcpl_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    }

    // Get the current size of the dataset.
    unsigned int RANK=2;
    hsize_t dims[RANK];
    hid_t dataspace_id;
    int result;
    HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(spheresDatasetHandle));
    HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspace_id, dims, NULL));
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));

    // Extend the dataset by the number of rows in the data set.
    dims[0] += model->sphere_type_size();
    HDF5_EXCEPTION_CHECK(H5Dset_extent(spheresDatasetHandle, dims));

    // Create the memory dataset.
    hid_t memspace_id;
    hsize_t memDims[RANK];
    memDims[0] = model->sphere_type_size();
    memDims[1] = 1;
    HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));

    // Write the values.
    HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(spheresDatasetHandle));
    hsize_t start[RANK], count[RANK];
    start[0] = dims[0]-model->sphere_type_size();
    count[0] = model->sphere_type_size();
    start[1] = 0;
    count[1] = 1;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, model->sphere_xc().data()));
    start[1] = 1;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, model->sphere_yc().data()));
    start[1] = 2;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, model->sphere_zc().data()));
    start[1] = 3;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, model->sphere_radius().data()));
    start[1] = 4;
    HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
    HDF5_EXCEPTION_CHECK(H5Dwrite(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, model->sphere_type().data()));

    // Cleanup some resources.
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
    HDF5_EXCEPTION_CHECK(H5Dclose(spheresDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(modelHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(spatialHandle));
}

void SimulationFile::getSpatialModelObjects(unsigned int replicate, lm::input::SpatialModel * model) throw(HDF5Exception)
{
    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Open the group.
    hid_t modelHandle, spatialHandle;
    HDF5_EXCEPTION_CALL(modelHandle,H5Gopen2(replicateHandles->group, "Model", H5P_DEFAULT));
    HDF5_EXCEPTION_CALL(spatialHandle,H5Gopen2(modelHandle, "Spatial", H5P_DEFAULT));

    // Open the datasets.
    hid_t spheresDatasetHandle;
    HDF5_EXCEPTION_CALL(spheresDatasetHandle,H5Dopen2(spatialHandle, "Spheres", H5P_DEFAULT));

    // Get the current size of the dataset.
    const uint RANK=2;
    const uint NUM_COLS=5;
    hsize_t dims[RANK];
    hid_t dataspace_id;
    int result;
    HDF5_EXCEPTION_CALL(dataspace_id,H5Dget_space(spheresDatasetHandle));
    HDF5_EXCEPTION_CALL(result,H5Sget_simple_extent_dims(dataspace_id, dims, NULL));

    // Create the memory dataset.
    hid_t memspace_id;
    hsize_t memDims[RANK];
    memDims[0] = TUNE_SPATIAL_MODEL_OJBECT_READ_BUFFER_SIZE;
    memDims[1] = NUM_COLS;
    HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));
    double * memData = new double[TUNE_SPATIAL_MODEL_OJBECT_READ_BUFFER_SIZE*NUM_COLS];

    // Loop over the data and read it in one section at a time.
    hsize_t start[RANK], count[RANK];
    start[0] = 0;
    start[1] = 0;
    count[0] = TUNE_SPATIAL_MODEL_OJBECT_READ_BUFFER_SIZE;
    count[1] = NUM_COLS;
    while (start[0] < dims[0])
    {
        // Figure out how many rows to read.
        if (start[0]+count[0] > dims[0])
        {
            // This is the last read, so modify its size.
            count[0] = dims[0]-start[0];

            // Modify the memory dataspace as well.
            HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
            memDims[0] = count[0];
            HDF5_EXCEPTION_CALL(memspace_id,H5Screate_simple(RANK, memDims, NULL));
        }

        // Read the rows.
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, start, NULL, count, NULL));
        HDF5_EXCEPTION_CHECK(H5Dread(spheresDatasetHandle, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, (void *)memData));

        // Add them to the model.
        for (uint index=0; index<count[0]*NUM_COLS; )
        {
            model->add_sphere_xc(memData[index++]);
            model->add_sphere_yc(memData[index++]);
            model->add_sphere_zc(memData[index++]);
            model->add_sphere_radius(memData[index++]);
            model->add_sphere_type(memData[index++]);
        }

        // Move to the next section.
        start[0] += count[0];
    }

    // Cleanup some resources.
    if (memData != NULL) delete []  memData; memData = NULL;
    HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    HDF5_EXCEPTION_CHECK(H5Sclose(memspace_id));
    HDF5_EXCEPTION_CHECK(H5Dclose(spheresDatasetHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(modelHandle));
    HDF5_EXCEPTION_CHECK(H5Gclose(spatialHandle));
}*/


vector<double> Hdf5File::getLatticeTimes(uint64_t replicate)
{
    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Read the initial species counts.
    int RANK=1, rank;
    hsize_t dims[RANK];
    H5T_class_t type;
    size_t size;
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(replicateHandles->group, "LatticeTimes", &rank));
    if (rank != RANK) throw Exception("Invalid dataset rank", filename.c_str(), "LatticeTimes");
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(replicateHandles->group, "LatticeTimes", dims, &type, &size));
    if (size != sizeof(double)) throw Exception("Invalid dataset type", filename.c_str(), "LatticeTimes");
    double * timesBuffer = new double[dims[0]];
    HDF5_EXCEPTION_CHECK(H5LTread_dataset_double(replicateHandles->group, "LatticeTimes", timesBuffer));
    vector<double> times;
    for (uint i=0; i<dims[0]; i++)
        times.push_back(timesBuffer[i]);
    delete [] timesBuffer;
    return times;
}

void Hdf5File::getLattice(uint64_t replicate, unsigned int latticeIndex, lm::rdme::Lattice * lattice)
{
    ReplicateHandles * replicateHandles = openReplicateHandles(replicate);

    // Read the lattice data.
    char latticeDatasetName[19];
    snprintf(latticeDatasetName, sizeof(latticeDatasetName), "Lattice/%010d", latticeIndex);
    int RANK=4, rank;
    hsize_t dims[RANK];
    H5T_class_t type;
    size_t size;
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_ndims(replicateHandles->group, latticeDatasetName, &rank));
    if (rank != RANK) throw Exception("Invalid dataset rank", filename.c_str(), latticeDatasetName);
    HDF5_EXCEPTION_CHECK(H5LTget_dataset_info(replicateHandles->group, latticeDatasetName, dims, &type, &size));
    if (size != sizeof(uint8_t)) throw Exception("Invalid dataset type", filename.c_str(), latticeDatasetName);
    if (lattice->getSize().x != dims[0] || lattice->getSize().y != dims[1] || lattice->getSize().z != dims[2] || lattice->getMaxOccupancy() != dims[3]) throw Exception("Invalid lattice dimensions", filename.c_str());
    size_t particlesBufferSize = dims[0]*dims[1]*dims[2]*dims[3];
    uint8_t * particlesBuffer = new uint8_t[particlesBufferSize];
    HDF5_EXCEPTION_CHECK(H5LTread_dataset(replicateHandles->group, latticeDatasetName, H5T_STD_U8LE, particlesBuffer));

    // Set the lattice object using the data.
    //TODO lattice->setFromRowMajorByteData(particlesBuffer, particlesBufferSize);

    // Free the intermediate lattice data.
    delete [] particlesBuffer;
}

void Hdf5File::closeReplicate(uint64_t replicate)
{
    ReplicateHandleMap::Key replicateKey(recordNamePrefix, replicate);
    ReplicateHandleMap::iterator it = openReplicates.find(replicateKey);
    if (it != openReplicates.end())
    {
        ReplicateHandles * handles = it->second;
        closeReplicateHandles(handles);
        delete handles;
        openReplicates.erase(it);
    }
}

void Hdf5File::closeAllReplicates()
{
    for (ReplicateHandleMap::iterator it=openReplicates.begin(); it != openReplicates.end(); it++)
    {
        ReplicateHandles * handles = it->second;
        closeReplicateHandles(handles);
        delete handles;
    }
    openReplicates.clear();
}


Hdf5File::ReplicateHandles* Hdf5File::openReplicateHandles(uint64_t replicate)
{
    // See if the replicate is already open.
    ReplicateHandleMap::Key replicateKey(recordNamePrefix, replicate);
    ReplicateHandleMap::iterator it = openReplicates.find(replicateKey);
    if (it == openReplicates.end())
    {
        ReplicateHandles * handles;

        // Construct a string representation of the replicate name.
        std::stringstream ss;
        ss.fill('0');
        ss.width(7);
        ss << replicate;
        string replicateString = ss.str();

        // Turn off exception handling again to work around odd behavior in hdf5 1.8.4.
        HDF5_EXCEPTION_CHECK(H5Eset_auto2(H5E_DEFAULT, NULL, NULL));

        // Open the replicate group.
        hid_t groupHandle;
        if ((groupHandle=H5Gopen2(simulationsGroup, replicateString.c_str(), H5P_DEFAULT)) >= 0)
        {
            // Open the rest of the handles.
            handles = new ReplicateHandles();
            handles->group = groupHandle;
            HDF5_EXCEPTION_CALL(handles->speciesCountsDataset,H5Dopen2(handles->group, "SpeciesCounts", H5P_DEFAULT));
            HDF5_EXCEPTION_CALL(handles->speciesCountTimesDataset,H5Dopen2(handles->group, "SpeciesCountTimes", H5P_DEFAULT));
        }
        else
        {
            // If we couldn't open it, try to create a new one.
            handles = createReplicateHandles(replicateString);
        }

        // Add it to the map.
        openReplicates[replicateKey] = handles;
        return handles;
    }
    else
    {
        // Use the already open replicate.
        return it->second;
    }

}

Hdf5File::ReplicateHandles* Hdf5File::createReplicateHandles(string replicateString)
{
    // Make sure the model is loaded, since need the number of species.
    loadModel();

    ReplicateHandles * handles = new ReplicateHandles();

    // Create the group.
    HDF5_EXCEPTION_CALL(handles->group,H5Gcreate2(simulationsGroup, replicateString.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));

    // Create the species counts dataset
    {
        const unsigned int RANK=2;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        dims[1] = numberSpecies;
        maxDims[0] = H5S_UNLIMITED;
        maxDims[1] = numberSpecies;
        chunkDims[0] = 100;
        chunkDims[1] = numberSpecies;
        hid_t dataspace_id, dcpl_id;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(dcpl_id,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcpl_id, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(handles->speciesCountsDataset,H5Dcreate2(handles->group, "SpeciesCounts", H5T_STD_I32LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcpl_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    }

    // Create the species count times dataset
    {
        const unsigned int RANK=1;
        hsize_t dims[RANK], maxDims[RANK], chunkDims[RANK];
        dims[0] = 0;
        maxDims[0] = H5S_UNLIMITED;
        chunkDims[0] = 100;
        hid_t dataspace_id, dcpl_id;
        HDF5_EXCEPTION_CALL(dataspace_id,H5Screate_simple(RANK, dims, maxDims));
        HDF5_EXCEPTION_CALL(dcpl_id,H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(dcpl_id, RANK, chunkDims));
        HDF5_EXCEPTION_CALL(handles->speciesCountTimesDataset,H5Dcreate2(handles->group, "SpeciesCountTimes", H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, dcpl_id, H5P_DEFAULT));
        HDF5_EXCEPTION_CHECK(H5Pclose(dcpl_id));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace_id));
    }

    return handles;
}

void Hdf5File::closeReplicateHandles(ReplicateHandles * handles)
{
    HDF5_EXCEPTION_CHECK(H5Gclose(handles->group));
    HDF5_EXCEPTION_CHECK(H5Dclose(handles->speciesCountsDataset));
    HDF5_EXCEPTION_CHECK(H5Dclose(handles->speciesCountTimesDataset));
    handles->group = H5I_INVALID_HID;
    handles->speciesCountsDataset = H5I_INVALID_HID;
    handles->speciesCountTimesDataset = H5I_INVALID_HID;
}

void Hdf5File::setRecordNamePrefix(const string& newRecordNamePrefix)
{
    if (newRecordNamePrefix!=recordNamePrefix)
    {
        recordNamePrefix.assign(newRecordNamePrefix);

        if (recordNamePrefix.empty())
        {
            simulationsGroup = initGroup("Simulations");
        }
        else
        {
            simulationsGroup = initGroup(pathJoin(recordNamePrefix, "Simulations"));
        }
    }
}

hsize_t Hdf5File::setDatasetFromNDArray(const std::string& groupPath, const std::string& datasetName, const robertslab::pbuf::NDArray& ndarrayRef, hid_t rootGroup)
{
    // a descriptor that we'll pass to the lower level output function
    DatasetDescriptor datasetDescriptor(groupPath, datasetName, ndarrayRef, rootGroup);

    // now that we have the data and the shape, call the generalized dataset writing function
    hsize_t rows = setDataset(datasetDescriptor);

    // return the number of rows written out
    return rows;
}

void Hdf5File::setDatasetFromNDArrayReplicate(uint64_t replicate, const std::string& groupRelativePath, const std::string& datasetName, const robertslab::pbuf::NDArray& ndarray, bool condensed)
{
    if (condensed)
    {
        setDatasetFromNDArrayReplicateCondensed(replicate, groupRelativePath, datasetName, ndarray);
    }
    else
    {
        ReplicateHandles * replicateHandles = openReplicateHandles(replicate);
        setDatasetFromNDArray(groupRelativePath, datasetName, ndarray, replicateHandles->group);
    }
}

// condensed version of the generalized NDArray hdf5 output. Condensed in the sense that it shoves all of the data into as few separate groups and datasets as possible
void Hdf5File::setDatasetFromNDArrayReplicateCondensed(uint64_t replicate, const std::string& groupRelativePath, const std::string& datasetName, const robertslab::pbuf::NDArray& ndarray)
{
    // write out the dataset directly to prefix/Simulations/groupRelativePath and get the number of rows written
    hsize_t rows = setDatasetFromNDArray(groupRelativePath, datasetName, ndarray, simulationsGroup);

    // create a 1D array containing one repeat of the trajectoryID for each row in ndarray
    std::vector<uint64_t> trajectoryIDs(rows, replicate);

    // write out the trajectoryID dataset we just created
    std::string trajectoryIDDatasetName = datasetName + "_-_TrajectoryIDs";
    setDatasetFromContainer(groupRelativePath, trajectoryIDDatasetName, trajectoryIDs, simulationsGroup);
}

hsize_t Hdf5File::setDataset(const DatasetDescriptor& dd)
{
    // initialize the group we'll be storing the NDArray dataset in
    hid_t group = initGroup(dd.groupPath, dd.rootGroup);

    // declare the HDF5 boilerplate variables
    uint RANK(dd.shape.len);
    hid_t dataspace, dataset, filespace, memspace, prop;
    hsize_t chunkdims[RANK], dims[RANK], dimsr[RANK], dimstotal[RANK], maxdims[RANK], offset[RANK];

    // We'd like to have the option to delete NDArray's dataset if it already exists, but HDF5 apparently can't really delete anything so leave it commented for now.
    /* if (H5Lexists(group, groupName.c_str(), H5P_DEFAULT))
    {
        HDF5_EXCEPTION_CHECK(H5Ldelete(group, groupName.c_str(), H5P_DEFAULT));
    } */

    // write or extend the NDArray dataset
    dims[0] = RANK > 0 ? dd.shape[0] : 0;
    chunkdims[0] = 10;
    maxdims[0] = H5S_UNLIMITED;
    for (int i=1; i<RANK; i++)
    {
        dims[i] = dd.shape[i];
        chunkdims[i] = dims[i];
        maxdims[i] = dims[i];
    }

    // if the dataset exists, extend it
    if ((dataset = H5Dopen2(group, dd.datasetName.c_str(), H5P_DEFAULT))>=0)
    {
        HDF5_EXCEPTION_CALL(prop, H5Dget_create_plist(dataset));

        HDF5_EXCEPTION_CALL(filespace, H5Dget_space(dataset));
        HDF5_EXCEPTION_CHECK(H5Sget_simple_extent_dims(filespace, dimsr, NULL));
        /* Extend the dataset */
        dimstotal[0] = dimsr[0] + dims[0];
        if (RANK==2) {dimstotal[1] = dimsr[1];}
        HDF5_EXCEPTION_CHECK(H5Dset_extent(dataset, dimstotal));
        // reopen the now-extended dataset's filespace
        HDF5_EXCEPTION_CALL(filespace, H5Dget_space(dataset));
        /* Select a hyperslab in extended portion of dataset  */
        offset[0] = dimsr[0];
        if (RANK==2) {offset[1] = 0;}
        HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(filespace, H5S_SELECT_SET, offset, NULL, dims, NULL));
        /* Define memory space */
        HDF5_EXCEPTION_CALL(memspace, H5Screate_simple(RANK, dims, NULL));
        HDF5_EXCEPTION_CHECK(H5Dwrite(dataset, dd.hdf5Type, memspace, filespace, H5P_DEFAULT, dd.data));

        HDF5_EXCEPTION_CHECK(H5Dclose(dataset));
        HDF5_EXCEPTION_CHECK(H5Sclose(memspace));
        HDF5_EXCEPTION_CHECK(H5Sclose(filespace));
    }
        // otherwise, create the dataset
    else
    {
        /* Create the dataField space with unlimited dimensions. */
        HDF5_EXCEPTION_CALL(dataspace, H5Screate_simple(RANK, dims, maxdims));
        /* Modify dataset creation properties, i.e. enable chunking  */
        HDF5_EXCEPTION_CALL(prop, H5Pcreate(H5P_DATASET_CREATE));
        HDF5_EXCEPTION_CHECK(H5Pset_chunk(prop, RANK, chunkdims));
        /* Create a new dataset within the file using chunk creation properties.  */
        dataset = H5Dcreate2(group, dd.datasetName.c_str(), dd.hdf5Type, dataspace, H5P_DEFAULT, prop, H5P_DEFAULT);
        /* Write dataField to dataset */
        HDF5_EXCEPTION_CHECK(H5Dwrite(dataset, dd.hdf5Type, H5S_ALL, H5S_ALL, H5P_DEFAULT, dd.data));

        HDF5_EXCEPTION_CHECK(H5Dclose(dataset));
        HDF5_EXCEPTION_CHECK(H5Pclose(prop));
        HDF5_EXCEPTION_CHECK(H5Sclose(dataspace));
    }

    return dims[0];
}

}
}
}

/* Generic code to read a dataset with a hyperslab.
 * int ndims;
hsize_t dims[4];
hsize_t memdims[1];
hsize_t start[4], count[4];
hid_t dataset, dataspace, memspace;
HDF5_EXCEPTION_CALL(dataset,H5Dopen2(file, "/Model/Diffusion/Lattice", H5P_DEFAULT));
HDF5_EXCEPTION_CALL(dataspace,H5Dget_space(dataset));
ndims=H5Sget_simple_extent_ndims(dataspace);
if (ndims != 4) throw Exception("Invalid dataset dimensions", filename.c_str(), "/Model/Diffusion/Lattice");
HDF5_EXCEPTION_CALL(ndims,H5Sget_simple_extent_dims(dataspace, dims, NULL));
if (dims[0]*dims[1]*dims[2]*dims[3] != latticeSize) throw InvalidArgException("lattice", "incorrect lattice size");
start[0] = 0;
start[1] = 0;
start[2] = 0;
start[3] = 0;
count[0] = dims[0];
count[1] = dims[1];
count[2] = dims[2];
count[3] = dims[3];
HDF5_EXCEPTION_CHECK(H5Sselect_hyperslab(dataspace, H5S_SELECT_SET, start, NULL, count, NULL));
memdims[0] = dims[0]*dims[1]*dims[2]*dims[3];
HDF5_EXCEPTION_CALL(memspace,H5Screate_simple(1, memdims, NULL));
HDF5_EXCEPTION_CHECK(H5Dread(dataset, H5T_NATIVE_UINT8, memspace, dataspace, H5P_DEFAULT, (void *)lattice));
HDF5_EXCEPTION_CHECK(H5Sclose(memspace));
HDF5_EXCEPTION_CHECK(H5Sclose(dataspace));
HDF5_EXCEPTION_CHECK(H5Dclose(dataset));
*/
