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

#include <string>

#include "lm/String.h"
#include "lm/ClassFactory.h"
#include "lm/input/Input.pb.h"
#include "lm/io/hdf5/Hdf5InputReader.h"
#include "lm/io/hdf5/SimulationFile.h"
#include "lm/types/BoundaryConditions.pb.h"


namespace lm {
namespace io {
namespace hdf5 {

bool Hdf5InputReader::registered=Hdf5InputReader::registerClass();

bool Hdf5InputReader::registerClass()
{
    lm::ClassFactory::getInstance().registerClass("lm::input::InputReader","lm::io::hdf5::Hdf5InputReader",&Hdf5InputReader::allocateObject);
    return true;
}

void* Hdf5InputReader::allocateObject()
{
    return new Hdf5InputReader();
}

Hdf5InputReader::Hdf5InputReader()
{
}

Hdf5InputReader::~Hdf5InputReader()
{
}

bool Hdf5InputReader::canReadFile(std::string filename)
{
    return (endsWith(filename, ".lm") || endsWith(filename, ".h5")) && lm::io::hdf5::Hdf5File::isValidFile(filename);
}

void Hdf5InputReader::readFileInto(std::string filename, lm::input::Input* input)
{
    lm::io::hdf5::Hdf5File hdf5File = lm::io::hdf5::Hdf5File(filename);
    readInput(hdf5File, input);

}

void Hdf5InputReader::readInput(const lm::io::hdf5::Hdf5File& file, lm::input::Input* input)
{
    // Read the options.
    readSimulationOptions(file, input);
    readOutputOptions(file, input);

    // Read the models.
    readReactionModel(file, input);
    readDiffusionModel(file, input);

    // Read the other data sets.
    readOrderParameters(file, input);
    readTilings(file, input);

    // Read the restart options.
    readRestart(file, input);

    // Read the ffpilot options.
    readFFPilotOptions(file, input);
}

void Hdf5InputReader::readReactionModel(const lm::io::hdf5::Hdf5File& file, lm::input::Input* input)
{
    if (file.hasReactionModel())
    {
        file.getReactionModel(input->mutable_reaction_model());
    }
}

void Hdf5InputReader::readDiffusionModel(const lm::io::hdf5::Hdf5File& file, lm::input::Input* input)
{
    if (file.hasDiffusionModel())
    {
        // Read the model.
        file.getDiffusionModel(input->mutable_diffusion_model());

        // See if we need to fill in the boundary conditions from the parameters.
        if (file.hasParameter("boundaryConditions"))
        {
            // Parse the boundary conditions string.
            lm::types::BoundaryConditions* bc = input->mutable_diffusion_model()->mutable_boundary_conditions();
            string bcString = file.getParameter("boundaryConditions");
            if (!parseBoundaryConditionsString(bcString, bc))
                THROW_EXCEPTION(InputException, "Could not parse boundaryConditions parameter: %s", bcString.c_str());

            // Parse some other parameters.
            if (file.hasParameter("boundarySite")) bc->set_boundary_site(atoi(file.getParameter("boundarySite").c_str()));
            if (file.hasParameter("boundarySpecies")) bc->set_boundary_site(atoi(file.getParameter("boundarySpecies").c_str()));
            if (file.hasParameter("boundaryConcentration")) bc->set_boundary_concentration(atof(file.getParameter("boundaryConcentration").c_str()));

            // Parse the boundary gradient, if there is one.
            if (file.hasBoundaryGradient()) file.getBoundaryGradient(input->mutable_diffusion_model()->mutable_boundary_conditions());
        }
    }
}

bool Hdf5InputReader::parseBoundaryConditionsString(string arg, lm::types::BoundaryConditions* bc)
{
    lm::types::BoundaryConditions::BoundaryConditionsType type;

    // See if it is a global boundary condition.
    if (lm::types::BoundaryConditions_BoundaryConditionsType_Parse(arg, &type))
    {
        bc->set_global(type);
        return true;
    }

    // See if there are axis specific boundary conditions.
    char* argbuf = new char[arg.size()+1];
    memset(argbuf,0,arg.size()+1);
    strcpy(argbuf,arg.c_str());
    char* pch = strtok(argbuf,",");
    while (pch != NULL)
    {
        if (strlen(pch) >= 3 && (pch[0] == 'x' || pch[0] == 'y' || pch[0] == 'z') && pch[1] == ':')
        {
            // Parse the axis-specific type.
            if (!lm::types::BoundaryConditions_BoundaryConditionsType_Parse(std::string(pch+2), &type))
            {
                delete[] argbuf;
                return false;
            }

            // Set the axis value.
            pch[1] = '\0';
            std::string axis=pch;
            if (axis == "x")
            {
                bc->set_axis_specific_boundaries(true);
                bc->set_x_plus(type);
                bc->set_x_minus(type);
            }
            else if (axis == "y")
            {
                bc->set_axis_specific_boundaries(true);
                bc->set_y_plus(type);
                bc->set_y_minus(type);
            }
            else if (axis == "z")
            {
                bc->set_axis_specific_boundaries(true);
                bc->set_z_plus(type);
                bc->set_z_minus(type);
            }
            else
            {
                delete[] argbuf;
                return false;
            }
        }
        else if (strlen(pch) >= 4 && ((pch[0] == '+' || pch[0] == '-') && (pch[1] == 'x' || pch[1] == 'y' || pch[1] == 'z')) && pch[2] == ':')
        {
            // Parse the axis-specific type.
            if (!lm::types::BoundaryConditions_BoundaryConditionsType_Parse(std::string(pch+3), &type))
            {
                delete[] argbuf;
                return false;
            }

            // Set the axis value.
            pch[2] = '\0';
            std::string axis=pch;
            if (axis == "+x" && type != lm::types::BoundaryConditions::PERIODIC)
            {
                bc->set_axis_specific_boundaries(true);
                bc->set_x_plus(type);
            }
            else if (axis == "-x" && type != lm::types::BoundaryConditions::PERIODIC)
            {
                bc->set_axis_specific_boundaries(true);
                bc->set_x_minus(type);
            }
            else if (axis == "+y" && type != lm::types::BoundaryConditions::PERIODIC)
            {
                bc->set_axis_specific_boundaries(true);
                bc->set_y_plus(type);
            }
            else if (axis == "-y" && type != lm::types::BoundaryConditions::PERIODIC)
            {
                bc->set_axis_specific_boundaries(true);
                bc->set_y_minus(type);
            }
            else if (axis == "+z" && type != lm::types::BoundaryConditions::PERIODIC)
            {
                bc->set_axis_specific_boundaries(true);
                bc->set_z_plus(type);
            }
            else if (axis == "-z")
            {
                bc->set_axis_specific_boundaries(true);
                bc->set_z_minus(type);
            }
            else
            {
                delete[] argbuf;
                return false;
            }
        }
        else
        {
            delete[] argbuf;
            return false;
        }
        pch = strtok(NULL,",");
    }
    delete[] argbuf;
    return bc->axis_specific_boundaries();
}

void Hdf5InputReader::readOrderParameters(const lm::io::hdf5::Hdf5File& file, lm::input::Input* input)
{
    if (file.hasOrderParameters())
    {
        file.getOrderParameters(input->mutable_order_parameters());
    }
}

// Get the tilings.
void Hdf5InputReader::readTilings(const lm::io::hdf5::Hdf5File& file, lm::input::Input* input)
{
    if (file.hasTilings())
    {
        file.getTilings(input->mutable_tilings());
    }
}

void Hdf5InputReader::readFFPilotOptions(const lm::io::hdf5::Hdf5File &file, lm::input::Input *input)
{
    if (file.hasFFPilotOptions())
    {
        file.getFFPilotOptions(input->mutable_ffpilot_options());
    }
}

void Hdf5InputReader::readSimulationOptions(const lm::io::hdf5::Hdf5File& file, lm::input::Input* input)
{
    // See if max time limit is set.
    if (file.hasParameter("maxTime")) input->mutable_simulation_options()->set_time_limit(atof(file.getParameter("maxTime").c_str()));

    // See if any species upper limits are set.
    if (file.hasParameter("speciesUpperLimitList"))
    {
        vector<string> limitStrings = parseStringList(file.getParameter("speciesUpperLimitList"), ",");
        for (size_t i=0; i<limitStrings.size(); i++)
        {
            vector<int32_t> entries = parseIntList(limitStrings[i], ":");
            if (entries.size() == 2 && entries[0] >= 0)
            {
                lm::input::SpeciesLimit* limit = input->mutable_simulation_options()->add_species_upper_limit();
                limit->set_species_index(uint32_t(entries[0]));
                limit->set_limit_value(entries[1]);
            }
            else
            {
                THROW_EXCEPTION(lm::InvalidArgException, "speciesUpperLimitList was invalid: %s", file.getParameter("speciesUpperLimitList").c_str());
            }
        }
    }

    // See if any species lower limits are set.
    if (file.hasParameter("speciesLowerLimitList"))
    {
        vector<string> limitStrings = parseStringList(file.getParameter("speciesLowerLimitList"), ",");
        for (size_t i=0; i<limitStrings.size(); i++)
        {
            vector<int32_t> entries = parseIntList(limitStrings[i], ":");
            if (entries.size() == 2 && entries[0] >= 0)
            {
                lm::input::SpeciesLimit* limit = input->mutable_simulation_options()->add_species_lower_limit();
                limit->set_species_index(uint32_t(entries[0]));
                limit->set_limit_value(entries[1]);
            }
            else
            {
                THROW_EXCEPTION(lm::InvalidArgException, "speciesLowerLimitList was invalid: %s", file.getParameter("speciesLowerLimitList").c_str());
            }
        }
    }

    // See if any species reflecting barriers are set.
    if (file.hasParameter("speciesReflectingBarrierList"))
    {
        vector<string> limitStrings = parseStringList(file.getParameter("speciesReflectingBarrierList"), ",");
        for (size_t i=0; i<limitStrings.size(); i++)
        {
            vector<int32_t> entries = parseIntList(limitStrings[i], ":");
            if (entries.size() == 2 && entries[0] >= 0)
            {
                lm::input::SpeciesLimit* limit = input->mutable_simulation_options()->add_species_reflecting_barrier();
                limit->set_species_index(uint32_t(entries[0]));
                limit->set_limit_value(entries[1]);
            }
            else
            {
                THROW_EXCEPTION(lm::InvalidArgException, "speciesReflectingBarrierList was invalid: %s", file.getParameter("speciesReflectingBarrierList").c_str());
            }
        }
    }

    // See if any order parameter reflecting barriers are set.
    if (file.hasParameter("orderParameterReflectingBarrierList"))
    {
        vector<string> limitStrings = parseStringList(file.getParameter("orderParameterReflectingBarrierList"), ",");
        for (size_t i=0; i<limitStrings.size(); i++)
        {
            vector<string> entries = parseStringList(limitStrings[i], ":");
            if (entries.size() == 2)
            {
                lm::input::OrderParameterLimit* limit = input->mutable_simulation_options()->add_order_parameter_reflecting_barrier();
                limit->set_order_parameter_index(static_cast<uint32_t>(atoi(entries[0].c_str())));
                limit->set_limit_value(atof(entries[1].c_str()));
            }
            else
            {
                THROW_EXCEPTION(lm::InvalidArgException, "orderParameterReflectingBarrierList was invalid: %s", file.getParameter("orderParameterReflectingBarrierList").c_str());
            }
        }
    }

    // See if any species tracking barriers are set.
    if (file.hasParameter("speciesTrackingBarrierList"))
    {
        vector<string> limitStrings = parseStringList(file.getParameter("speciesTrackingBarrierList"), ",");
        for (size_t i=0; i<limitStrings.size(); i++)
        {
            vector<int32_t> entries = parseIntList(limitStrings[i], ":");
            if (entries.size() == 2 && entries[0] >= 0)
            {
                lm::input::SpeciesLimit* limit = input->mutable_simulation_options()->add_species_tracking_barrier();
                limit->set_species_index(uint32_t(entries[0]));
                limit->set_limit_value(entries[1]);
            }
            else
            {
                THROW_EXCEPTION(lm::InvalidArgException, "speciesTrackingBarrierList was invalid: %s", file.getParameter("speciesTrackingBarrierList").c_str());
            }
        }
    }

    // See if any order parameter tracking barriers are set.
    if (file.hasParameter("orderParameterTrackingBarrierList"))
    {
        vector<string> limitStrings = parseStringList(file.getParameter("orderParameterTrackingBarrierList"), ",");
        for (size_t i=0; i<limitStrings.size(); i++)
        {
            vector<string> entries = parseStringList(limitStrings[i], ":");
            if (entries.size() == 2)
            {
                lm::input::OrderParameterLimit* limit = input->mutable_simulation_options()->add_order_parameter_tracking_barrier();
                limit->set_order_parameter_index(static_cast<uint32_t>(atoi(entries[0].c_str())));
                limit->set_limit_value(atof(entries[1].c_str()));
            }
            else
            {
                THROW_EXCEPTION(lm::InvalidArgException, "orderParameterTrackingBarrierList was invalid: %s", file.getParameter("orderParameterTrackingBarrierList").c_str());
            }
        }
    }

    // See if any tracking barrier limits are set.
    if (file.hasParameter("trackingBarrierLimitList"))
    {
        vector<string> limitStrings = parseStringList(file.getParameter("trackingBarrierLimitList"), ",");
        for (size_t i=0; i<limitStrings.size(); i++)
        {
            vector<int32_t> entries = parseIntList(limitStrings[i], ":");
            if (entries.size() == 2 && entries[0] >= 0 && entries[0] < (input->simulation_options().species_tracking_barrier().size() + input->simulation_options().order_parameter_tracking_barrier().size()))
            {
                lm::input::BarrierLimit* limit = input->mutable_simulation_options()->add_tracking_barrier_crossing_limit();
                limit->set_barrier_index(static_cast<uint32_t>(entries[0]));
                limit->set_limit_value(static_cast<uint32_t>(entries[1]));
            }
            else
            {
                THROW_EXCEPTION(lm::InvalidArgException, "trackingBarrierLimitList was invalid: %s", file.getParameter("trackingBarrierLimitList").c_str());
            }
        }
    }


    // See if the steps per work unit is set.
    if (file.hasParameter("stepsPerWorkUnit")) input->mutable_simulation_options()->set_steps_per_work_unit_part(static_cast<uint64_t>(atoi(file.getParameter("stepsPerWorkUnit").c_str())));
}

void Hdf5InputReader::readOutputOptions(const lm::io::hdf5::Hdf5File& file, lm::input::Input* input)
{
    if (file.hasParameter("writeInterval")) input->mutable_output_options()->set_species_write_interval(atof(file.getParameter("writeInterval").c_str()));
    if (file.hasParameter("latticeWriteInterval")) input->mutable_output_options()->set_lattice_write_interval(atof(file.getParameter("latticeWriteInterval").c_str()));
    if (file.hasParameter("writeInterval")) input->mutable_output_options()->set_concentrations_write_interval(atof(file.getParameter("writeInterval").c_str()));
    if (file.hasParameter("degreeAdvancementWriteInterval")) input->mutable_output_options()->set_degree_advancement_write_interval(atof(file.getParameter("degreeAdvancementWriteInterval").c_str()));
    if (file.hasParameter("orderParameterWriteInterval")) input->mutable_output_options()->set_order_parameter_write_interval(atof(file.getParameter("orderParameterWriteInterval").c_str()));
    if (file.hasParameter("noiseWriteInterval")) input->mutable_output_options()->set_stochastic_process_write_interval(atof(file.getParameter("noiseWriteInterval").c_str()));

    // Add the species FPT tracking list.
    if (file.hasParameter("fptTrackingList"))
    {
        vector<uint32_t> indices = parseUintList(file.getParameter("fptTrackingList"));
        for (size_t i=0; i<indices.size(); i++)
            input->mutable_output_options()->add_fpt_species_to_track(indices[i]);
    }

    // Add the order parameter FPT tracking list.
    if (file.hasParameter("fptOrderParameterTrackingList"))
    {
        vector<uint32_t> indices = parseUintList(file.getParameter("fptOrderParameterTrackingList"));
        for (size_t i=0; i<indices.size(); i++)
            input->mutable_output_options()->add_fpt_order_parameter_to_track(indices[i]);
    }
}

void Hdf5InputReader::readRestart(const lm::io::hdf5::Hdf5File& file, lm::input::Input* input)
{
    if (file.hasReactionRestart())
    {
        file.getReactionRestart(input->mutable_cme_restart());
        if (file.hasDiffusionRestart())
        {
            file.getDiffusionRestart(input->mutable_rdme_restart(), input->cme_restart().restart_trajectory_ids().shape(0));
        }
    }
}

}
}
}
