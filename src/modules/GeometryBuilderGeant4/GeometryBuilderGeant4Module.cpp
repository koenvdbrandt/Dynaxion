/**
 * @file
 * @brief Implementation of Geant4 geometry construction module
 * @copyright Copyright (c) 2017-2020 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "GeometryBuilderGeant4Module.hpp"

#include <cassert>
#include <memory>
#include <string>
#include <utility>

#include <G4RunManager.hh>
#include <G4UImanager.hh>
#include <G4UIterminal.hh>
#include <G4Version.hh>
#include <G4VisManager.hh>

#include <Math/Vector3D.h>

#include "GeometryConstructionG4.hpp"

#include "tools/ROOT.h"
#include "tools/geant4.h"

#include "DetectorConstructionG4.hpp"
#include "PassiveMaterialConstructionG4.hpp"

#include "core/config/ConfigReader.hpp"
#include "core/config/exceptions.h"
#include "core/geometry/GeometryManager.hpp"
#include "core/utils/log.h"

using namespace allpix;
using namespace ROOT;

GeometryBuilderGeant4Module::GeometryBuilderGeant4Module(Configuration& config, Messenger*, GeometryManager* geo_manager)
    : Module(config), geo_manager_(geo_manager), run_manager_g4_(nullptr) {
    // Check for passive materials
    bool passive_materials = config_.get<bool>("passive_materials", false);

    if(passive_materials) {
        // Reading passive material file
        std::string passive_material_file_name = config_.getPath("passive_materials_file", true);
        LOG(TRACE) << "Reading passive material configuration";

        std::ifstream passive_material_file(passive_material_file_name);
        ConfigReader passive_material_reader(passive_material_file, passive_material_file_name);
        passive_material_configs_ = passive_material_reader.getConfigurations();

        std::set<std::string> passive_material_names;
        for(auto& passive_material_section : passive_material_configs_) {
            auto name = passive_material_section.getName();
            if(passive_material_names.find(name) != passive_material_names.end()) {
                throw ModuleError("Passive Material with name '" + name +
                                  "' is already registered, Passive Material names should be unique");
            }
            passive_material_names.insert(name);

            // Add the min and max points to the world volume
            for(auto& point : PassiveMaterialConstructionG4(passive_material_section, geo_manager_).addPoints()) {
                LOG(TRACE) << "adding point " << Units::display(point, {"mm", "um"}) << "to the geometry";
                geo_manager_->addPoint(point);
            }
        }
    }
}

/**
 * @brief Checks if a particular Geant4 dataset is available in the environment
 * @throws ModuleError If a certain Geant4 dataset is not set or not available
 */
static void check_dataset_g4(const std::string& env_name) {
    const char* file_name = std::getenv(env_name.c_str());
    if(file_name == nullptr) {
        throw ModuleError("Geant4 environment variable " + env_name +
                          " is not set, make sure to source a Geant4 "
                          "environment with all datasets");
    }
    std::ifstream file(file_name);
    if(!file.good()) {
        throw ModuleError("Geant4 environment variable " + env_name +
                          " does not point to existing dataset, the Geant4 "
                          "environment is invalid");
    }
    // FIXME: check if file does actually contain a correct dataset
}

void GeometryBuilderGeant4Module::init() {
    // Check if all the required geant4 datasets are defined
    LOG(DEBUG) << "Checking Geant4 datasets";
    check_dataset_g4("G4LEVELGAMMADATA");
    check_dataset_g4("G4RADIOACTIVEDATA");
    check_dataset_g4("G4PIIDATA");
    check_dataset_g4("G4SAIDXSDATA");
    check_dataset_g4("G4ABLADATA");
    check_dataset_g4("G4REALSURFACEDATA");
    check_dataset_g4("G4NEUTRONHPDATA");
    check_dataset_g4("G4ENSDFSTATEDATA");
    check_dataset_g4("G4LEDATA");

// Check for Neutron XS data only for Geant4 version prior to 10.5, deprecated dataset from 10.5
#if G4VERSION_NUMBER < 1050
    check_dataset_g4("G4NEUTRONXSDATA");
#endif

    // Suppress all output (also stdout due to a part in Geant4 where G4cout is not used)
    SUPPRESS_STREAM(std::cout);
    SUPPRESS_STREAM(G4cout);

    // Create the G4 run manager
    run_manager_g4_ = std::make_unique<G4RunManager>();

    // Release stdout again
    RELEASE_STREAM(std::cout);

    // Set the geometry construction to use
    auto geometry_construction = new GeometryConstructionG4(geo_manager_, config_, passive_material_configs_);
    run_manager_g4_->SetUserInitialization(geometry_construction);

    // Run the geometry construct function in GeometryConstructionG4
    LOG(TRACE) << "Building Geant4 geometry";
    run_manager_g4_->InitializeGeometry();

    // Release output from G4
    RELEASE_STREAM(G4cout);
}
