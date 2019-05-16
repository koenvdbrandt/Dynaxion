/**
 * @file
 * @brief Defines the internal Geant4 geometry construction
 * @copyright Copyright (c) 2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_MODULE_TARGET_CONSTRUCTION_H
#define ALLPIX_MODULE_TARGET_CONSTRUCTION_H

#include <G4LogicalVolume.hh>
#include <memory>
#include <utility>
#include "G4Material.hh"
#include "G4VSolid.hh"
#include "G4LogicalVolume.hh"
#include "core/config/ConfigReader.hpp"
#include "core/geometry/BaseBuilder.hpp"
#include "core/geometry/GeometryManager.hpp"

namespace allpix {
    /**
     * @brief Constructs the Geant4 geometry during Geant4 initialization
     */
    //template<typename WorldVolume, typename Materials>    
    class TargetConstructionG4 : public BaseBuilder<G4LogicalVolume, G4Material> {
    public:
        /**
         * @brief Constructs geometry construction module
         * @param config Configuration object of the geometry builder module
         */
        TargetConstructionG4(Configuration& config);

        /**
         * @brief Constructs the world geometry with all detectors
         * @return Physical volume representing the world
         */
        void build(G4LogicalVolume* world_log, std::map<std::string, G4Material*> materials_) override;

    private:
        Configuration& config_;

        // Storage of internal objects
        std::vector<std::shared_ptr<G4VSolid>> solids_;
    };

} // namespace allpix

#endif /* ALLPIX_MODULE_TARGET_CONSTRUCTION_H */