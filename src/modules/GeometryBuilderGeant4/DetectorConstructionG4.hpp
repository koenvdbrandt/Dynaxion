/**
 * @file
 * @brief Defines the internal Geant4 geometry construction
 * @copyright Copyright (c) 2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_MODULE_DETECTOR_CONSTRUCTION_H
#define ALLPIX_MODULE_DETECTOR_CONSTRUCTION_H

#include <memory>
#include <utility>

#include "G4LogicalVolume.hh"
#include "G4Material.hh"
#include "G4VSolid.hh"

#include "G4OpticalSurface.hh"

#include "core/geometry/GeometryBuilder.hpp"

#include "core/geometry/GeometryManager.hpp"

namespace allpix {
    /**
     * @brief Constructs the Geant4 geometry during Geant4 initialization
     */
    class DetectorConstructionG4 : public GeometryBuilder<G4LogicalVolume, G4Material> {
    public:
        /**
         * @brief Constructs geometry construction module
         * @param geo_manager Pointer to the geometry manager, containing the detectors
         * @param config Configuration object of the geometry builder module
         */
        DetectorConstructionG4(GeometryManager* geo_manager, Configuration& config);

        /**
         * @brief Constructs the world geometry with all detectors
         * @return Physical volume representing the world
         */
        void build(G4LogicalVolume* world_log, std::map<std::string, G4Material*> materials_) override;
        std::vector<G4ThreeVector> GetPMTPositions() { return PMT_positions; }

    private:
        // PMT stuff
        void PlacePMTs(G4LogicalVolume* pmt_Log,
                       G4RotationMatrix* rot,
                       G4double& a,
                       G4double& b,
                       G4double da,
                       G4double db,
                       G4double amin,
                       G4double bmin,
                       G4int na,
                       G4int nb,
                       G4double& x,
                       G4double& y,
                       G4double& z,
                       G4int& k);
        // void SurfaceProperties();
        std::shared_ptr<G4LogicalVolume> housing_log;
        std::shared_ptr<G4LogicalVolume> sensor_log;
        std::shared_ptr<G4VPhysicalVolume> sensor_phys;
        std::shared_ptr<G4VPhysicalVolume> housing_phys;
        std::shared_ptr<G4VPhysicalVolume> scint_phys;
        
        std::vector<G4ThreeVector> PMT_positions;
        G4double housing_reflectivity;

        // Logical volumes
        //
        std::shared_ptr<G4LogicalVolume> PMT_log;
        std::shared_ptr<G4LogicalVolume> photo_cath_log;

        GeometryManager* geo_manager_;
        Configuration& config_;

        // Storage of internal objects
        std::vector<std::shared_ptr<G4VSolid>> solids_;
        G4Material* world_material_{};
    };
} // namespace allpix

#endif /* ALLPIX_MODULE_DETECTOR_CONSTRUCTION_H */
