/**
 * @file
 * @brief Parameters of a box passive material model
 *
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_PASSIVE_MATERIAL_DYNAXION_H
#define ALLPIX_PASSIVE_MATERIAL_DYNAXION_H

#include <string>
#include <utility>

#include <Math/Cartesian2D.h>
#include <Math/DisplacementVector2D.h>
#include <Math/Point3D.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>

#include <G4Box.hh>
#include <G4Tubs.hh>
#include <G4SubtractionSolid.hh>
#include <G4VSolid.hh>
#include "PassiveMaterialModel.hpp"

namespace allpix {

    /**
     * @ingroup PassiveMaterialModel
     * @brief Model of a rectangular box.
     */
    class DynaxionModel : public PassiveMaterialModel {
    public:
        /**
         * @brief Constructs the box passive material model
         * @param config Configuration with description of the model
         */
        explicit DynaxionModel(Configuration& config) : PassiveMaterialModel(config), config_(config) {

            // Set the box specifications
            std::string name = config_.getName();
            auto sub_type = config.get<std::string>("sub_type");
            auto large_tube_inner_radius = config_.get<double>("large_tube_inner_radius", 0);
            auto large_tube_outer_radius = config_.get<double>("large_tube_outer_radius", 0);
            auto large_tube_lenght = config_.get<double>("large_tube_lenght", 0);
            auto small_tube_inner_radius = config_.get<double>("small_tube_inner_radius", 0);
            auto small_tube_outer_radius = config_.get<double>("small_tube_outer_radius", 0);
            auto small_tube_lenght = config_.get<double>("small_tube_lenght", 0);
            auto large_tube_volume = new G4Tubs(name + "large_tube_volume",
                                                large_tube_inner_radius,
                                                large_tube_outer_radius,
                                                large_tube_lenght / 2,
                                                0,
                                                2 * CLHEP::pi);

            auto large_tube_inner_volume = new G4Tubs(
                name + "large_tube_inner_volume", 0, large_tube_inner_radius, large_tube_lenght / 2, 0, 2 * CLHEP::pi);

            auto small_tube_volume = new G4Tubs(name + "small_tube_volume",
                                                small_tube_inner_radius,
                                                small_tube_outer_radius,
                                                small_tube_lenght / 2,
                                                0,
                                                2 * CLHEP::pi);

            auto small_tube_outer_volume = new G4Tubs(
                name + "small_tube_inner_volume", 0, small_tube_outer_radius, small_tube_lenght / 2, 0, 2 * CLHEP::pi); 
            G4RotationMatrix* yRot = new G4RotationMatrix; // Rotates X and Z axes only
            yRot->rotateY(CLHEP::pi / 2. * CLHEP::rad); // Rotates 45 degrees

            if(sub_type == "wrapper") {
                // Create Wrapper for the tube
                solid_ = new G4Box(
                    name + "wrapper_volume", large_tube_lenght / 2, large_tube_outer_radius, small_tube_lenght);

                max_size_ = std::max(large_tube_lenght, std::max(large_tube_outer_radius, 2*small_tube_lenght));
            }
            else if(sub_type == "large_tube") {
                solid_ = new G4SubtractionSolid(
                    name + "large_volume", large_tube_volume, small_tube_outer_volume, yRot, G4ThreeVector(-0.5 * small_tube_lenght, 0, 0));
                max_size_ = std::max(large_tube_lenght, large_tube_outer_radius);
            }
            else if(sub_type == "small_tube") {
                solid_ = new G4SubtractionSolid(name + "small_volume_final",
                                                            small_tube_volume,
                                                            large_tube_inner_volume,
                                                            yRot,
                                                            G4ThreeVector(0, 0, -0.5 * small_tube_lenght));
                max_size_ = std::max(small_tube_lenght, small_tube_outer_radius);

            }

            else if(sub_type == "conveyer_support") {
                auto new_large_tube_volume = new G4Tubs(name + "large_tube_volume",
                                                    large_tube_inner_radius,
                                                    5 * large_tube_outer_radius,
                                                    large_tube_lenght,
                                                    0,
                                                    2 * CLHEP::pi);

                auto support_width = config_.get<double>("support_width", 0);
                auto support_thickness = config_.get<double>("support_thickness", 0);
                auto support_volume =
                    new G4Box(name + "support_volume", large_tube_lenght / 2, support_thickness / 2, support_width / 2);

                auto luggage_height = config_.get<double>("luggage_height", 0);
                auto conveyer_height = config_.get<double>("conveyer_height", 0);
                std::cout << "luggage_height = " << luggage_height << std::endl;
                std::cout << "conveyer_height = " << conveyer_height << std::endl;
                std::cout << "support_thickness = " << support_thickness << std::endl;
                solid_ = new G4SubtractionSolid(
                    name + "_support_volume",
                    support_volume,
                    new_large_tube_volume,
                    yRot,
                    G4ThreeVector(0, +(luggage_height / 2 + conveyer_height + support_thickness / 2), 0));

            max_size_ = std::max(large_tube_lenght, std::max(support_thickness,support_width));
            }
            else if(sub_type == "concrete_sides") {
                auto concrete_size = config_.get<ROOT::Math::XYZVector>("concrete_size", {0., 0., 0.});
                auto concrete_thickness = config_.get<double>("concrete_thickness", 0);
                auto concrete_volume =
                    new G4Box(name + "concrete_volume", concrete_size.x() / 2, concrete_size.y() / 2, concrete_thickness / 2);
                auto large_tube_outer_volume =
                    new G4Tubs(name + "large_tube_volume", 0, large_tube_outer_radius, concrete_size.x(), 0, 2 * CLHEP::pi);
                auto allignment = config_.get<ROOT::Math::XYVector>("allignment", {0., 0.});
                solid_ = new G4SubtractionSolid(name + "concrete_wall_volume",
                                                                                concrete_volume,
                                                                                large_tube_outer_volume,
                                                                                yRot,
                                                                                G4ThreeVector(0, allignment.y(), allignment.x()));
            max_size_ = std::max(concrete_size.x(), std::max(concrete_size.y(),concrete_size.z()));
            }
            else{
                throw ModuleError("Dynaxion sub type '" + sub_type + "' is incorrect.");
            }                                                
        }
        // Set the override functions of PassiveMaterialModel
        G4VSolid* getSolid() override { return solid_; }
        double getMaxSize() override { return max_size_; }

    private:
        Configuration& config_;

        // G4VSolid returnables
        G4VSolid* solid_;

        double max_size_;
    };
} // namespace allpix

#endif /* ALLPIX_PASSIVE_MATERIAL_DYNAXION_H */
