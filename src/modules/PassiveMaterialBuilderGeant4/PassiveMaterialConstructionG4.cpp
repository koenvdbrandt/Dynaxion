/**
 * @file
 * @brief Implements the Geant4 geometry construction process
 * @remarks Code is based on code from Mathieu Benoit
 * @copyright Copyright (c) 2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "PassiveMaterialConstructionG4.hpp"
#include <memory>
#include <string>
#include <utility>

#include <Math/RotationX.h>
#include <Math/RotationY.h>
#include <Math/RotationZ.h>
#include <Math/RotationZYX.h>
#include <Math/Vector3D.h>

#include <G4Box.hh>
#include <G4IntersectionSolid.hh>
#include <G4RotationMatrix.hh>
#include <G4Sphere.hh>
#include <G4SubtractionSolid.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include "CLHEP/Vector/Rotation.h"

#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4Material.hh>
#include <G4PVDivision.hh>
#include <G4PVPlacement.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4StepLimiterPhysics.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4UserLimits.hh>
#include <G4VSolid.hh>
#include <G4VisAttributes.hh>
#include "core/module/exceptions.h"
#include "core/utils/log.h"
#include "tools/ROOT.h"
#include "tools/geant4.h"

#include "GB01BOptrMultiParticleChangeCrossSection.hpp"

using namespace allpix;

PassiveMaterialConstructionG4::PassiveMaterialConstructionG4(Configuration& config, bool bias)
    : config_(config), bias_(bias) {}

/**
 * @brief Version of std::make_shared that does not delete the pointer
 *
 * This version is needed because some pointers are deleted by Geant4 internally, but they are stored as std::shared_ptr in
 * the framework.
 */
template <typename T, typename... Args> static std::shared_ptr<T> make_shared_no_delete(Args... args) {
    return std::shared_ptr<T>(new T(args...), [](T*) {});
}

void PassiveMaterialConstructionG4::build(G4LogicalVolume* world_log, std::map<std::string, G4Material*> materials_) {
    /*
    Get the name of the Passive Material
    */
    std::string name = config_.getName();
    /*
    Get the world_material
    */
    std::string world_material = world_log->GetMaterial()->GetName();
    /*
    Get the information for the passive materials
    */
    auto passive_material_location = config_.get<ROOT::Math::XYZPoint>("position", {0., 0., 0.});
    auto passive_material_pos = toG4Vector(passive_material_location);
    auto passive_material = config_.get<std::string>("material", world_material);

    std::transform(passive_material.begin(), passive_material.end(), passive_material.begin(), ::tolower);
    auto orientation_vector = config_.get<ROOT::Math::XYZVector>("orientation", {0., 0., 0.});

    /*
        // Calculate possible detector misalignment to be added !! Have to figure out the random_generater from
       GeoManager.cpp
        auto misalignment = [&](auto residuals) {
        double dx = std::normal_distribution<double>(0, residuals.x())(random_generator);
        double dy = std::normal_distribution<double>(0, residuals.y())(random_generator);
        double dz = std::normal_distribution<double>(0, residuals.z())(random_generator);
        return DisplacementVector3D<Cartesian3D<double>>(dx, dy, dz);
        };
        orientation_vector += misalignment(config_.get<ROOT::Math::XYZVector>("alignment_precision_orientation", {0., 0.,
       0.}));
    */

    ROOT::Math::Rotation3D orientation;

    auto orientation_mode = config_.get<std::string>("orientation_mode", "xyz");
    if(orientation_mode == "zyx") {
        // First angle given in the configuration file is around z, second around y, last around x:
        LOG(DEBUG) << "Interpreting Euler angles as ZYX rotation";
        orientation = ROOT::Math::RotationZYX(orientation_vector.x(), orientation_vector.y(), orientation_vector.z());
    } else if(orientation_mode == "xyz") {
        LOG(DEBUG) << "Interpreting Euler angles as XYZ rotation";
        // First angle given in the configuration file is around x, second around y, last around z:
        orientation = ROOT::Math::RotationZ(orientation_vector.z()) * ROOT::Math::RotationY(orientation_vector.y()) *
                      ROOT::Math::RotationX(orientation_vector.x());
    } else if(orientation_mode == "zxz") {
        LOG(DEBUG) << "Interpreting Euler angles as ZXZ rotation";
        // First angle given in the configuration file is around z, second around x, last around z:
        orientation = ROOT::Math::EulerAngles(orientation_vector.x(), orientation_vector.y(), orientation_vector.z());
    } else {
        throw InvalidValueError(config_, "orientation_mode", "orientation_mode should be either 'zyx', xyz' or 'zxz'");
    }

    std::vector<double> copy_vec(9);
    orientation.GetComponents(copy_vec.begin(), copy_vec.end());
    ROOT::Math::XYZPoint vx, vy, vz;
    orientation.GetComponents(vx, vy, vz);
    auto rotWrapper = std::make_shared<G4RotationMatrix>(copy_vec.data());
    // auto wrapperGeoTranslation = toG4Vector((0.,0.,0.));
    // wrapperGeoTranslation *= *rotWrapper;
    G4ThreeVector posWrapper = toG4Vector(passive_material_location); // - wrapperGeoTranslation;
    G4Transform3D transform_phys(*rotWrapper, posWrapper);

    if(config_.get<std::string>("type") == "box") {
        auto box_size = config_.get<ROOT::Math::XYVector>("size", {0, 0});
        auto box_thickness = config_.get<double>("thickness", 0);
        auto box_volume = std::make_shared<G4Box>(name + "_volume", box_size.x() / 2, box_size.y() / 2, box_thickness / 2);
        solids_.push_back(box_volume);

        // Place the logical volume of the box
        auto box_log = make_shared_no_delete<G4LogicalVolume>(box_volume.get(), materials_[passive_material], name + "_log");

        // Place the physical volume of the box
        auto box_phys_ =
            make_shared_no_delete<G4PVPlacement>(transform_phys, box_log.get(), name + "_phys", world_log, false, 0, true);
    }

    if(config_.get<std::string>("type") == "cylinder") {
        auto cylinder_inner_radius = config_.get<double>("inner_radius", 0);
        auto cylinder_outer_radius = config_.get<double>("outer_radius", 0);
        auto cylinder_height = config_.get<double>("height", 0);
        auto cylinder_starting_angle = config_.get<double>("starting_angle", 0);
        auto cylinder_arc_length = config_.get<double>("arc_length", 0);
        auto cylinder_volume = std::make_shared<G4Tubs>(name + "_volume",
                                                        cylinder_inner_radius,
                                                        cylinder_outer_radius,
                                                        cylinder_height / 2,
                                                        cylinder_starting_angle * CLHEP::pi,
                                                        cylinder_arc_length * CLHEP::pi);
        solids_.push_back(cylinder_volume);

        // Place the logical volume of the cylinder
        auto cylinder_log =
            make_shared_no_delete<G4LogicalVolume>(cylinder_volume.get(), materials_[passive_material], name + "_log");

        // Place the physical volume of the cylinder
        auto cylinder_phys_ = make_shared_no_delete<G4PVPlacement>(
            transform_phys, cylinder_log.get(), name + "_phys", world_log, false, 0, true);

        auto filling_material = config_.get<std::string>("filling_material", "");

        if(filling_material != "") {
            // if(cylinder_length == 2){
            auto cylinder_filling_volume = std::make_shared<G4Tubs>(name + "_filling_volume",
                                                                    0,
                                                                    cylinder_inner_radius,
                                                                    cylinder_height,
                                                                    cylinder_starting_angle * CLHEP::pi,
                                                                    cylinder_arc_length * CLHEP::pi);
            solids_.push_back(cylinder_filling_volume);

            // Place the logical volume of the filling material
            auto cylinder_filling_log = make_shared_no_delete<G4LogicalVolume>(
                cylinder_filling_volume.get(), materials_[filling_material], name + "_filling_log");

            // Place the physical volume of the filling material
            auto cylinder_filling_phys_ = make_shared_no_delete<G4PVPlacement>(
                transform_phys, cylinder_filling_log.get(), name + "_filling_phys", world_log, false, 0, true);
            //}
            // else{ throw ModuleError("Cylinder '" + name + "' is not closed! Can't fill it with material");}
        }
    }

    if(config_.get<std::string>("type") == "tube") {
        auto tube_outer_diameter = config_.get<ROOT::Math::XYVector>("outer_diameter", {0, 0});
        auto tube_inner_diameter = config_.get<ROOT::Math::XYVector>("inner_diameter", {0, 0});
        auto tube_length = config_.get<double>("length", 0);
        if(tube_inner_diameter.x() >= tube_outer_diameter.x() || tube_inner_diameter.y() >= tube_outer_diameter.y()) {
            throw ModuleError("Inner diameter of '" + name +
                              "' is larget than its outer diameter! Can't construct the tube");
        }

        auto tube_outer_volume =
            new G4Box(name + "_outer_volume", tube_outer_diameter.x() / 2, tube_outer_diameter.y() / 2, tube_length / 2);

        auto tube_inner_volume = new G4Box(
            name + "_inner_volume", tube_inner_diameter.x() / 2, tube_inner_diameter.y() / 2, 1.1 * tube_length / 2);

        auto tube_final_volume =
            std::make_shared<G4SubtractionSolid>(name + "_final_volume", tube_outer_volume, tube_inner_volume);
        solids_.push_back(tube_final_volume);

        // Place the logical volume of the tube
        auto tube_log =
            make_shared_no_delete<G4LogicalVolume>(tube_final_volume.get(), materials_[passive_material], name + "_log");

        // Place the physical volume of the tube
        auto tube_phys_ =
            make_shared_no_delete<G4PVPlacement>(transform_phys, tube_log.get(), name + "_phys", world_log, false, 0, true);

        // Fill the material with the filling_material if needed
        auto filling_material = config_.get<std::string>("filling_material", "");
        if(filling_material != "") {
            auto tube_filling_volume = std::make_shared<G4Box>(
                name + "_filling_volume", tube_inner_diameter.x() / 2, tube_inner_diameter.y() / 2, tube_length / 2);
            solids_.push_back(tube_filling_volume);
            auto tube_filling_log = make_shared_no_delete<G4LogicalVolume>(
                tube_filling_volume.get(), materials_[filling_material], name + "filling_log");
            auto tube_filling_phys_ = make_shared_no_delete<G4PVPlacement>(
                transform_phys, tube_filling_log.get(), name + "filling_phys", world_log, false, 0, true);
        }
    }

    if(config_.get<std::string>("type") == "sphere") {
        auto sphere_inner_radius = config_.get<double>("inner_radius", 0);
        auto sphere_outer_radius = config_.get<double>("outer_radius", 0);
        auto sphere_starting_angle_phi = config_.get<double>("starting_angle_phi", 0);
        auto sphere_arc_length_phi = config_.get<double>("arc_length_phi", 2);
        auto sphere_starting_angle_theta = config_.get<double>("starting_angle_theta", 0);
        auto sphere_arc_length_theta = config_.get<double>("arc_length_theta", 1);
        auto sphere_volume = std::make_shared<G4Sphere>(name + "_volume",
                                                        sphere_inner_radius,
                                                        sphere_outer_radius,
                                                        sphere_starting_angle_phi * CLHEP::pi,
                                                        sphere_arc_length_phi * CLHEP::pi,
                                                        sphere_starting_angle_theta * CLHEP::pi,
                                                        sphere_arc_length_theta * CLHEP::pi);
        solids_.push_back(sphere_volume);

        // Place the logical volume of the sphere
        auto sphere_log =
            make_shared_no_delete<G4LogicalVolume>(sphere_volume.get(), materials_[passive_material], name + "_log");

        // Place the physical volume of the sphere
        auto sphere_phys_ = make_shared_no_delete<G4PVPlacement>(
            transform_phys, sphere_log.get(), name + "_phys", world_log, false, 0, true);

        auto filling_material = config_.get<std::string>("filling_material", "");

        if(filling_material != "") {
            auto sphere_filling_volume = std::make_shared<G4Sphere>(name + "_filling_volume",
                                                                    0,
                                                                    sphere_inner_radius,
                                                                    sphere_starting_angle_phi * CLHEP::pi,
                                                                    sphere_arc_length_phi * CLHEP::pi,
                                                                    sphere_starting_angle_theta * CLHEP::pi,
                                                                    sphere_arc_length_theta * CLHEP::pi);
            solids_.push_back(sphere_filling_volume);

            // Place the logical volume of the filling material
            auto sphere_filling_log = make_shared_no_delete<G4LogicalVolume>(
                sphere_filling_volume.get(), materials_[filling_material], name + "_filling_log");

            // Place the physical volume of the filling material
            auto sphere_filling_phys_ = make_shared_no_delete<G4PVPlacement>(
                transform_phys, sphere_filling_log.get(), name + "_filling_phys", world_log, false, 0, true);
        }
    }
    if(config_.get<std::string>("type") == "dynaxion_tube") {

        auto large_tube_inner_radius = config_.get<double>("large_tube_inner_radius", 0);
        auto large_tube_outer_radius = config_.get<double>("large_tube_outer_radius", 0);
        auto large_tube_lenght = config_.get<double>("large_tube_lenght", 0);
        auto small_tube_inner_radius = config_.get<double>("small_tube_inner_radius", 0);
        auto small_tube_outer_radius = config_.get<double>("small_tube_outer_radius", 0);
        auto small_tube_lenght = config_.get<double>("small_tube_lenght", 0);

        // Create Wrapper for the tube
        auto dynaxion_tube_wrapper = std::make_shared<G4Box>(
            name + "wrapper_volume", large_tube_lenght / 2, large_tube_outer_radius, small_tube_lenght / 2);
        solids_.push_back(dynaxion_tube_wrapper);

        // Place the logical volume of the box
        auto dynaxion_tube_wrapper_log =
            make_shared_no_delete<G4LogicalVolume>(dynaxion_tube_wrapper.get(), materials_["vacuum"], "Dynaxion_log");

        // Place the physical volume of the box
        // G4ThreeVector posDynaxionWrapper = toG4Vector(ROOT::Math::XYZVector(passive_material_location.x(),
        // passive_material_location.y(), passive_material_location.z() + 0.25*small_tube_lenght));
        // G4Transform3D transform_phys_Dynaxion(*rotWrapper, posDynaxionWrapper);
        auto dynaxion_tube_wrapper_phys_ = make_shared_no_delete<G4PVPlacement>(
            transform_phys, dynaxion_tube_wrapper_log.get(), name + "wrapper_phys", world_log, false, 0, true);

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
        // yRot = orientation;
        yRot->rotateY(CLHEP::pi / 2. * CLHEP::rad); // Rotates 45 degrees
        // yRot *= *orientation;
        auto dynaxion_tube_large_volume = std::make_shared<G4SubtractionSolid>(
            name + "large_volume", large_tube_volume, small_tube_outer_volume, yRot, G4ThreeVector(0, 0, 0));
        solids_.push_back(dynaxion_tube_large_volume);

        // Place the logical volume of the dynaxion large tube
        auto dynaxion_large_volume_log = make_shared_no_delete<G4LogicalVolume>(
            dynaxion_tube_large_volume.get(), materials_[passive_material], name + "large_log");

        // Place the physical volume of the dynaxion large tube
        auto dynaxion_large_volume_phys =
            make_shared_no_delete<G4PVPlacement>(yRot,
                                                 G4ThreeVector(0, 0, -0. /*25 * small_tube_lenght*/),
                                                 dynaxion_large_volume_log.get(),
                                                 name + "large_phys",
                                                 dynaxion_tube_wrapper_log.get(),
                                                 false,
                                                 0,
                                                 true);

        auto dynaxion_tube_small_volume_final =
            std::make_shared<G4SubtractionSolid>(name + "small_volume_final",
                                                 small_tube_volume,
                                                 large_tube_inner_volume,
                                                 yRot,
                                                 G4ThreeVector(0, 0, -0.25 * small_tube_lenght));
        solids_.push_back(dynaxion_tube_small_volume_final);

        // Place the logical volume of the dynaxion large tube
        auto dynaxion_small_volume_log = make_shared_no_delete<G4LogicalVolume>(
            dynaxion_tube_small_volume_final.get(), materials_[passive_material], name + "small_log");

        // Place the physical volume of the dynaxion large tube
        auto dynaxion_tube_small_volume_phys =
            make_shared_no_delete<G4PVPlacement>(nullptr,
                                                 G4ThreeVector(0, 0, +0.25 * small_tube_lenght),
                                                 dynaxion_small_volume_log.get(),
                                                 name + "small_phys",
                                                 dynaxion_tube_wrapper_log.get(),
                                                 false,
                                                 0,
                                                 true);
    }

    if(config_.get<std::string>("type") == "dynaxion_conveyer_support") {

        auto large_tube_inner_radius = config_.get<double>("large_tube_inner_radius", 0);
        auto large_tube_outer_radius = config_.get<double>("large_tube_outer_radius", 0);
        auto large_tube_lenght = config_.get<double>("large_tube_lenght", 0);
        auto large_tube_volume = new G4Tubs(name + "large_tube_volume",
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
        auto conveyer_height = config_.get<double>("conveyer_height ", 0);

        G4RotationMatrix* yRot = new G4RotationMatrix; // Rotates X and Z axes only
        yRot->rotateY(CLHEP::pi / 2. * CLHEP::rad);    // Rotates 45 degrees
        auto dynaxion_conveyer_support_volume = std::make_shared<G4SubtractionSolid>(
            name + "_support_volume",
            support_volume,
            large_tube_volume,
            yRot,
            G4ThreeVector(0, +(luggage_height / 2 + conveyer_height + support_thickness), 0));
        solids_.push_back(dynaxion_conveyer_support_volume);

        // Place the logical volume of the dynaxion_conveyer_support
        auto dynaxion_conveyer_support_volume_log = make_shared_no_delete<G4LogicalVolume>(
            dynaxion_conveyer_support_volume.get(), materials_[passive_material], name + "dynaxion_conveyer_support_log");

        // Place the physical volume of the dynaxion large tube
        G4ThreeVector posDynaxionSupport = toG4Vector(
            ROOT::Math::XYZVector(passive_material_location.x(),
                                  passive_material_location.y() - (luggage_height / 2 + conveyer_height + support_thickness),
                                  passive_material_location.z()));
        G4Transform3D transform_phys_Support(*rotWrapper, posDynaxionSupport);
        auto dynaxion_conveyer_support_volume_phys =
            make_shared_no_delete<G4PVPlacement>(transform_phys_Support,
                                                 dynaxion_conveyer_support_volume_log.get(),
                                                 name + "dynaxion_conveyer_support_phys",
                                                 world_log,
                                                 false,
                                                 0,
                                                 true);
    }
    if(config_.get<std::string>("type") == "concrete_sides") {
        auto concrete_size = config_.get<ROOT::Math::XYVector>("concrete_size", {0., 0.});
        auto concrete_thickness = config_.get<double>("concrete_thickness", 0);
        auto concrete_volume =
            new G4Box(name + "concrete_volume", concrete_size.x() / 2, concrete_size.y() / 2, concrete_thickness / 2);

        auto large_tube_outer_radius = config_.get<double>("large_tube_outer_radius", 0);
        auto large_tube_volume =
            new G4Tubs(name + "large_tube_volume", 0, large_tube_outer_radius, concrete_size.x(), 0, 2 * CLHEP::pi);
        G4RotationMatrix* yRot = new G4RotationMatrix; // Rotates X and Z axes only
        yRot->rotateY(CLHEP::pi / 2. * CLHEP::rad);    // Rotates 45 degrees
        auto allignment = config_.get<ROOT::Math::XYVector>("allignment", {0., 0.});
        auto concrete_wall_volume = std::make_shared<G4SubtractionSolid>(name + "concrete_wall_volume",
                                                                         concrete_volume,
                                                                         large_tube_volume,
                                                                         yRot,
                                                                         G4ThreeVector(0, allignment.y(), allignment.x()));
        solids_.push_back(concrete_wall_volume);

        // Place the logical volume of the dynaxion_conveyer_support
        auto concrete_wall_volume_log = make_shared_no_delete<G4LogicalVolume>(
            concrete_wall_volume.get(), materials_[passive_material], name + "concrete_wall_volume_log");

        // Place the physical volume of the dynaxion large tube
        auto concrete_wall_volume_phys = make_shared_no_delete<G4PVPlacement>(
            transform_phys, concrete_wall_volume_log.get(), name + "concrete_wall_volume_phys", world_log, false, 0, true);
    }
    if(bias_ == true) {
        LOG(TRACE) << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX";
        ConstructSDandField();
    }
}

std::vector<ROOT::Math::XYZPoint> PassiveMaterialConstructionG4::addPoints() {
    ROOT::Math::XYZPoint passive_material_location = config_.get<ROOT::Math::XYZPoint>("position", {0., 0., 0.});
    std::array<int, 8> offset_x = {{1, 1, 1, 1, -1, -1, -1, -1}};
    std::array<int, 8> offset_y = {{1, 1, -1, -1, 1, 1, -1, -1}};
    std::array<int, 8> offset_z = {{1, -1, 1, -1, 1, -1, 1, -1}};
    if(config_.get<std::string>("type") == "box") {
        ROOT::Math::XYVector box_size = config_.get<ROOT::Math::XYVector>("size", {0, 0});
        auto box_thickness = config_.get<double>("thickness", 0);
        for(size_t i = 0; i < 8; ++i) {
            points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + offset_x.at(i) * box_size.x() / 2,
                                                      passive_material_location.y() + offset_y.at(i) * box_size.y() / 2,
                                                      passive_material_location.z() + offset_z.at(i) * box_thickness / 2));
        }
    }
    if(config_.get<std::string>("type") == "tube") {
        auto tube_outer_diameter = config_.get<ROOT::Math::XYVector>("outer_diameter", {0, 0});
        auto tube_length = config_.get<double>("length", 0);
        for(size_t i = 0; i < 8; ++i) {
            points_.emplace_back(
                ROOT::Math::XYZPoint(passive_material_location.x() + offset_x.at(i) * tube_outer_diameter.x() / 2,
                                     passive_material_location.y() + offset_y.at(i) * tube_outer_diameter.y() / 2,
                                     passive_material_location.z() + offset_z.at(i) * tube_length / 2));
        }
    }

    if(config_.get<std::string>("type") == "cylinder") {
        auto cylinder_outer_radius = config_.get<double>("outer_radius", 0);
        auto cylinder_height = config_.get<double>("height", 0);
        for(size_t i = 0; i < 8; ++i) {
            points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + offset_x.at(i) * cylinder_outer_radius,
                                                      passive_material_location.y() + offset_y.at(i) * cylinder_outer_radius,
                                                      passive_material_location.z() + offset_z.at(i) * cylinder_height / 2));
        }
    }

    if(config_.get<std::string>("type") == "sphere") {
        auto sphere_outer_radius = config_.get<double>("outer_radius", 0);
        for(size_t i = 0; i < 8; ++i) {
            points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + offset_x.at(i) * sphere_outer_radius,
                                                      passive_material_location.y() + offset_y.at(i) * sphere_outer_radius,
                                                      passive_material_location.z() + offset_z.at(i) * sphere_outer_radius));
        }
    }
    if(config_.get<std::string>("type") == "dynaxion_tube") {
        auto large_tube_outer_radius = config_.get<double>("large_tube_outer_radius", 0);
        auto large_tube_lenght = config_.get<double>("large_tube_lenght", 0);
        auto small_tube_lenght = config_.get<double>("small_tube_lenght", 0);
        for(size_t i = 0; i < 8; ++i) {
            points_.emplace_back(ROOT::Math::XYZPoint(
                passive_material_location.x() + offset_x.at(i) * large_tube_lenght,
                passive_material_location.y() + offset_y.at(i) * large_tube_outer_radius,
                passive_material_location.z() + 0.5 * small_tube_lenght + offset_z.at(i) * small_tube_lenght));
        }
    }
    if(config_.get<std::string>("type") == "dynaxion_conveyer_support") {
        auto large_tube_lenght = config_.get<double>("large_tube_lenght", 0);
        auto support_width = config_.get<double>("support_width", 0);
        auto support_thickness = config_.get<double>("support_thickness", 0);
        auto luggage_height = config_.get<double>("luggage_height", 0);
        auto conveyer_height = config_.get<double>("conveyer_height ", 0);
        for(size_t i = 0; i < 8; ++i) {
            points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + offset_x.at(i) * large_tube_lenght,
                                                      passive_material_location.y() -
                                                          (luggage_height / 2 + conveyer_height / 2) +
                                                          offset_y.at(i) * support_thickness,
                                                      passive_material_location.z() + offset_z.at(i) * support_width));
        }
    }
    if(config_.get<std::string>("type") == "concrete_sides") {
        auto concrete_size = config_.get<ROOT::Math::XYVector>("concrete_size", {0., 0.});
        auto concrete_thickness = config_.get<double>("concrete_thickness", 0);

        for(size_t i = 0; i < 8; ++i) {
            points_.emplace_back(
                ROOT::Math::XYZPoint(passive_material_location.x() + offset_x.at(i) * concrete_size.x() / 2,
                                     passive_material_location.y() + offset_y.at(i) * concrete_size.y() / 2,
                                     passive_material_location.z() + offset_z.at(i) * concrete_thickness / 2));
        }
    }

    return points_;
}

void PassiveMaterialConstructionG4::ConstructSDandField() {
    std::string name = config_.getName();

    // -- Fetch volume for biasing:
    G4LogicalVolume* logicTest = G4LogicalVolumeStore::GetInstance()->GetVolume(name + "_log");

    // ----------------------------------------------
    // -- operator creation and attachment to volume:
    // ----------------------------------------------
    auto testMany = new GB01BOptrMultiParticleChangeCrossSection();
    LOG(TRACE) << "Adding Particles";

    testMany->AddParticle("deuteron");
    testMany->AddParticle("neutron");
    testMany->AttachTo(logicTest);
    LOG(TRACE) << " Attaching biasing operator " << testMany->GetName() << " to logical volume " << logicTest->GetName();
}
