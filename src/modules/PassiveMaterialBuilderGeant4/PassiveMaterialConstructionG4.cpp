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
#include <G4NistManager.hh>
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

using namespace allpix;

PassiveMaterialConstructionG4::PassiveMaterialConstructionG4(Configuration& config) : config_(config) {}
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
    auto world_material = world_log->GetMaterial()->GetName();
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
                                                        cylinder_height,
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
}

std::vector<ROOT::Math::XYZPoint> PassiveMaterialConstructionG4::addPoints() {
    ROOT::Math::XYZPoint passive_material_location = config_.get<ROOT::Math::XYZPoint>("position", {0., 0., 0.});

    if(config_.get<std::string>("type") == "box") {
        ROOT::Math::XYVector box_size = config_.get<ROOT::Math::XYVector>("size", {0, 0});
        auto box_thickness = config_.get<double>("thickness", 0);

        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + box_size.x() / 2,
                                                  passive_material_location.y() + box_size.y() / 2,
                                                  passive_material_location.z() + box_thickness / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + box_size.x() / 2,
                                                  passive_material_location.y() + box_size.y() / 2,
                                                  passive_material_location.z() - box_thickness / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + box_size.x() / 2,
                                                  passive_material_location.y() - box_size.y() / 2,
                                                  passive_material_location.z() + box_thickness / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + box_size.x() / 2,
                                                  passive_material_location.y() - box_size.y() / 2,
                                                  passive_material_location.z() - box_thickness / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() - box_size.x() / 2,
                                                  passive_material_location.y() + box_size.y() / 2,
                                                  passive_material_location.z() + box_thickness / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() - box_size.x() / 2,
                                                  passive_material_location.y() + box_size.y() / 2,
                                                  passive_material_location.z() - box_thickness / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() - box_size.x() / 2,
                                                  passive_material_location.y() - box_size.y() / 2,
                                                  passive_material_location.z() + box_thickness / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() - box_size.x() / 2,
                                                  passive_material_location.y() - box_size.y() / 2,
                                                  passive_material_location.z() - box_thickness / 2));
    }

    if(config_.get<std::string>("type") == "cylinder") {
        auto cylinder_outer_radius = config_.get<double>("outer_radius", 0);
        auto cylinder_height = config_.get<double>("height", 0);

        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + cylinder_outer_radius,
                                                  passive_material_location.y() + cylinder_outer_radius,
                                                  passive_material_location.z() + cylinder_height / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + cylinder_outer_radius,
                                                  passive_material_location.y() + cylinder_outer_radius,
                                                  passive_material_location.z() - cylinder_height / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + cylinder_outer_radius,
                                                  passive_material_location.y() - cylinder_outer_radius,
                                                  passive_material_location.z() + cylinder_height / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() + cylinder_outer_radius,
                                                  passive_material_location.y() - cylinder_outer_radius,
                                                  passive_material_location.z() - cylinder_height / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() - cylinder_outer_radius,
                                                  passive_material_location.y() + cylinder_outer_radius,
                                                  passive_material_location.z() + cylinder_height / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() - cylinder_outer_radius,
                                                  passive_material_location.y() + cylinder_outer_radius,
                                                  passive_material_location.z() - cylinder_height / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() - cylinder_outer_radius,
                                                  passive_material_location.y() - cylinder_outer_radius,
                                                  passive_material_location.z() + cylinder_height / 2));
        points_.emplace_back(ROOT::Math::XYZPoint(passive_material_location.x() - cylinder_outer_radius,
                                                  passive_material_location.y() - cylinder_outer_radius,
                                                  passive_material_location.z() - cylinder_height / 2));
    }

    return points_;
}
