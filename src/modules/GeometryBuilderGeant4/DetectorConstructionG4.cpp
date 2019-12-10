/**
 * @file
 * @brief Implements the Geant4 geometry construction process
 * @remarks Code is based on code from Mathieu Benoit
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "DetectorConstructionG4.hpp"

#include <memory>
#include <string>
#include <utility>

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4LogicalVolumeStore.hh>
#include <G4NistManager.hh>
#include <G4PVDivision.hh>
#include <G4PVPlacement.hh>
#include <G4PhysicalVolumeStore.hh>
#include <G4Sphere.hh>
#include <G4StepLimiterPhysics.hh>
#include <G4SubtractionSolid.hh>
#include <G4ThreeVector.hh>
#include <G4Tubs.hh>
#include <G4UnionSolid.hh>
#include <G4UserLimits.hh>
#include <G4VSolid.hh>
#include <G4VisAttributes.hh>
#include "G4Material.hh"

#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4MultiFunctionalDetector.hh"
#include "G4OpticalSurface.hh"
#include "G4SurfaceProperty.hh"
#include "G4SystemOfUnits.hh"
#include "G4UImanager.hh"
#include "G4UnitsTable.hh"
#include "G4VPrimitiveScorer.hh"

#include "core/geometry/HybridPixelDetectorModel.hpp"
#include "core/geometry/ScintillatorModel.hpp"
#include "core/module/exceptions.h"
#include "core/utils/log.h"
#include "tools/ROOT.h"
#include "tools/geant4.h"

#include "GeometryConstructionG4.hpp"
#include "Parameterization2DG4.hpp"

using namespace allpix;

DetectorConstructionG4::DetectorConstructionG4(GeometryManager* geo_manager) : geo_manager_(geo_manager) {}

/**
 * @brief Version of std::make_shared that does not delete the pointer
 *
 * This version is needed because some pointers are deleted by Geant4 internally, but they are stored as std::shared_ptr in
 * the framework.
 */
template <typename T, typename... Args> static std::shared_ptr<T> make_shared_no_delete(Args... args) {
    return std::shared_ptr<T>(new T(args...), [](T*) {});
}

void GeometryConstructionG4::build_detectors() {
    // Loop through all detectors and construct them
    std::vector<std::shared_ptr<Detector>> detectors = geo_manager_->getDetectors();
    LOG(TRACE) << "Building " << detectors.size() << " device(s)";

    for(auto& detector : detectors) {
        // Get pointer to the model of the detector
        auto model = detector->getModel();

        std::string name = detector->getName();
        auto sensor_material = model->getSensorMaterial();
        LOG(DEBUG) << "Creating Geant4 model for " << name;
        LOG(DEBUG) << " Wrapper dimensions of model: " << Units::display(model->getSize(), {"mm", "um"});
        LOG(TRACE) << " Sensor dimensions: " << Units::display(model->getSensorSize(), {"mm", "um"});
        LOG(TRACE) << " Sensor material: " << sensor_material;
        LOG(DEBUG) << " Global position and orientation of the detector:";

        // Get position and orientation
        auto position = detector->getPosition();
        LOG(DEBUG) << " - Position\t\t:\t" << Units::display(position, {"mm", "um"});
        ROOT::Math::Rotation3D orientation = detector->getOrientation();
        std::vector<double> copy_vec(9);
        orientation.GetComponents(copy_vec.begin(), copy_vec.end());
        ROOT::Math::XYZPoint vx, vy, vz;
        orientation.GetComponents(vx, vy, vz);
        auto rotWrapper = std::make_shared<G4RotationMatrix>(copy_vec.data());
        auto wrapperGeoTranslation = toG4Vector(model->getCenter() - model->getGeometricalCenter());
        wrapperGeoTranslation *= *rotWrapper;
        G4ThreeVector posWrapper = toG4Vector(position) - wrapperGeoTranslation;
        detector->setExternalObject("rotation_matrix", rotWrapper);
        G4Transform3D transform_phys(*rotWrapper, posWrapper);

        auto scint_model = std::dynamic_pointer_cast<ScintillatorModel>(model);
        if(scint_model != nullptr) {
            /*  Scintillator
                    Creates a housing box with a scintillator and a sensor in it
                    Support layer is optional
                 */

            // Get parameters from model
            auto sensor_size = scint_model->getSensorSize();
            auto scint_shape = scint_model->getScintShape();
            auto scint_size = scint_model->getScintSize();
            auto scint_material = scint_model->getScintMaterial();
            auto housing_material = scint_model->getHousingMaterial();
            auto housing_thickness = scint_model->getHousingThickness();
            auto housing_reflectivity = scint_model->getHousingReflectivity();
            auto housing_surface_model = static_cast<G4OpticalSurfaceModel>(scint_model->getHousingSurfaceModel());
            auto housing_surface_type = static_cast<G4SurfaceType>(scint_model->getHousingSurfaceType());
            auto housing_surface_finish = static_cast<G4OpticalSurfaceFinish>(scint_model->getHousingSurfaceFinish());
            auto housing_surface_value = scint_model->getHousingSurfaceValue();
            auto photocathode_surface_model = static_cast<G4OpticalSurfaceModel>(scint_model->getPhotocathodeSurfaceModel());
            auto photocathode_surface_type = static_cast<G4SurfaceType>(scint_model->getPhotocathodeSurfaceType());
            auto photocathode_surface_finish =
                static_cast<G4OpticalSurfaceFinish>(scint_model->getPhotocathodeSurfaceFinish());
            auto photocathode_surface_value = scint_model->getHousingSurfaceValue();

            // Check is scintillator has the correct properties
            auto scint_prop_table = materials_[scint_material]->GetMaterialPropertiesTable();
            if(scint_prop_table == nullptr) {
                throw ModuleError("Cannot construct a scintillator of material '" + scint_material +
                                  "'. Material doesn't have a property table");
            } else if(scint_prop_table->GetProperty("FASTCOMPONENT") == nullptr ||
                      scint_prop_table->GetProperty("ABSLENGTH") == nullptr ||
                      scint_prop_table->GetProperty("RINDEX") == nullptr ||
                      !scint_prop_table->ConstPropertyExists("SCINTILLATIONYIELD")) {
                throw ModuleError("Cannot construct a scintillator of material '" + scint_material +
                                  "'. Material misses one of the following Material Properties: 'FASTCOMPONENT', "
                                  "'ABSLENGTH', 'RINDEX' , 'SCINTILLATIONYIELD' ");
            }

            // Check if the optical surface options match
            if(housing_surface_model > 1 || photocathode_surface_model > 1) {
                throw ModuleError("Optical surface model must be either 0 or 1");
            } else if(housing_surface_type > 2 || photocathode_surface_type > 2) {
                throw ModuleError("Optical surface type must be 0,1 or 2");
            } else if(housing_surface_finish > 29 || photocathode_surface_finish > 29) {
                throw ModuleError("Optical surface can't be larger than 29");
            }

            else if((housing_surface_type == 2 && housing_surface_finish < 6) ||
                    (photocathode_surface_type == 2 && photocathode_surface_finish < 6)) {
                throw ModuleError("For optical surface type Dielectric-LUT the surface finish must be between 6 and 29");
            } else if((housing_surface_type != 2 && housing_surface_finish > 5) ||
                      (photocathode_surface_type != 2 && photocathode_surface_finish > 5)) {
                throw ModuleError("For optical surface type Dielectric-Metal and Dielectric-Dielectric the surface finish "
                                  "must be between 0 and 5");
            }

            // Mode specific parameters
            LOG(TRACE) << " Scintillator Model Parameters for detector " << name;
            LOG(TRACE) << " Scintillator dimensions: " << Units::display(scint_size, {"mm", "um"});
            LOG(TRACE) << " Scintillator material: " << scint_material;
            LOG(TRACE) << " Housing dimensions: "
                       << Units::display(
                              ROOT::Math::XYZVector(scint_size.x() / 2.0 + housing_thickness,
                                                    scint_size.y() / 2.0 + housing_thickness,
                                                    scint_size.z() / 2.0 + housing_thickness + sensor_size.z() / 2.0),
                              {"mm", "um"});
            LOG(TRACE) << " Housing material: " << housing_material;
            LOG(DEBUG) << " Center of the geometry parts relative to the detector wrapper geometric center:";

            // Create the solids of the housing and the scintillator
            if(scint_shape == "box") {
                housing_solid_ = std::make_shared<G4Box>("housing_" + name + "_solid",
                                                         scint_size.x() / 2.0 + housing_thickness,
                                                         scint_size.y() / 2.0 + housing_thickness,
                                                         scint_size.z() / 2.0 + housing_thickness + sensor_size.z() / 2.0);
                scint_solid_ = std::make_shared<G4Box>(
                    "scint_" + name, scint_size.x() / 2.0, scint_size.y() / 2.0, scint_size.z() / 2.0);
            } else if(scint_shape == "cylinder") {
                housing_solid_ = std::make_shared<G4Tubs>("housing_" + name + "_solid",
                                                          0,
                                                          scint_size.x() / 2.0 + housing_thickness,
                                                          scint_size.z() / 2.0 + sensor_size.z() / 2.0 + housing_thickness,
                                                          0,
                                                          2 * CLHEP::pi);
                scint_solid_ = std::make_shared<G4Tubs>(
                    "scint_" + name, 0, scint_size.x() / 2.0, scint_size.z() / 2.0, 0, 2 * CLHEP::pi);
            }
            solids_.push_back(housing_solid_);
            solids_.push_back(scint_solid_);

            /*   Housing
                   the housing of the scintillator
                   housing works like the wrapper. There was some issue whith placing the housing in a wrapper
                   where the scintillation effect stopped working when not hit at exactly the middle of the scintillator
           */
            // Create the housing logical volume
            auto housing_log = make_shared_no_delete<G4LogicalVolume>(
                housing_solid_.get(), materials_[housing_material], "housing_" + name + "_log");
            detector->setExternalObject("housing_log", housing_log);
            housing_phys_ = make_shared_no_delete<G4PVPlacement>(
                transform_phys, housing_log.get(), "housing_" + name + "_phys", world_log_.get(), false, 0, true);
            detector->setExternalObject("housing_phys", housing_phys_);

            /* Scintillator
            * the scintillator is the part that creates the optical photons
            */
            auto scint_log = make_shared_no_delete<G4LogicalVolume>(
                scint_solid_.get(), materials_[scint_material], "scint_" + name + "_log");
            detector->setExternalObject("scint_log", scint_log);

            // Place the scintillator box inside the housing
            ROOT::Math::XYZVector scint_displacement = {0, 0, sensor_size.z() / 2.0};
            auto scint_pos = toG4Vector(scint_displacement);
            LOG(DEBUG) << "  - Scintillator\t\t:\t" << Units::display(scint_pos, {"mm", "um"});
            scint_phys_ = make_shared_no_delete<G4PVPlacement>(
                nullptr, scint_pos, scint_log.get(), "scint_" + name + "_phys", housing_log.get(), false, 0, true);
            detector->setExternalObject("scint_phys", scint_phys_);

            /* SENSOR
            * the sensitive photocathode is the part that collects the optical photons
            */
            // Create the sensor box and logical volume
            auto sensor_box = std::make_shared<G4Box>(
                "sensor_" + name, sensor_size.x() / 2.0, sensor_size.y() / 2.0, sensor_size.z() / 2.0);
            solids_.push_back(sensor_box);
            auto sensor_log = make_shared_no_delete<G4LogicalVolume>(
                sensor_box.get(), materials_[sensor_material], "sensor_" + name + "_log");
            detector->setExternalObject("sensor_log", sensor_log);

            // Place the sensor box
            ROOT::Math::XYZVector sensor_displacement = {0, 0, -scint_size.z() / 2.0};
            auto sensor_pos = toG4Vector(sensor_displacement);
            LOG(DEBUG) << "  - Sensor\t\t:\t" << Units::display(sensor_pos, {"mm", "um"});
            sensor_phys_ = make_shared_no_delete<G4PVPlacement>(
                nullptr, sensor_pos, sensor_log.get(), "sensor_" + name + "_phys", housing_log.get(), false, 0, true);
            detector->setExternalObject("sensor_phys", sensor_phys_);

            // FIXME:: Find somewhere to store this information for multiple scintillators

            // General Info
            G4double cebr3_Energy[] = {3.96 * CLHEP::eV,
                                       3.81 * CLHEP::eV,
                                       3.67 * CLHEP::eV,
                                       3.54 * CLHEP::eV,
                                       3.42 * CLHEP::eV,
                                       3.31 * CLHEP::eV,
                                       3.20 * CLHEP::eV,
                                       3.10 * CLHEP::eV,
                                       3.0 * CLHEP::eV,
                                       2.92 * CLHEP::eV,
                                       2.83 * CLHEP::eV,
                                       2.76 * CLHEP::eV,
                                       2.68 * CLHEP::eV,
                                       1.77 * CLHEP::eV};
            const G4int num = sizeof(cebr3_Energy) / sizeof(G4double);

            // Housing Properties
            G4double reflectivity[] = {housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity,
                                       housing_reflectivity};
            assert(sizeof(reflectivity) == sizeof(cebr3_Energy));
            G4double efficiency[] = {
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            };
            assert(sizeof(efficiency) == sizeof(cebr3_Energy));
            auto scintHsngPT = new G4MaterialPropertiesTable();
            scintHsngPT->AddProperty("REFLECTIVITY", cebr3_Energy, reflectivity, num);
            scintHsngPT->AddProperty("EFFICIENCY", cebr3_Energy, efficiency, num);

            auto OpScintHousingSurface = new G4OpticalSurface("HousingSurface",
                                                              housing_surface_model,
                                                              housing_surface_finish,
                                                              housing_surface_type,
                                                              housing_surface_value);
            OpScintHousingSurface->SetMaterialPropertiesTable(scintHsngPT);

            // Photocathode Properties
            G4double photocath_EFF[] = {
                1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.}; // Enables 'detection' of photons
            assert(sizeof(photocath_EFF) == sizeof(cebr3_Energy));
            G4double photocath_ReR[] = {1.92, 1.92, 1.92, 1.92, 1.92, 1.92, 1.92, 1.92, 1.92, 1.92, 1.92, 1.92, 1.92, 1.92};
            assert(sizeof(photocath_ReR) == sizeof(cebr3_Energy));
            G4double photocath_ImR[] = {1.69, 1.69, 1.69, 1.69, 1.69, 1.69, 1.69, 1.69, 1.69, 1.69, 1.69, 1.69, 1.69, 1.69};
            assert(sizeof(photocath_ImR) == sizeof(cebr3_Energy));
            auto photocath_mt = new G4MaterialPropertiesTable();
            photocath_mt->AddProperty("EFFICIENCY", cebr3_Energy, photocath_EFF, num);
            photocath_mt->AddProperty("REALRINDEX", cebr3_Energy, photocath_ReR, num);
            photocath_mt->AddProperty("IMAGINARYRINDEX", cebr3_Energy, photocath_ImR, num);
            auto photocath_opsurf = new G4OpticalSurface("photocath_opsurf",
                                                         photocathode_surface_model,
                                                         photocathode_surface_finish,
                                                         photocathode_surface_type,
                                                         photocathode_surface_value);
            photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

            //  Create logical skin surfaces
            new G4LogicalBorderSurface("photocath_surf", scint_phys_.get(), sensor_phys_.get(), photocath_opsurf);
            new G4LogicalBorderSurface("housing_surf", scint_phys_.get(), housing_phys_.get(), OpScintHousingSurface);

        } else {
            /*  Detector
                    Creates a detector with a sensitive sensor with pixels, a chip and an optional support layer
                */

            LOG(TRACE) << " Pixel Detector Model Parameters for detector " << name;
            LOG(TRACE) << " Chip dimensions: " << model->getChipSize();
            LOG(DEBUG) << " Center of the geometry parts relative to the detector wrapper geometric center:";

            // Create the wrapper box and logical volume
            auto wrapper_box = std::make_shared<G4Box>(
                "wrapper_" + name, model->getSize().x() / 2.0, model->getSize().y() / 2.0, model->getSize().z() / 2.0);
            solids_.push_back(wrapper_box);
            auto wrapper_log =
                make_shared_no_delete<G4LogicalVolume>(wrapper_box.get(), world_material_, "wrapper_" + name + "_log");
            detector->setExternalObject("wrapper_log", wrapper_log);

            // Place the wrapper
            auto wrapper_phys = make_shared_no_delete<G4PVPlacement>(
                transform_phys, wrapper_log.get(), "wrapper_" + name + "_phys", world_log_.get(), false, 0, true);
            detector->setExternalObject("wrapper_phys", wrapper_phys);

            /* SENSOR
            * the sensitive detector is the part that collects the deposits
            */

            // Create the sensor box and logical volume
            auto sensor_box = std::make_shared<G4Box>("sensor_" + name,
                                                      model->getSensorSize().x() / 2.0,
                                                      model->getSensorSize().y() / 2.0,
                                                      model->getSensorSize().z() / 2.0);
            solids_.push_back(sensor_box);
            auto sensor_log = make_shared_no_delete<G4LogicalVolume>(
                sensor_box.get(), materials_[sensor_material], "sensor_" + name + "_log");
            detector->setExternalObject("sensor_log", sensor_log);

            // Place the sensor box
            auto sensor_pos = toG4Vector(model->getSensorCenter() - model->getGeometricalCenter());
            LOG(DEBUG) << "  - Sensor\t\t:\t" << Units::display(sensor_pos, {"mm", "um"});
            sensor_phys_ = make_shared_no_delete<G4PVPlacement>(
                nullptr, sensor_pos, sensor_log.get(), "sensor_" + name + "_phys", wrapper_log.get(), false, 0, true);
            detector->setExternalObject("sensor_phys", sensor_phys_);

            // Create the pixel box and logical volume
            auto pixel_box = std::make_shared<G4Box>("pixel_" + name,
                                                     model->getPixelSize().x() / 2.0,
                                                     model->getPixelSize().y() / 2.0,
                                                     model->getSensorSize().z() / 2.0);
            solids_.push_back(pixel_box);
            auto pixel_log =
                make_shared_no_delete<G4LogicalVolume>(pixel_box.get(), materials_["silicon"], "pixel_" + name + "_log");
            detector->setExternalObject("pixel_log", pixel_log);

            // Create the parameterization for the pixel grid
            std::shared_ptr<G4VPVParameterisation> pixel_param =
                std::make_shared<Parameterization2DG4>(model->getNPixels().x(),
                                                       model->getPixelSize().x(),
                                                       model->getPixelSize().y(),
                                                       -model->getGridSize().x() / 2.0,
                                                       -model->getGridSize().y() / 2.0,
                                                       0);
            detector->setExternalObject("pixel_param", pixel_param);

            // WARNING: do not place the actual parameterization, only use it if we need it

            /* CHIP
            * the chip connected to the bumps bond and the support
            */

            // Construct the chips only if necessary
            if(model->getChipSize().z() > 1e-9) {
                // Create the chip box
                auto chip_box = std::make_shared<G4Box>("chip_" + name,
                                                        model->getChipSize().x() / 2.0,
                                                        model->getChipSize().y() / 2.0,
                                                        model->getChipSize().z() / 2.0);
                solids_.push_back(chip_box);

                // Create the logical volume for the chip
                auto chip_log =
                    make_shared_no_delete<G4LogicalVolume>(chip_box.get(), materials_["silicon"], "chip_" + name + "_log");
                detector->setExternalObject("chip_log", chip_log);

                // Place the chip
                auto chip_pos = toG4Vector(model->getChipCenter() - model->getGeometricalCenter());
                LOG(DEBUG) << "  - Chip\t\t:\t" << Units::display(chip_pos, {"mm", "um"});
                auto chip_phys = make_shared_no_delete<G4PVPlacement>(
                    nullptr, chip_pos, chip_log.get(), "chip_" + name + "_phys", wrapper_log.get(), false, 0, true);
                detector->setExternalObject("chip_phys", chip_phys);
            }

            /*
            * SUPPORT
            * optional layers of support
            */
            auto supports_log = std::make_shared<std::vector<std::shared_ptr<G4LogicalVolume>>>();
            auto supports_phys = std::make_shared<std::vector<std::shared_ptr<G4PVPlacement>>>();
            int support_idx = 0;
            for(auto& layer : model->getSupportLayers()) {
                // Create the box containing the support
                auto support_box = std::make_shared<G4Box>("support_" + name + "_" + std::to_string(support_idx),
                                                           layer.getSize().x() / 2.0,
                                                           layer.getSize().y() / 2.0,
                                                           layer.getSize().z() / 2.0);
                solids_.push_back(support_box);

                std::shared_ptr<G4VSolid> support_solid = support_box;
                if(layer.hasHole()) {
                    // NOTE: Double the hole size in the z-direction to ensure no fake surfaces are created
                    auto hole_box = std::make_shared<G4Box>("support_" + name + "_hole_" + std::to_string(support_idx),
                                                            layer.getHoleSize().x() / 2.0,
                                                            layer.getHoleSize().y() / 2.0,
                                                            layer.getHoleSize().z());
                    solids_.push_back(hole_box);

                    G4Transform3D transform(G4RotationMatrix(), toG4Vector(layer.getHoleCenter() - layer.getCenter()));
                    auto subtraction_solid = std::make_shared<G4SubtractionSolid>("support_" + name + "_subtraction_" +
                                                                                      std::to_string(support_idx),
                                                                                  support_box.get(),
                                                                                  hole_box.get(),
                                                                                  transform);
                    solids_.push_back(subtraction_solid);
                    support_solid = subtraction_solid;
                }

                // Create the logical volume for the support
                auto support_material_iter = materials_.find(layer.getMaterial());
                if(support_material_iter == materials_.end()) {
                    throw ModuleError("Cannot construct a support layer of material '" + layer.getMaterial() + "'");
                }
                auto support_log =
                    make_shared_no_delete<G4LogicalVolume>(support_solid.get(),
                                                           support_material_iter->second,
                                                           "support_" + name + "_log_" + std::to_string(support_idx));
                supports_log->push_back(support_log);

                // Place the support
                auto support_pos = toG4Vector(layer.getCenter() - model->getGeometricalCenter());
                LOG(DEBUG) << "  - Support\t\t:\t" << Units::display(support_pos, {"mm", "um"});
                auto support_phys =
                    make_shared_no_delete<G4PVPlacement>(nullptr,
                                                         support_pos,
                                                         support_log.get(),
                                                         "support_" + name + "_phys_" + std::to_string(support_idx),
                                                         wrapper_log.get(),
                                                         false,
                                                         0,
                                                         true);
                supports_phys->push_back(support_phys);

                ++support_idx;
            }
            detector->setExternalObject("supports_log", supports_log);
            detector->setExternalObject("supports_phys", supports_phys);

            // Build the bump bonds only for hybrid pixel detectors
            auto hybrid_model = std::dynamic_pointer_cast<HybridPixelDetectorModel>(model);
            if(hybrid_model != nullptr) {
                /* BUMPS
                * the bump bonds connect the sensor to the readout chip
                */

                // Get parameters from model
                auto bump_height = hybrid_model->getBumpHeight();
                auto bump_sphere_radius = hybrid_model->getBumpSphereRadius();
                auto bump_cylinder_radius = hybrid_model->getBumpCylinderRadius();

                // Create the volume containing the bumps
                auto bump_box = std::make_shared<G4Box>("bump_box_" + name,
                                                        hybrid_model->getSensorSize().x() / 2.0,
                                                        hybrid_model->getSensorSize().y() / 2.0,
                                                        bump_height / 2.);
                solids_.push_back(bump_box);

                // Create the logical wrapper volume
                auto bumps_wrapper_log = make_shared_no_delete<G4LogicalVolume>(
                    bump_box.get(), world_material_, "bumps_wrapper_" + name + "_log");
                detector->setExternalObject("bumps_wrapper_log", bumps_wrapper_log);

                // Place the general bumps volume
                G4ThreeVector bumps_pos = toG4Vector(hybrid_model->getBumpsCenter() - hybrid_model->getGeometricalCenter());
                LOG(DEBUG) << "  - Bumps\t\t:\t" << Units::display(bumps_pos, {"mm", "um"});
                auto bumps_wrapper_phys = make_shared_no_delete<G4PVPlacement>(nullptr,
                                                                               bumps_pos,
                                                                               bumps_wrapper_log.get(),
                                                                               "bumps_wrapper_" + name + "_phys",
                                                                               wrapper_log.get(),
                                                                               false,
                                                                               0,
                                                                               true);
                detector->setExternalObject("bumps_wrapper_phys", bumps_wrapper_phys);

                // Create the individual bump solid
                auto bump_sphere = std::make_shared<G4Sphere>(
                    "bumps_" + name + "_sphere", 0, bump_sphere_radius, 0, 360 * CLHEP::deg, 0, 360 * CLHEP::deg);
                solids_.push_back(bump_sphere);
                auto bump_tube = std::make_shared<G4Tubs>(
                    "bumps_" + name + "_tube", 0., bump_cylinder_radius, bump_height / 2., 0., 360 * CLHEP::deg);
                solids_.push_back(bump_tube);
                auto bump = std::make_shared<G4UnionSolid>("bumps_" + name, bump_sphere.get(), bump_tube.get());
                solids_.push_back(bump);

                // Create the logical volume for the individual bumps
                auto bumps_cell_log =
                    make_shared_no_delete<G4LogicalVolume>(bump.get(), materials_["solder"], "bumps_" + name + "_log");
                detector->setExternalObject("bumps_cell_log", bumps_cell_log);

                // Place the bump bonds grid
                std::shared_ptr<G4VPVParameterisation> bumps_param = std::make_shared<Parameterization2DG4>(
                    hybrid_model->getNPixels().x(),
                    hybrid_model->getPixelSize().x(),
                    hybrid_model->getPixelSize().y(),
                    -(hybrid_model->getNPixels().x() * hybrid_model->getPixelSize().x()) / 2.0 +
                        (hybrid_model->getBumpsCenter().x() - hybrid_model->getCenter().x()),
                    -(hybrid_model->getNPixels().y() * hybrid_model->getPixelSize().y()) / 2.0 +
                        (hybrid_model->getBumpsCenter().y() - hybrid_model->getCenter().y()),
                    0);
                detector->setExternalObject("bumps_param", bumps_param);

                std::shared_ptr<G4PVParameterised> bumps_param_phys =
                    std::make_shared<ParameterisedG4>("bumps_" + name + "_phys",
                                                      bumps_cell_log.get(),
                                                      bumps_wrapper_log.get(),
                                                      kUndefined,
                                                      hybrid_model->getNPixels().x() * hybrid_model->getNPixels().y(),
                                                      bumps_param.get(),
                                                      false);
                detector->setExternalObject("bumps_param_phys", bumps_param_phys);
            }

            // ALERT: NO COVER LAYER YET

            LOG(TRACE) << " Constructed detector " << detector->getName() << " successfully";
        }
    }
}
