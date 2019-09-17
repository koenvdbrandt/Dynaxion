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
#include "G4OpticalSurface.hh"
#include "G4SystemOfUnits.hh"

#include "core/geometry/HybridPixelDetectorModel.hpp"
#include "core/geometry/ScintillatorModel.hpp"
#include "core/module/exceptions.h"
#include "core/utils/log.h"
#include "tools/ROOT.h"
#include "tools/geant4.h"

#include "GeometryConstructionG4.hpp"
#include "Parameterization2DG4.hpp"

using namespace allpix;

DetectorConstructionG4::DetectorConstructionG4(GeometryManager* geo_manager, Configuration& config)
    : geo_manager_(geo_manager), config_(config) {}

/**
 * @brief Version of std::make_shared that does not delete the pointer
 *
 * This version is needed because some pointers are deleted by Geant4 internally, but they are stored as std::shared_ptr in
 * the framework.
 */
template <typename T, typename... Args> static std::shared_ptr<T> make_shared_no_delete(Args... args) {
    return std::shared_ptr<T>(new T(args...), [](T*) {});
}

void DetectorConstructionG4::build(G4LogicalVolume* world_log, std::map<std::string, G4Material*> materials_) {

    /*
    Get the world material
    */
    std::string world_material = config_.get<std::string>("world_material", "air");
    world_material_ = materials_[world_material];
    /*
    Build the individual detectors
    */
    std::vector<std::shared_ptr<Detector>> detectors = geo_manager_->getDetectors();
    LOG(TRACE) << "Building " << detectors.size() << " device(s)";

    for(auto& detector : detectors) {
        // Get pointer to the model of the detector
        auto model = detector->getModel();

        std::string name = detector->getName();
        LOG(DEBUG) << "Creating Geant4 model for " << name;
        LOG(DEBUG) << " Wrapper dimensions of model: " << Units::display(model->getSize(), {"mm", "um"});
        LOG(TRACE) << " Sensor dimensions: " << model->getSensorSize();
        LOG(TRACE) << " Chip dimensions: " << model->getChipSize();
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

        LOG(DEBUG) << " Center of the geometry parts relative to the detector wrapper geometric center:";

        std::map<std::string, std::string> type = geo_manager_->getDetectorType();
        if(type[detector->getType()] == "scintillator") {
            auto scint_model = std::dynamic_pointer_cast<ScintillatorModel>(model);

            // Get parameters from model
            auto sensor_size = scint_model->getSensorSize();
            auto sensor_material = scint_model->getActiveMaterial();
            auto scint_size = scint_model->getScintSize();
            auto scint_material = scint_model->getScintMaterial();
            auto housing_thickness = scint_model->getHousingThickness();
            auto housing_material = scint_model->getHousingMaterial();
            housing_reflectivity = scint_model->getHousingReflectivity();
            /*   Housing
                    the housing of the scintillator
                */

            // Create the volume containing the housing
            auto housing_box = std::make_shared<G4Box>("housing_box_" + name,
                                                       scint_size.x() / 2.0 + 2 * housing_thickness,
                                                       scint_size.y() / 2.0 + sensor_size.y() / 2.0 + 2 * housing_thickness,
                                                       scint_size.z() / 2.0 + 2 * housing_thickness);
            solids_.push_back(housing_box);

            // Create the housing logical volume
            housing_log = make_shared_no_delete<G4LogicalVolume>(
                housing_box.get(), materials_[housing_material], "housing_" + name + "_log");
            detector->setExternalObject("housing_log", housing_log);
            ROOT::Math::XYZVector housing_displacement = {0, scint_size.y() / 2.0, 0};
            auto housing_pos = posWrapper + toG4Vector(housing_displacement);
            G4Transform3D transform_phys_housing(*rotWrapper, housing_pos);
            auto housing_phys = make_shared_no_delete<G4PVPlacement>(
                transform_phys_housing, housing_log.get(), "housing_" + name + "_phys", world_log, false, 0, true);
            detector->setExternalObject("housing_phys", housing_phys);

            /* Scintillator
             * the scintillator is the part that creates the optical photons
            */
            // Create the scintillator box and logical volume
            auto scint_box =
                std::make_shared<G4Box>("scint_" + name, scint_size.x() / 2.0, scint_size.y() / 2.0, scint_size.z() / 2.0);
            solids_.push_back(scint_box);
            auto scint_log = make_shared_no_delete<G4LogicalVolume>(
                scint_box.get(), materials_[scint_material], "scint_" + name + "_log");
            detector->setExternalObject("scint_log", scint_log);

            // Place the scintillator box
            ROOT::Math::XYZVector scint_displacement = {0, sensor_size.y() / 2.0, 0};
            auto scint_pos = toG4Vector(scint_displacement);
            LOG(DEBUG) << "  - Scintillator\t\t:\t" << Units::display(scint_pos, {"mm", "um"});
            auto scint_phys = make_shared_no_delete<G4PVPlacement>(
                nullptr, scint_pos, scint_log.get(), "scint_" + name + "_phys", housing_log.get(), false, 0, true);
            detector->setExternalObject("scint_phys", scint_phys);

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
            ROOT::Math::XYZVector sensor_displacement = {0, -scint_size.y() / 2.0, 0};
            LOG(WARNING) << "  - Sensor\t\t:\t" << Units::display(sensor_displacement, {"mm", "um"});
            LOG(WARNING) << "  - Scint_size_y" << scint_size.y();
            LOG(WARNING) << "  - Sensor_size_y" << sensor_size.y();

            auto sensor_pos = toG4Vector(sensor_displacement);
            LOG(DEBUG) << "  - Sensor\t\t:\t" << Units::display(sensor_pos, {"mm", "um"});
            auto sensor_phys = make_shared_no_delete<G4PVPlacement>(
                nullptr, sensor_pos, sensor_log.get(), "sensor_" + name + "_phys", housing_log.get(), false, 0, true);
            detector->setExternalObject("sensor_phys", sensor_phys);

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
                                       2.61 * CLHEP::eV};
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
            G4OpticalSurface* OpScintHousingSurface =
                new G4OpticalSurface("HousingSurface", unified, polished, dielectric_metal);
            OpScintHousingSurface->SetMaterialPropertiesTable(scintHsngPT);

            // Scintillator Properties
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
            G4OpticalSurface* photocath_opsurf =
                new G4OpticalSurface("photocath_opsurf", glisur, polished, dielectric_metal);
            photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

            //**Create logical skin surfaces
            new G4LogicalSkinSurface("photocath_surf", sensor_log.get(), photocath_opsurf);
            new G4LogicalSkinSurface("housing_surf", housing_log.get(), OpScintHousingSurface);

            /*   PMT
                    the PMT
                                  
            if(scint_model->getPMType() == "PMT"){
                // Get parameters from model
                auto PMT_outer_radius = scint_model->getPMTOuterRadius();
                auto PMT_height = scint_model->getPMTHeight();
                auto PMT_material = scint_model->getPMTMaterial();
                auto photo_cath_material = scint_model->getPhotoCathMaterial();
                // Create the volume containing the housing
                auto PMT_tube = std::make_shared<G4Tubs>("pmt_tube",0,PMT_outer_radius,
                        PMT_height/2,0,2*CLHEP::pi);

                solids_.push_back(PMT_tube);
                auto photo_cath_tube = std::make_shared<G4Tubs>("photo_cath_tube",0,PMT_outer_radius,
                        PMT_height/4,0,2*CLHEP::pi);
                solids_.push_back(photo_cath_tube);

                // Create the PMT logical volume
                PMT_log =
                    make_shared_no_delete<G4LogicalVolume>(PMT_tube.get(), materials_[PMT_material], "PMT_" + name + "_log");
                detector->setExternalObject("PMT_log", PMT_log);
                photo_cath_log =
                    make_shared_no_delete<G4LogicalVolume>(photo_cath_tube.get(), materials_[photo_cath_material], "photo_cath_" + name + "_log");
                detector->setExternalObject("photo_cath_log", photo_cath_log);
                new G4PVPlacement(0,G4ThreeVector(0,0,-PMT_height/4),
                                        photo_cath_log.get(),"photocath",
                                        PMT_log.get(),false,0);
                detector->setExternalObject("photo_cath_log", photo_cath_log);
        
                //Place the PMT's around the Volume
                auto Nx =scint_model->getNx();
                auto Ny =scint_model->getNy();
                auto Nz =scint_model->getNz();
                G4double dx = scint_size.x()/Nx;
                G4double dy = scint_size.y()/Ny;
                G4double dz = scint_size.z()/Nz;
                
                G4double x,y,z;
                G4double xmin = -scint_size.x()/2. - dx/2.;
                G4double ymin = -scint_size.y()/2. - dy/2.;
                G4double zmin = -scint_size.z()/2. - dz/2.;
                G4int k=0;
                
                z = -scint_size.z()/2. - PMT_height/2;      //front
                PlacePMTs(PMT_log.get(),0,x,y,dx,dy,xmin,ymin,Nx,Ny,x,y,z,k);

                G4RotationMatrix* rm_z = new G4RotationMatrix();
                rm_z->rotateY(180*CLHEP::deg);
                z = scint_size.z()/2. + PMT_height/2;       //back
                PlacePMTs(PMT_log.get(),rm_z,x,y,dx,dy,xmin,ymin,Nx,Ny,x,y,z,k);
                
                G4RotationMatrix* rm_y1 = new G4RotationMatrix();
                rm_y1->rotateY(-90*CLHEP::deg);
                x = -scint_size.x()/2. - PMT_height/2;      //left
                PlacePMTs(PMT_log.get(),rm_y1,y,z,dy,dz,ymin,zmin,Ny,Nz,x,y,z,k);

                G4RotationMatrix* rm_y2 = new G4RotationMatrix();
                rm_y2->rotateY(90*CLHEP::deg);
                x = scint_size.x()/2. + PMT_height/2;      //right
                PlacePMTs(PMT_log.get(),rm_y2,y,z,dy,dz,ymin,zmin,Ny,Nz,x,y,z,k);
                
                G4RotationMatrix* rm_x1 = new G4RotationMatrix();
                rm_x1->rotateX(90*CLHEP::deg);
                y = -scint_size.y()/2. - PMT_height/2;     //bottom
                PlacePMTs(PMT_log.get(),rm_x1,x,z,dx,dz,xmin,zmin,Nx,Nz,x,y,z,k);

                G4RotationMatrix* rm_x2 = new G4RotationMatrix();
                rm_x2->rotateX(-90*CLHEP::deg);
                y = scint_size.y()/2. + PMT_height/2;      //top
                PlacePMTs(PMT_log.get(),rm_x2,x,z,dx,dz,xmin,zmin,Nx,Nz,x,y,z,k);
                
                //VisAttributes();
                //SurfaceProperties();

                //SetLogicalVolume(housing_log.get());
            }
             */ /*
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
                                                         housing_log.get(),
                                                         false,
                                                         0,
                                                         true);
                supports_phys->push_back(support_phys);

                ++support_idx;
            }
            detector->setExternalObject("supports_log", supports_log);
            detector->setExternalObject("supports_phys", supports_phys);
        } else {
            // Create the wrapper box and logical volume
            auto wrapper_box = std::make_shared<G4Box>(
                "wrapper_" + name, model->getSize().x() / 2.0, model->getSize().y() / 2.0, model->getSize().z() / 2.0);
            solids_.push_back(wrapper_box);
            auto wrapper_log =
                make_shared_no_delete<G4LogicalVolume>(wrapper_box.get(), world_material_, "wrapper_" + name + "_log");
            detector->setExternalObject("wrapper_log", wrapper_log);
            // Place the wrapper
            auto wrapper_phys = make_shared_no_delete<G4PVPlacement>(
                transform_phys, wrapper_log.get(), "wrapper_" + name + "_phys", world_log, false, 0, true);
            detector->setExternalObject("wrapper_phys", wrapper_phys);

            /*
                SENSOR
                * the sensitive detector is the part that collects the deposits
            */
            // Create the sensor box and logical volume
            ROOT::Math::XYZVector sensor_size = model->getSensorSize();
            std::string sensor_material = model->getActiveMaterial();

            auto sensor_box = std::make_shared<G4Box>(
                "sensor_" + name, sensor_size.x() / 2.0, sensor_size.y() / 2.0, sensor_size.z() / 2.0);
            solids_.push_back(sensor_box);
            auto sensor_log = make_shared_no_delete<G4LogicalVolume>(
                sensor_box.get(), materials_[sensor_material], "sensor_" + name + "_log");
            detector->setExternalObject("sensor_log", sensor_log);

            // Place the sensor box
            auto sensor_pos = toG4Vector(model->getSensorCenter() - model->getGeometricalCenter());
            LOG(DEBUG) << "  - Sensor\t\t:\t" << Units::display(sensor_pos, {"mm", "um"});
            auto sensor_phys = make_shared_no_delete<G4PVPlacement>(
                nullptr, sensor_pos, sensor_log.get(), "sensor_" + name + "_phys", wrapper_log.get(), false, 0, true);
            detector->setExternalObject("sensor_phys", sensor_phys);

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
            /*
                    CHIP
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

                /*   BUMPS
                        the bump bonds connect the sensor to the readout chip
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
        }
        // ALERT: NO COVER LAYER YET

        LOG(TRACE) << " Constructed detector " << detector->getName() << " successfully";
    }
}

void DetectorConstructionG4::PlacePMTs(G4LogicalVolume* pmt_log,
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
                                       G4int& k) {
    /*PlacePMTs : a different way to parameterize placement that does not depend on
    calculating the position from the copy number

    pmt_log = logical volume for pmts to be placed
    rot = rotation matrix to apply
    a,b = coordinates to vary(ie. if varying in the xy plane then pass x,y)
    da,db = value to increment a,b by
    amin,bmin = start values for a,b
    na,nb = number of repitions in a and b
    x,y,z = just pass x,y, and z by reference (the same ones passed for a,b)
    k = copy number to start with
    sd = sensitive detector for pmts
    */
    a = amin;
    for(G4int j = 1; j <= na; j++) {
        a += da;
        b = bmin;
        for(G4int i = 1; i <= nb; i++) {
            b += db;
            new G4PVPlacement(rot, G4ThreeVector(x, y, z), pmt_log, "pmt", housing_log.get(), false, k);
            PMT_positions.emplace_back(G4ThreeVector(x, y, z));
            k++;
        }
    }
}

/* void DetectorConstructionG4::SurfaceProperties(){
    G4double ephoton[] = {7.0*CLHEP::eV, 7.14*CLHEP::eV};
    const G4int num = sizeof(ephoton)/sizeof(G4double);

    // **Scintillator housing properties
    G4double reflectivity[] = {housing_reflectivity, housing_reflectivity};
    assert(sizeof(reflectivity) == sizeof(ephoton));
    G4double efficiency[] = {0.0, 0.0};
    assert(sizeof(efficiency) == sizeof(ephoton));
    G4MaterialPropertiesTable* scintHsngPT = new G4MaterialPropertiesTable();
    scintHsngPT->AddProperty("REFLECTIVITY", ephoton, reflectivity, num);
    scintHsngPT->AddProperty("EFFICIENCY", ephoton, efficiency, num);
    G4OpticalSurface* scint_housing_opsurf =
        new G4OpticalSurface("housing_opsurf",unified,polished,dielectric_metal);
    scint_housing_opsurf->SetMaterialPropertiesTable(scintHsngPT);

    // **Photocathode surface properties
    G4double photocath_EFF[]={1.,1.}; //Enables 'detection' of photons
    assert(sizeof(photocath_EFF) == sizeof(ephoton));
    G4double photocath_ReR[]={1.92,1.92};
    assert(sizeof(photocath_ReR) == sizeof(ephoton));
    G4double photocath_ImR[]={1.69,1.69};
    assert(sizeof(photocath_ImR) == sizeof(ephoton));
    G4MaterialPropertiesTable* photocath_mt = new G4MaterialPropertiesTable();
    photocath_mt->AddProperty("EFFICIENCY",ephoton,photocath_EFF,num);
    photocath_mt->AddProperty("REALRINDEX",ephoton,photocath_ReR,num);
    photocath_mt->AddProperty("IMAGINARYRINDEX",ephoton,photocath_ImR,num);
    G4OpticalSurface* photocath_opsurf=
        new G4OpticalSurface("photocath_opsurf",glisur,polished,
                            dielectric_metal);
    photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

    // **Create logical skin surfaces
    new G4LogicalSkinSurface("housign_surf", housing_log.get(),
                            scint_housing_opsurf);
    new G4LogicalSkinSurface("photocath_surf", photo_cath_log.get(), photocath_opsurf);
}*/
