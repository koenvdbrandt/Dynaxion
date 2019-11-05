/**
 * @file
 * @brief Implements the Geant4 geometry construction process
 * @remarks Code is based on code from Mathieu Benoit
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "GeometryConstructionG4.hpp"

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

#include "Parameterization2DG4.hpp"

using namespace allpix;

GeometryConstructionG4::GeometryConstructionG4(GeometryManager* geo_manager, Configuration& config)
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

/**
 * First initializes all the materials. Then constructs the world from the internally calculated world size with a certain
 * margin. Finally builds all the individual detectors.
 */
G4VPhysicalVolume* GeometryConstructionG4::Construct() {
    // Initialize materials
    init_materials();

    // Set world material
    std::string world_material = config_.get<std::string>("world_material", "air");
    std::transform(world_material.begin(), world_material.end(), world_material.begin(), ::tolower);
    if(materials_.find(world_material) == materials_.end()) {
        throw InvalidValueError(config_, "world_material", "material does not exists, use 'air' or 'vacuum'");
    }

    world_material_ = materials_[world_material];
    LOG(TRACE) << "Material of world is " << world_material_->GetName();

    // Calculate world size
    ROOT::Math::XYZVector half_world_size;
    ROOT::Math::XYZPoint min_coord = geo_manager_->getMinimumCoordinate();
    ROOT::Math::XYZPoint max_coord = geo_manager_->getMaximumCoordinate();
    half_world_size.SetX(std::max(std::abs(min_coord.x()), std::abs(max_coord.x())));
    half_world_size.SetY(std::max(std::abs(min_coord.y()), std::abs(max_coord.y())));
    half_world_size.SetZ(std::max(std::abs(min_coord.z()), std::abs(max_coord.z())));

    // Calculate and apply margins to world size
    auto margin_percentage = config_.get<double>("world_margin_percentage", 0.1);
    auto minimum_margin = config_.get<ROOT::Math::XYZPoint>("world_minimum_margin", {0, 0, 0});
    double add_x = half_world_size.x() * margin_percentage;
    if(add_x < minimum_margin.x()) {
        add_x = minimum_margin.x();
    }
    double add_y = half_world_size.y() * margin_percentage;
    if(add_y < minimum_margin.y()) {
        add_y = minimum_margin.y();
    }
    double add_z = half_world_size.z() * margin_percentage;
    if(add_z < minimum_margin.z()) {
        add_z = minimum_margin.z();
    }
    half_world_size.SetX(half_world_size.x() + add_x);
    half_world_size.SetY(half_world_size.y() + add_y);
    half_world_size.SetZ(half_world_size.z() + add_z);

    LOG(DEBUG) << "World size is " << Units::display(2 * half_world_size, {"mm"});

    // Build the world
    auto world_box = std::make_shared<G4Box>("World", half_world_size.x(), half_world_size.y(), half_world_size.z());
    solids_.push_back(world_box);
    world_log_ = std::make_unique<G4LogicalVolume>(world_box.get(), world_material_, "World", nullptr, nullptr, nullptr);

    // Set the world to invisible in the viewer
    world_log_->SetVisAttributes(G4VisAttributes::GetInvisible());

    // Place the world at the center
    world_phys_ =
        std::make_unique<G4PVPlacement>(nullptr, G4ThreeVector(0., 0., 0.), world_log_.get(), "World", nullptr, false, 0);

    // Build all the detectors in the world
    build_detectors();

    // Check for overlaps:
    check_overlaps();

    return world_phys_.get();
}

/**
 * Initializes all the internal materials. The following materials are supported by this module:
 * - vacuum
 * - air
 * - silicon
 * - epoxy
 * - kapton
 * - solder
 * - CeBr3 scintillator
 */
void GeometryConstructionG4::init_materials() {
    G4NistManager* nistman = G4NistManager::Instance();

    // Build table of materials from database
    materials_["silicon"] = nistman->FindOrBuildMaterial("G4_Si");
    materials_["plexiglass"] = nistman->FindOrBuildMaterial("G4_PLEXIGLASS");
    materials_["kapton"] = nistman->FindOrBuildMaterial("G4_KAPTON");
    materials_["copper"] = nistman->FindOrBuildMaterial("G4_Cu");
    materials_["aluminum"] = nistman->FindOrBuildMaterial("G4_Al");
    materials_["air"] = nistman->FindOrBuildMaterial("G4_AIR");
    materials_["lead"] = nistman->FindOrBuildMaterial("G4_Pb");
    materials_["tungsten"] = nistman->FindOrBuildMaterial("G4_W");

    // Create required elements:
    G4Element* H = new G4Element("Hydrogen", "H", 1., 1.01 * CLHEP::g / CLHEP::mole);
    G4Element* C = new G4Element("Carbon", "C", 6., 12.01 * CLHEP::g / CLHEP::mole);
    G4Element* O = new G4Element("Oxygen", "O", 8., 16.0 * CLHEP::g / CLHEP::mole);
    G4Element* Cl = new G4Element("Chlorine", "Cl", 17., 35.45 * CLHEP::g / CLHEP::mole);
    G4Element* Sn = new G4Element("Tin", "Sn", 50., 118.710 * CLHEP::g / CLHEP::mole);
    G4Element* Pb = new G4Element("Lead", "Pb", 82., 207.2 * CLHEP::g / CLHEP::mole);
    G4Element* Ce = new G4Element("Cerium", "Ce", 58., 140.12 * CLHEP::g / CLHEP::mole);
    G4Element* Br = new G4Element("Bromine", "Br", 35., 79.904 * CLHEP::g / CLHEP::mole);

    // Add vacuum
    materials_["vacuum"] = new G4Material("Vacuum", 1, 1.008 * CLHEP::g / CLHEP::mole, CLHEP::universe_mean_density);

    // Create Epoxy material
    G4Material* Epoxy = new G4Material("Epoxy", 1.3 * CLHEP::g / CLHEP::cm3, 3);
    Epoxy->AddElement(H, 44);
    Epoxy->AddElement(C, 15);
    Epoxy->AddElement(O, 7);
    materials_["epoxy"] = Epoxy;

    // Create Carbon Fiber material:
    G4Material* CarbonFiber = new G4Material("CarbonFiber", 1.5 * CLHEP::g / CLHEP::cm3, 2);
    CarbonFiber->AddMaterial(Epoxy, 0.4);
    CarbonFiber->AddElement(C, 0.6);
    materials_["carbonfiber"] = CarbonFiber;

    // Create PCB G-10 material
    G4Material* GTen = new G4Material("G10", 1.7 * CLHEP::g / CLHEP::cm3, 3);
    GTen->AddMaterial(nistman->FindOrBuildMaterial("G4_SILICON_DIOXIDE"), 0.773);
    GTen->AddMaterial(Epoxy, 0.147);
    GTen->AddElement(Cl, 0.08);
    materials_["g10"] = GTen;

    // Create solder material
    G4Material* Solder = new G4Material("Solder", 8.4 * CLHEP::g / CLHEP::cm3, 2);
    Solder->AddElement(Sn, 0.63);
    Solder->AddElement(Pb, 0.37);
    materials_["solder"] = Solder;

    // Create CeBr3
    G4Material* CeBr3 = new G4Material("CeBr3", 5.18 * CLHEP::g / CLHEP::cm3, 2);
    CeBr3->AddElement(Ce, 1);
    CeBr3->AddElement(Br, 3);
    materials_["cebr3"] = CeBr3;
    //***Material properties tables or CeBr3

    // Info CeBr3 from https://www.advatech-uk.co.uk/CeBr3%20Enission%20Spectrum.jpg
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
    const G4int cebr3num = sizeof(cebr3_Energy) / sizeof(G4double);
    G4double cebr3_SCINT[] = {0.01, 0.02, 0.04, 0.13, 0.36, 0.98, 0.98, 0.34, 0.20, 0.12, 0.08, 0.03, 0.02, 0.01};
    assert(sizeof(cebr3_SCINT) == sizeof(cebr3_Energy));

    // info from
    // https://www.researchgate.net/figure/The-calculated-optical-constants-of-CeCl-3-and-CeBr-3-a-refractive-index-b_fig7_256855384
    G4double cebr3_RIND[] = {2.40, 2.38, 2.35, 2.33, 2.30, 2.28, 2.25, 2.23, 2.20, 2.17, 2.10, 2.05, 2.10, 2.10};
    assert(sizeof(cebr3_RIND) == sizeof(cebr3_Energy));
    // Info from http://www.freepatentsonline.com/7405404.pdf
    G4double cebr3_ABSL[] = {2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm,
                             2.1 * CLHEP::cm};
    assert(sizeof(cebr3_ABSL) == sizeof(cebr3_Energy));
    auto CeBr3_mt = new G4MaterialPropertiesTable();
    CeBr3_mt->AddProperty("FASTCOMPONENT", cebr3_Energy, cebr3_SCINT, cebr3num);
    CeBr3_mt->AddProperty("RINDEX", cebr3_Energy, cebr3_RIND, cebr3num);
    CeBr3_mt->AddProperty("ABSLENGTH", cebr3_Energy, cebr3_ABSL, cebr3num);
    // Info https://www.ipen.br/biblioteca/cd/ieee/2004/DATA/3R04-1.PDF
    CeBr3_mt->AddConstProperty("SCINTILLATIONYIELD", 68000. / CLHEP::MeV);
    // ??
    CeBr3_mt->AddConstProperty("RESOLUTIONSCALE", 1.0);
    // Info  https://www.degruyter.com/downloadpdf/j/nuka.2017.62.issue-3/nuka-2017-0032/nuka-2017-0032.pdf
    CeBr3_mt->AddConstProperty("FASTTIMECONSTANT", 18. * CLHEP::ns);
    // https://www.berkeleynucleonics.com/sites/default/files/products/resources/cebr3_fast_timing_study_below_120_ps.pdf
    CeBr3_mt->AddConstProperty("FASTSCINTILLATIONRISETIME", 0.7 * CLHEP::ns);
    // Only fast component so yield = 1
    CeBr3_mt->AddConstProperty("YIELDRATIO", 1.0);
    CeBr3->SetMaterialPropertiesTable(CeBr3_mt);
    // ??
    // CeBr3->GetIonisation()->SetBirksConstant(0.126 * CLHEP::mm / CLHEP::MeV);
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
        LOG(TRACE) << " Sensor dimensions: " << model->getSensorSize();
        LOG(TRACE) << " Sensor material: " << sensor_material;
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
            auto housing_shape = scint_model->getHousingShape();
            auto housing_thickness = scint_model->getHousingThickness();
            auto housing_material = scint_model->getHousingMaterial();
            housing_reflectivity_ = scint_model->getHousingReflectivity();

            /*   Housing
                   the housing of the scintillator
                   housing works like the wrapper. There was some issue whith placing the housing in a wrapper
                   where the scintillation effect stopped working when not hit at exactly the middle of the scintillator
           */

            // Create the volume containing the housing
            if(housing_shape == "box") {
                housing_solid_ = std::make_shared<G4Box>("housing_" + name + "_solid",
                                                         scint_size.x() / 2.0 + housing_thickness,
                                                         scint_size.y() / 2.0 + housing_thickness,
                                                         scint_size.z() / 2.0 + housing_thickness + sensor_size.z() / 2.0);
            } else if(housing_shape == "cylinder") {
                housing_solid_ = std::make_shared<G4Tubs>("housing_" + name + "_solid",
                                                          0,
                                                          scint_size.x() / 2.0 + housing_thickness,
                                                          scint_size.z() / 2.0 + sensor_size.z() / 2.0 + housing_thickness,
                                                          0,
                                                          2 * CLHEP::pi);
            }
            solids_.push_back(housing_solid_);

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
            if(scint_shape == "box") {
                // Create the scintillator box and logical volume
                scint_solid_ = std::make_shared<G4Box>(
                    "scint_" + name, scint_size.x() / 2.0, scint_size.y() / 2.0, scint_size.z() / 2.0);
            } else if(scint_shape == "cylinder") {
                scint_solid_ = std::make_shared<G4Tubs>(
                    "scint_" + name, 0, scint_size.x() / 2.0, scint_size.z() / 2.0, 0, 2 * CLHEP::pi);
            }
            solids_.push_back(scint_solid_);
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
                                       2.61 * CLHEP::eV};
            const G4int num = sizeof(cebr3_Energy) / sizeof(G4double);

            // Housing Properties
            G4double reflectivity[] = {housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_,
                                       housing_reflectivity_};
            assert(sizeof(reflectivity) == sizeof(cebr3_Energy));
            G4double efficiency[] = {
                0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
            };
            assert(sizeof(efficiency) == sizeof(cebr3_Energy));
            auto scintHsngPT = new G4MaterialPropertiesTable();
            scintHsngPT->AddProperty("REFLECTIVITY", cebr3_Energy, reflectivity, num);
            scintHsngPT->AddProperty("EFFICIENCY", cebr3_Energy, efficiency, num);
            G4OpticalSurface* OpScintHousingSurface =
                new G4OpticalSurface("HousingSurface", unified, polishedteflonair, dielectric_metal);
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
            G4OpticalSurface* photocath_opsurf =
                new G4OpticalSurface("photocath_opsurf", glisur, polished, dielectric_metal);
            photocath_opsurf->SetMaterialPropertiesTable(photocath_mt);

            //  Create logical skin surfaces
            new G4LogicalBorderSurface("photocath_surf", scint_phys_.get(), sensor_phys_.get(), photocath_opsurf);
            new G4LogicalBorderSurface("housing_surf", scint_phys_.get(), housing_phys_.get(), OpScintHousingSurface);

        } else {
            /*  Detector
                    Creates a detector with a sensitive sensor with pixels, a chip and an optional support layer
                */
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

void GeometryConstructionG4::check_overlaps() {
    G4PhysicalVolumeStore* phys_volume_store = G4PhysicalVolumeStore::GetInstance();
    LOG(DEBUG) << phys_volume_store->size() << " physical volumes are defined";

    bool overlapFlag = false;
    for(auto volume : (*phys_volume_store)) {
        LOG(TRACE) << "Checking overlaps for physical volume \"" << volume->GetName() << "\"";
        overlapFlag = volume->CheckOverlaps(1000, 0., false) || overlapFlag;
    }
    if(overlapFlag) {
        LOG(ERROR) << "Overlapping volumes detected.";
    } else {
        LOG(INFO) << "No overlapping volumes detected.";
    }
}