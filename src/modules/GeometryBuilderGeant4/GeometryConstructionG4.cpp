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
#include "core/geometry/GeometryBuilder.hpp"
#include "core/geometry/HybridPixelDetectorModel.hpp"
#include "core/geometry/ScintillatorModel.hpp"
#include "core/module/exceptions.h"
#include "core/utils/log.h"
#include "tools/ROOT.h"
#include "tools/geant4.h"

#include "DetectorConstructionG4.hpp"
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
    // Add world_material to the list of materials to be called from other modules
    materials_["world_material"] = world_material_;
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

    // Build all the geometries that have been added to the GeometryBuilder vector, including Detectors and Passive Materials

    auto builders = geo_manager_->getBuilders();

    for(const auto& builder : builders) {
        auto g4_builder = std::dynamic_pointer_cast<GeometryBuilder<G4Material>>(builder);
        if(g4_builder != nullptr) {
            g4_builder->build(materials_);
        }
    }

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
 * - lithium
 * - beryllium
 * - CeBr3 scintillator
 */
void GeometryConstructionG4::init_materials() {
    G4NistManager* nistman = G4NistManager::Instance();

    // Build table of materials from database
    materials_["hydrogen"] = nistman->FindOrBuildMaterial("G4_H");
    materials_["carbon"] = nistman->FindOrBuildMaterial("G4_C");
    materials_["silicon"] = nistman->FindOrBuildMaterial("G4_Si");
    materials_["plexiglass"] = nistman->FindOrBuildMaterial("G4_PLEXIGLASS");
    materials_["kapton"] = nistman->FindOrBuildMaterial("G4_KAPTON");
    materials_["copper"] = nistman->FindOrBuildMaterial("G4_Cu");
    materials_["aluminum"] = nistman->FindOrBuildMaterial("G4_Al");
    materials_["air"] = nistman->FindOrBuildMaterial("G4_AIR");
    materials_["lead"] = nistman->FindOrBuildMaterial("G4_Pb");
    materials_["tungsten"] = nistman->FindOrBuildMaterial("G4_W");
    materials_["lithium"] = nistman->FindOrBuildMaterial("G4_Li");
    materials_["beryllium"] = nistman->FindOrBuildMaterial("G4_Be");
    materials_["al2o3"] = nistman->FindOrBuildMaterial("G4_ALUMINUM_OXIDE");
    materials_["teflon"] = nistman->FindOrBuildMaterial("G4_TEFLON");

    // Create required elements:
    G4Element* H = new G4Element("Hydrogen", "H", 1., 1.01 * CLHEP::g / CLHEP::mole);
    G4Element* C = new G4Element("Carbon", "C", 6., 12.01 * CLHEP::g / CLHEP::mole);
    G4Element* N = new G4Element("Nitrogen", "N", 7., 14.0 * CLHEP::g / CLHEP::mole);
    G4Element* O = new G4Element("Oxygen", "O", 8., 16.0 * CLHEP::g / CLHEP::mole);
    G4Element* Cl = new G4Element("Chlorine", "Cl", 17., 35.45 * CLHEP::g / CLHEP::mole);
    G4Element* Sn = new G4Element("Tin", "Sn", 50., 118.710 * CLHEP::g / CLHEP::mole);
    G4Element* Pb = new G4Element("Lead", "Pb", 82., 207.2 * CLHEP::g / CLHEP::mole);
    G4Element* Ce = new G4Element("Cerium", "Ce", 58., 140.12 * CLHEP::g / CLHEP::mole);
    G4Element* Br = new G4Element("Bromine", "Br", 35., 79.904 * CLHEP::g / CLHEP::mole);

    // Add vacuum
    materials_["vacuum"] = new G4Material("Vacuum", 1, 1.008 * CLHEP::g / CLHEP::mole, CLHEP::universe_mean_density);
   
    // Create g4_cocaine
    G4Material* g4_cocaine = new G4Material("g4_cocaine", 1.2 * CLHEP::g / CLHEP::cm3, 4);
    g4_cocaine->AddElement(C, 14);
    g4_cocaine->AddElement(H, 21);
    g4_cocaine->AddElement(N, 1);
    g4_cocaine->AddElement(O, 4);
    materials_["g4_cocaine"] = g4_cocaine;    

    // Create g4_tnt
    G4Material* g4_tnt = new G4Material("g4_tnt", 1.6 * CLHEP::g / CLHEP::cm3, 4);
    g4_tnt->AddElement(C, 7);
    g4_tnt->AddElement(H, 5);
    g4_tnt->AddElement(N, 3);
    g4_tnt->AddElement(O, 6);
    materials_["g4_tnt"] = g4_tnt;
  
    // Create g4_aspirin
    G4Material* g4_aspirin = new G4Material("g4_aspirin", 1.3 * CLHEP::g / CLHEP::cm3, 3);
    g4_aspirin->AddElement(C, 9);
    g4_aspirin->AddElement(H, 8);
    g4_aspirin->AddElement(O, 4);
    materials_["g4_aspirin"] = g4_aspirin;

    // Create g4_sucrose_sugar
    G4Material* g4_sucrose_sugar = new G4Material("g4_sucrose_sugar", 1.8 * CLHEP::g / CLHEP::cm3, 3);
    g4_sucrose_sugar->AddElement(C, 12);
    g4_sucrose_sugar->AddElement(H, 22);
    g4_sucrose_sugar->AddElement(O, 11);
    materials_["g4_sucrose_sugar"] = g4_sucrose_sugar;

    // Create g4_wood
    G4Material* g4_wood = new G4Material("g4_wood", 0.62 * CLHEP::g / CLHEP::cm3, 4);
    g4_wood->AddElement(H, 0.06);
    g4_wood->AddElement(C, 0.47);
    g4_wood->AddElement(O, 0.44);
    g4_wood->AddElement(N, 0.03);
    materials_["g4_wood"] = g4_wood;

    // Dense Hydrodgen
    G4Material* g4_hydrogen = new G4Material("g4_hydrogen", 1.0*CLHEP::g / CLHEP::cm3, 1);  
    g4_hydrogen->AddElement(H, 1);  
    materials_["g4_hydrogen"] = g4_hydrogen;  

    // Dense Nitrogen
    G4Material* g4_nitrogen = new G4Material("g4_hydrogen", 1.0*CLHEP::g / CLHEP::cm3, 1);  
    g4_nitrogen->AddElement(N, 1);  
    materials_["g4_nitrogen"] = g4_nitrogen;  

    // Dense Oxygen
    G4Material* g4_oxygen = new G4Material("g4_hydrogen", 1.0*CLHEP::g / CLHEP::cm3, 1);  
    g4_oxygen->AddElement(O, 1);  
    materials_["g4_oxygen"] = g4_oxygen;  

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
    // Material properties tables of CeBr3
    // Provided by Scionix
    G4double cebr3[] = {320, 330, 340, 345, 350, 355, 360, 365, 370, 371, 372, 375, 377.5, 380, 385,
                        389, 390, 391, 392, 395, 400, 405, 410, 415, 420, 425, 430, 435,   440};
    G4double conv_factor = 0.0012398;
    const G4int size = sizeof(cebr3) / sizeof(G4double);
    G4double cebr3_ch[size] = {0};
    for(auto i = 0; i < size; i++) {
        cebr3_ch[i] = (conv_factor) / cebr3[i];
    }
    // Provided by Scionix
    G4double cebr3_SC[] = {0.000, 0.000, 0.001, 0.030, 0.098, 0.199, 0.434, 0.779, 0.957, 0.961,
                           0.955, 0.906, 0.866, 0.849, 0.823, 0.804, 0.796, 0.781, 0.760, 0.696,
                           0.541, 0.370, 0.224, 0.120, 0.059, 0.029, 0.012, 0.004, 0.001};
    assert(sizeof(cebr3_SC) == sizeof(cebr3_ch));

    // info from
    // https://www.researchgate.net/figure/The-calculated-optical-constants-of-CeCl-3-and-CeBr-3-a-refractive-index-b_fig7_256855384
    G4double cebr3_RIND[] = {2.5,  2.48, 2.46, 2.55, 2.42, 2.40, 2.41, 2.42, 2.42, 2.42, 2.42, 2.43, 2.43, 2.43, 2.43,
                             2.44, 2.44, 2.44, 2.44, 2.44, 2.44, 2.45, 2.46, 2.47, 2.48, 2.49, 2.50, 2.51, 2.52};
    assert(sizeof(cebr3_RIND) == sizeof(cebr3_ch));
    // Info from https://www.advatech-uk.co.uk/cebr3.html -> Mass attenuation coeff = abs_length
    G4double cebr3_ABSL[] = {100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm,
                             100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm,
                             100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm,
                             100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm,
                             100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm,
                             100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm, 100 * CLHEP::cm};
    assert(sizeof(cebr3_ABSL) == sizeof(cebr3_ch));
    auto CeBr3_mt = new G4MaterialPropertiesTable();
    CeBr3_mt->AddProperty("FASTCOMPONENT", cebr3_ch, cebr3_SC, size);
    CeBr3_mt->AddProperty("RINDEX", cebr3_ch, cebr3_RIND, size);
    CeBr3_mt->AddProperty("ABSLENGTH", cebr3_ch, cebr3_ABSL, size);
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

void GeometryConstructionG4::check_overlaps() {
    G4PhysicalVolumeStore* phys_volume_store = G4PhysicalVolumeStore::GetInstance();
    LOG(TRACE) << "Checking overlaps";
    bool overlapFlag = false;
    // Release Geant4 output for better error description
    RELEASE_STREAM(G4cout);
    for(auto volume : (*phys_volume_store)) {
        overlapFlag = volume->CheckOverlaps(1000, 0., false) || overlapFlag;
    }
    // Supress again to prevent further complications
    SUPPRESS_STREAM(G4cout);
    if(overlapFlag) {
        LOG(ERROR) << "Overlapping volumes detected.";
    } else {
        LOG(INFO) << "No overlapping volumes detected.";
    }
}
