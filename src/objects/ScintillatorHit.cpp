/**
 * @file
 * @brief Implementation of deposited charge object
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "ScintillatorHit.hpp"

#include "exceptions.h"

using namespace allpix;

ScintillatorHit::ScintillatorHit(ROOT::Math::XYZPoint local_position,
                                 ROOT::Math::XYZPoint global_position,
                                 CarrierType type,
                                 double energy,
                                 double wavelength,
                                 double emission_time,
                                 double detection_time,
                                 const MCParticle* mc_particle)

    : local_position_(std::move(local_position)), global_position_(std::move(global_position)), type_(type), energy_(energy),
      wavelength_(wavelength), emission_time_(emission_time), detection_time_(detection_time) {
    setMCParticle(mc_particle);
}

ROOT::Math::XYZPoint ScintillatorHit::getLocalPosition() const {
    return local_position_;
}

ROOT::Math::XYZPoint ScintillatorHit::getGlobalPosition() const {
    return global_position_;
}

CarrierType ScintillatorHit::getType() const {
    return type_;
}

double ScintillatorHit::getEnergy() const {
    return energy_;
}

double ScintillatorHit::getWavelength() const {
    return wavelength_;
}

double ScintillatorHit::getEmissionTime() const {
    return emission_time_;
}

double ScintillatorHit::getDetectionTime() const {
    return detection_time_;
}

void ScintillatorHit::print(std::ostream& out) const {
    out << "Type: " << (type_ == CarrierType::ELECTRON ? "\"e\"" : "\"h\"") << "\nEnergy " << energy_
        << " eV (Wavelength = " << wavelength_ << "nm"
        << "\nLocal Position: (" << local_position_.X() << ", " << local_position_.Y() << ", " << local_position_.Z()
        << ") mm\n"
        << "Global Position: (" << global_position_.X() << ", " << global_position_.Y() << ", " << global_position_.Z()
        << ") mm\n";
}

const MCParticle* ScintillatorHit::getMCParticle() const {
    auto mc_particle = dynamic_cast<MCParticle*>(mc_particle_.GetObject());
    if(mc_particle == nullptr) {
        throw MissingReferenceException(typeid(*this), typeid(MCParticle));
    }
    return mc_particle;
}

void ScintillatorHit::setMCParticle(const MCParticle* mc_particle) {
    mc_particle_ = const_cast<MCParticle*>(mc_particle); // NOLINT
}
