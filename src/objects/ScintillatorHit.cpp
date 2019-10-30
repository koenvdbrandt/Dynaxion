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
                                 unsigned int charge,
                                 double event_time,
                                 double charge_deposit,
                                 const MCParticle* mc_particle)
    : SensorCharge(std::move(local_position), std::move(global_position), type, charge, event_time),
      charge_deposit_(charge_deposit) {
    setMCParticle(mc_particle);
}

/**
 * @throws MissingReferenceException If the pointed object is not in scope
 *
 * Object is stored as TRef and can only be accessed if pointed object is in scope
 */
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

double ScintillatorHit::getChargeDeposit() const {
    return charge_deposit_;
}

void ScintillatorHit::print(std::ostream& out) const {
    out << "--- Deposited charge information\n";
    SensorCharge::print(out);
}
