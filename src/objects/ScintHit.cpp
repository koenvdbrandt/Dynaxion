/**
 * @file
 * @brief Implementation of object with digitized pixel hit
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "ScintHit.hpp"

#include <set>

#include "DepositedCharge.hpp"
#include "PropagatedCharge.hpp"
#include "exceptions.h"

using namespace allpix;

ScintHit::ScintHit(std::string sensor, double time, double signal, const DepositedCharge* scint_hits)
    : sensor_(std::move(sensor)), time_(time), signal_(signal) {
    scint_hits_ = const_cast<DepositedCharge*>(scint_hits); // NOLINT
    if(scint_hits != nullptr) {
        mc_particle_ = scint_hits->mc_particle_;
    }
}

const std::string& ScintHit::getSensor() const {
    return sensor_;
}

/**
 * @throws MissingReferenceException If the pointed object is not in scope
 *
 * Object is stored as TRef and can only be accessed if pointed object is in scope
 */
const DepositedCharge* ScintHit::getScintHits() const {
    auto scint_hit = dynamic_cast<DepositedCharge*>(scint_hits_.GetObject());
    if(scint_hit == nullptr) {
        throw MissingReferenceException(typeid(*this), typeid(DepositedCharge));
    }
    return scint_hit;
}

/**
 * @throws MissingReferenceException If the pointed object is not in scope
 *
 * MCParticles can only be fetched if the full history of objects are in scope and stored
 */
const MCParticle* ScintHit::getMCParticles() const {
    auto mc_particle = dynamic_cast<MCParticle*>(mc_particle_.GetObject());
    if(mc_particle == nullptr) {
        throw MissingReferenceException(typeid(*this), typeid(MCParticle));
    }
    return mc_particle;
}

void ScintHit::print(std::ostream& out) const {
    out << "ScintHit " << this->getSensor() << ", " << this->getSignal() << ", " << this->getTime();
}
