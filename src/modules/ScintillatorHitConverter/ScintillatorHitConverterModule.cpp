/**
 * @file
 * @brief Definition of default digitization module
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "ScintillatorHitConverterModule.hpp"

#include <fstream>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>

#include "core/config/exceptions.h"
#include "core/geometry/ScintillatorModel.hpp"
#include "core/messenger/Messenger.hpp"
#include "core/utils/log.h"
#include "objects/DepositedCharge.hpp"
#include "objects/PixelCharge.hpp"
#include "objects/ScintillatorHit.hpp"

using namespace allpix;

ScintillatorHitConverterModule::ScintillatorHitConverterModule(Configuration& config,
                                                               Messenger* messenger,
                                                               std::shared_ptr<Detector> detector)
    : Module(config, detector), messenger_(messenger), detector_(std::move(detector)) {
    // Save detector model
    model_ = detector_->getModel();

    random_generator_.seed(getRandomSeed());

    // Require deposits message for single detector
    messenger_->bindSingle(this, &ScintillatorHitConverterModule::scint_hit_message_, MsgFlags::REQUIRED);

    // Set default value for config variables
    config_.setDefault<double>("quantum_efficiency", 1.0);
}

void ScintillatorHitConverterModule::init() {
    //
    // DO I STILL NEED THIS?
    //
    // if(std::dynamic_pointer_cast<ScintillatorModel>(detector_->getModel()) == nullptr) {
    //    throw ModuleError("The detector " + detector_->getName() +
    //                     " is not a scintillator. Use other method of propagation");
    //}
}

void ScintillatorHitConverterModule::run(unsigned int) {

    // Create vector of propagated charges to output
    std::vector<DepositedCharge> deposited_charges;

    // Loop over all scint_hits
    for(auto& scint_hits : scint_hit_message_->getData()) {

        auto charge = scint_hits.getChargeDeposit();
        (void)charge;
        // FIXME:: Funtion of QE dependant on wavelength
        auto quantum_efficiency = config_.get<double>("quantum_efficiency");
        if(quantum_efficiency > 1.0) {
            throw ModuleError("Quantum efficiency must be between 0.0 and 1.0");
        }
        auto rand = std::uniform_real_distribution<double>(0, 1)(random_generator_);
        LOG(DEBUG) << "Random variable = " << rand;

        if(quantum_efficiency < rand) {
            LOG(DEBUG) << "Optical Photon did not get recognized by the Scintillator";
            continue;
        }
        unsigned int transferred_charge = 1;
        deposited_charges.emplace_back(scint_hits.getLocalPosition(),
                                       scint_hits.getGlobalPosition(),
                                       scint_hits.getType(),
                                       transferred_charge,
                                       scint_hits.getEventTime());
    }
    LOG(INFO) << "Total scintillator hits: " << scint_hit_message_->getData().size()
              << " (lost: " << scint_hit_message_->getData().size() - scint_hit_message_->getData().size() << ", "
              << (deposited_charges.size() * 100 / scint_hit_message_->getData().size()) << "%)";
    LOG(DEBUG) << "Total count of propagated charge carriers: " << deposited_charges.size();
    // Create a new message with deposited charges
    auto deposited_charge_message = std::make_shared<DepositedChargeMessage>(std::move(deposited_charges), detector_);
    // Dispatch the message with deposited charges
    messenger_->dispatchMessage(this, deposited_charge_message);
}
void ScintillatorHitConverterModule::finalize() {}