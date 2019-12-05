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

#include "core/config/ConfigReader.hpp"
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
}

void ScintillatorHitConverterModule::init() {
    // Reading detector file
    std::ifstream file(config_.getPath("quantum_efficiency", true));
    std::string line;
    LOG(INFO) << "Reading Quantum Efficiency File '" << config_.get<std::string>("quantum_efficiency") << "'";
    while(std::getline(file, line)) {
        line = allpix::trim(line);
        if(isdigit(line[0]) == 0) {
            continue;
        }
        auto values = allpix::from_string<ROOT::Math::XYVector>(line);
        LOG(TRACE) << "Input(energy , efficiency)= " << values;
        wavelength_.push_back(values.x());
        if(!std::is_sorted(wavelength_.begin(), wavelength_.end())) {
            throw ModuleError("Quantum efficiency input is not sorted! Please sort the list of given wavelengths.");
        }
        efficiency_.push_back(values.y());
    }
}

void ScintillatorHitConverterModule::run(unsigned int) {

    if(efficiency_.empty()) {
        throw ModuleError("Quantum efficiency of specific detector is not given!");
    }

    // Create vector of propagated charges to output
    std::vector<DepositedCharge> deposited_charges;

    // Loop over all scint_hits
    for(auto& scint_hits : scint_hit_message_->getData()) {
        auto charge = scint_hits.getCharge();
        auto wavelength = (0.0012398) / charge;
        auto low = wavelength_.size();
        auto high = wavelength_.size();

        for(std::vector<double>::size_type i = 0; i != wavelength_.size(); i++) {
            if(wavelength_[i] <= wavelength && (low == wavelength_.size() || wavelength_[low] < wavelength_[i])) {
                low = i;
            } else if(wavelength <= wavelength_[i] && (high == wavelength_.size() || wavelength_[i] < wavelength_[high])) {
                high = i;
            }
        }
        // Give an error if the provided quantum efficiency file does not include the wavelength range of the photons
        // produced
        if(low == wavelength_.size()) {
            throw ModuleError("Wavelength of incident photon ('" + to_string(wavelength) +
                              "nm') is lower than specified in the quantum efficiency file!");
        } else if(high == wavelength_.size()) {
            throw ModuleError("Wavelength of incident photon ('" + to_string(wavelength) +
                              "nm') is higher than specified in the quantum efficiency file!");
        }

        auto quantum_efficiency = (wavelength_[low] * efficiency_[low] + wavelength_[high] * efficiency_[high]) /
                                  (wavelength_[high] + wavelength_[low]);

        LOG(TRACE) << "quantum_efficiency of the photocathode = " << quantum_efficiency << " for wavelenght " << wavelength
                   << "nm";

        // FIXME:: Funtion of QE dependant on wavelength
        // auto quantum_efficiency = 1;
        if(quantum_efficiency > 1.0) {
            throw ModuleError("Quantum efficiency must be between 0.0 and 1.0");
        }
        auto rand = std::uniform_real_distribution<double>(0, 1)(random_generator_);
        LOG(DEBUG) << "Random variable = " << rand;

        if(quantum_efficiency < rand) {
            LOG(DEBUG) << "Optical Photon did not get recognized by the Photocathode";
            continue;
        }
        unsigned int transferred_charge = 1;
        deposited_charges.emplace_back(scint_hits.getLocalPosition(),
                                       scint_hits.getGlobalPosition(),
                                       scint_hits.getType(),
                                       transferred_charge,
                                       scint_hits.getEventTime());
    }
    total_received_ += scint_hit_message_->getData().size();
    total_propagated_ += deposited_charges.size();
    LOG(TRACE) << "Photocathode Hits for this event: " << scint_hit_message_->getData().size()
               << " (propagated: " << deposited_charges.size() << ", "
               << (deposited_charges.size() * 100 / scint_hit_message_->getData().size()) << "%)";
    // Create a new message with deposited charges
    auto deposited_charge_message = std::make_shared<DepositedChargeMessage>(std::move(deposited_charges), detector_);
    // Dispatch the message with deposited charges
    messenger_->dispatchMessage(this, deposited_charge_message);
}
void ScintillatorHitConverterModule::finalize() {
    if(total_received_ > 0) {
        LOG(INFO) << "Total Photocathode Hits: " << total_received_ << " (propagated: " << total_propagated_ << ", "
                  << (total_propagated_ * 100 / total_received_) << "%)";
    } else {
        LOG(WARNING) << "No Photocathode hits registered";
    }
}
