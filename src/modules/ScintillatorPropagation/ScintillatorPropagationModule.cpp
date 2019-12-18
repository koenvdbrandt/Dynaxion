/**
 * @file
 * @brief Definition of default digitization module
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "ScintillatorPropagationModule.hpp"

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

ScintillatorPropagationModule::ScintillatorPropagationModule(Configuration& config,
                                                             Messenger* messenger,
                                                             std::shared_ptr<Detector> detector)
    : Module(config, detector), messenger_(messenger), detector_(std::move(detector)) {
    // Save detector model
    model_ = detector_->getModel();

    random_generator_.seed(getRandomSeed());

    // Require deposits message for single detector
    // messenger_->bindSingle(this, &ScintillatorPropagationModule::scintillator_message_, MsgFlags::REQUIRED);
    messenger_->bindSingle(this, &ScintillatorPropagationModule::deposits_message_, MsgFlags::REQUIRED);

    // Set default value for config variables
    config_.setDefault<int>("gain", 10);
    config_.setDefault<double>("gain_smearing", 0.0);
    config_.setDefault<int>("gain_stages", 10);

    config_.setDefault<double>("transit_time", Units::get(50, "ns"));
    config_.setDefault<double>("transit_time_spread", Units::get(10, "ns"));
    config_.setDefault<double>("rise_time", Units::get(10, "ns"));
    config_.setDefault<double>("rise_time_spread", Units::get(2, "ns"));
    config_.setDefault<ROOT::Math::XYZVector>("pm_prop", ROOT::Math::XYZVector());

    output_plots_ = config_.get<bool>("output_plots");
}

void ScintillatorPropagationModule::init() {
    if(std::dynamic_pointer_cast<ScintillatorModel>(detector_->getModel()) == nullptr) {
        throw ModuleError("The detector " + detector_->getName() +
                          " is not a scintillator. Use other method of propagation");
    }

    if(output_plots_) {
        auto max_electrons = static_cast<int>(config_.get<double>("output_scale_electrons", 10000));
        auto nbins = 5 * max_electrons;
        photo_electrons_before_ =
            new TH1D("photo_electrons_before", "Photo electrons; Photo electrons;events", nbins, 0., max_electrons);
        gain_stages_ = config_.get<int>("gain_stages");
        photo_electrons_after_ = new TH1D("photo_electrons_after",
                                          "Photo electrons; Photo electrons;events",
                                          nbins,
                                          0.,
                                          max_electrons * gain_stages_ * 10);
    }
}

void ScintillatorPropagationModule::run(unsigned int) {

    // Create vector of propagated charges to output
    std::vector<PropagatedCharge> propagated_charges;

    auto initial_deposits = deposits_message_->getData().size();
    unsigned int total_charge = 0;

    // Loop over all deposits for propagation
    for(auto& deposit : deposits_message_->getData()) {

        // FIXME DARK CURRENT
        auto charge = deposit.getCharge();
        // Gain stages
        std::normal_distribution<double> gain_smearing(config_.get<double>("gain"), config_.get<double>("gain_smearing"));
        for(int i = 1; i <= gain_stages_; i++) {
            auto gain = static_cast<unsigned int>(std::round(gain_smearing(random_generator_)));
            charge *= gain;
            LOG(DEBUG) << "photo-electrons after amplifying stage " << i << " with gain " << gain << " is " << charge;
        }

        // NEED TO ADD A FUNCTION WHICH INCREASES TIME BY THE PROPER AMOUNT
        auto time = deposit.getEventTime();
        std::normal_distribution<double> transit_time_spread(config_.get<double>("transit_time"),
                                                             config_.get<double>("transit_time_spread"));
        std::normal_distribution<double> rise_time_spread(config_.get<double>("rise_time"),
                                                          config_.get<double>("rise_time_spread"));
        auto transit_time = transit_time_spread(random_generator_);
        //// Find a proper way to include rise time with transit time!!
        auto rise_time = rise_time_spread(random_generator_);
        (void)rise_time;
        time += transit_time;
        // Position Propagation
        // FIXME:: HAVE TO FIND SOMETHING WITH ORIENTATION OF PM
        auto pm_prop = config_.get<ROOT::Math::XYZVector>("pm_prop");
        LOG(TRACE) << "Location propagation in PM = " << pm_prop;
        deposit.getLocalPosition() += pm_prop;
        deposit.getGlobalPosition() += pm_prop;
        /////
        propagated_charges.emplace_back(
            deposit.getLocalPosition(), deposit.getGlobalPosition(), deposit.getType(), charge, time, &deposit);
        total_charge += charge;
        LOG(TRACE) << "prop charge info:";
        LOG(TRACE) << "propagated charges = " << charge;
        LOG(TRACE) << "transit time = " << transit_time;
        LOG(TRACE) << "arival time charge time = " << time;
    }

    LOG(INFO) << "Total photo electrons at the end of the Photo Multiplier for this event = " << total_charge
              << " (started with  " << initial_deposits << " photons)";

    if(output_plots_) {
        photo_electrons_before_->Fill(initial_deposits);
        photo_electrons_after_->Fill(total_charge);
    }

    // Create a new message with propagated charges
    auto prop_chare_message = std::make_shared<PropagatedChargeMessage>(std::move(propagated_charges), detector_);
    // Dispatch the message with propagated charges
    messenger_->dispatchMessage(this, prop_chare_message);
}

void ScintillatorPropagationModule::finalize() {
    if(output_plots_) {
        // Write output plot
        photo_electrons_before_->Write();
        photo_electrons_after_->Write();
    }
}