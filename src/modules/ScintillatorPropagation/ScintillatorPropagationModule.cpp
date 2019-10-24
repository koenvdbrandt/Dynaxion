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
    config_.setDefault<double>("quantum_efficiency", 1.0);

    config_.setDefault<int>("charge_per_step", 10);
    config_.setDefault<double>("integration_time", Units::get(25, "ns"));
    config_.setDefault<bool>("output_plots", false);

    integration_time_ = config_.get<double>("integration_time");
    output_plots_ = config_.get<bool>("output_plots");

    config_.setDefault<double>("transit_time", Units::get(50, "ns"));
    config_.setDefault<double>("transit_time_spread", Units::get(10, "ns"));
    config_.setDefault<double>("rise_time", Units::get(10, "ns"));
    config_.setDefault<double>("rise_time_spread", Units::get(5, "ns"));

    /*     // Set default for charge carrier propagation:
        config_.setDefault<bool>("propagate_holes", false);
        if(config_.get<bool>("propagate_holes")) {
            propagate_type_ = CarrierType::HOLE;
            LOG(INFO) << "Holes are chosen for propagation. Electrons are therefore not propagated.";
        } else {
            propagate_type_ = CarrierType::ELECTRON;
        }*/

    // Parameterization variables from https://doi.org/10.1016/0038-1101(77)90054-5 (section 5.2)
    auto temperature = config_.get<double>("temperature", 293);
    electron_Vm_ = Units::get(1.53e9 * std::pow(temperature, -0.87), "cm/s");
    electron_Ec_ = Units::get(1.01 * std::pow(temperature, 1.55), "V/cm");
    electron_Beta_ = 2.57e-2 * std::pow(temperature, 0.66);

    hole_Vm_ = Units::get(1.62e8 * std::pow(temperature, -0.52), "cm/s");
    hole_Ec_ = Units::get(1.24 * std::pow(temperature, 1.68), "V/cm");
    hole_Beta_ = 0.46 * std::pow(temperature, 0.17);

    boltzmann_kT_ = Units::get(8.6173e-5, "eV/K") * temperature;

    config_.setDefault<bool>("ignore_magnetic_field", false);
}

void ScintillatorPropagationModule::init() {
    if(std::dynamic_pointer_cast<ScintillatorModel>(detector_->getModel()) == nullptr) {
        throw ModuleError("The detector " + detector_->getName() +
                          " is not a scintillator. Use other method of propagation");
    }

    /*   if(output_plots_) {
          auto time_bins =
              static_cast<int>(config_.get<double>("output_plots_range") / config_.get<double>("output_plots_step"));
          drift_time_histo = new TH1D("drift_time_histo",
                                      "Charge carrier arrival time;t[ns];charge carriers",
                                      time_bins,
                                      0.,
                                      config_.get<double>("output_plots_range"));
      if(output_plots_) {
          // Initialize output plot
          drift_time_histo_ = new TH1D("drift_time_histo",
                                       "Drift time;Drift time [ns];charge carriers",
                                       static_cast<int>(Units::convert(integration_time_, "ns") * 5),
                                       0,
                                       static_cast<double>(Units::convert(integration_time_, "ns")));
      }*/
}

void ScintillatorPropagationModule::run(unsigned int) {

    // Create vector of propagated charges to output
    std::vector<PixelCharge> pixel_charges;
    std::vector<PropagatedCharge> propagated_charges;
    std::vector<const PropagatedCharge*> prop_charges;

    double initial_charge = 0;
    double total_charge = 0;
    // double total_projected_charge = 0;
    // Loop over all deposits for propagation
    // for(auto& deposit : scintillator_message_->getData()) {
    for(auto& deposit : deposits_message_->getData()) {

        // auto position = deposit.getLocalPosition();
        auto charge = deposit.getCharge();
        initial_charge += charge;
        // DARK CURRENT

        auto time = deposit.getEventTime();
        // NEED TO ADD A FUNCTION WHICH INCREASES TIME BY THE PROPER AMOUNT

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
        double photo_electrons = 1;
        // Gain stages
        std::normal_distribution<double> gain_smearing(config_.get<double>("gain"), config_.get<double>("gain_smearing"));
        for(int i = 1; i <= config_.get<double>("gain_stages"); i++) {
            auto gain = gain_smearing(random_generator_);
            photo_electrons *= gain;
            LOG(DEBUG) << "photo-electrons after amplifying stage " << i << " with gain " << gain << " is "
                       << photo_electrons;
        }
        std::normal_distribution<double> transit_time_spread(config_.get<double>("transit_time"),
                                                             config_.get<double>("transit_time_spread"));
        std::normal_distribution<double> rise_time_spread(config_.get<double>("rise_time"),
                                                          config_.get<double>("rise_time_spread"));
        auto transit_time = transit_time_spread(random_generator_);
        //// Find a proper way to include rise time with transit time!!
        auto rise_time = rise_time_spread(random_generator_);
        (void)rise_time;
        time += transit_time;
        ///////
        auto end_of_sensor_local =
            (deposit.getLocalPosition().x(), deposit.getLocalPosition().y(), deposit.getLocalPosition().z());
        auto end_of_sensor_global =
            (deposit.getGlobalPosition().x(), deposit.getGlobalPosition().y(), deposit.getGlobalPosition().z());
        (void)end_of_sensor_local;
        (void)end_of_sensor_global;
        /////
        propagated_charges.emplace_back(
            deposit.getLocalPosition(), deposit.getGlobalPosition(), deposit.getType(), photo_electrons, time, &deposit);
        total_charge += photo_electrons;
        LOG(TRACE) << "prop photo_electrons info:";
        LOG(TRACE) << "Initial photon of energy : " << Units::display(charge, {"eV"});
        LOG(TRACE) << "propagated photo-electrons = " << photo_electrons;
        LOG(TRACE) << "transit time = " << transit_time;
        LOG(TRACE) << "arival time photo_electrons time = " << time;
    }

    // std::map<Pixel::Index, std::vector<const PropagatedCharge*>> pixel_map;
    Pixel::Index pixel_index(1, 1);
    unsigned long charge = 0;
    for(auto& prop_charge : propagated_charges) {
        charge += prop_charge.getCharge();
        prop_charges.emplace_back(&prop_charge);
    }
    if(charge != total_charge) {
        LOG(WARNING) << "charges not equal, look at code";
    }
    auto pixel = detector_->getPixel(pixel_index);
    pixel_charges.emplace_back(pixel, charge, prop_charges);
    /*
            if(output_plots_) {
                drift_time_histo_->Fill(drift_time, deposit.getCharge());
            }
    */

    LOG(INFO) << "Total photo electrons at the end of the Photo Multiplier for this event = " << total_charge
              << " (started with  " << deposits_message_->getData().size() << " photons)";

    // Create a new message with propagated charges
    auto pixel_message = std::make_shared<PixelChargeMessage>(std::move(pixel_charges), detector_);
    // Dispatch the message with propagated charges
    messenger_->dispatchMessage(this, pixel_message);
}

void ScintillatorPropagationModule::finalize() {
    if(output_plots_) {
        // Write output plot
        drift_time_histo_->Write();
    }
}