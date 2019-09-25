/**
 * @file
 * @brief Implementation of default digitization module
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "ScintillatorDigitizerModule.hpp"

#include "core/utils/unit.h"
#include "tools/ROOT.h"

#include <TFile.h>
#include <TH1D.h>
#include <TProfile.h>

#include "core/geometry/GeometryManager.hpp"
#include "objects/DepositedCharge.hpp"
#include "objects/ScintHit.hpp"

using namespace allpix;

ScintillatorDigitizerModule::ScintillatorDigitizerModule(Configuration& config,
                                                         Messenger* messenger,
                                                         std::shared_ptr<Detector> detector)
    : Module(config, std::move(detector)), messenger_(messenger), deposits_message_(nullptr) {
    // Enable parallelization of this module if multithreading is enabled
    enable_parallelization();

    // Require DepositedCharge message for single detector
    messenger_->bindSingle(this, &ScintillatorDigitizerModule::deposits_message_, MsgFlags::REQUIRED);

    // Seed the random generator with the global seed
    random_generator_.seed(getRandomSeed());

    //!!!!!!!!!!!!!!!!!!!!! SET DEFAULTS

    // Set defaults for config variables
    config_.setDefault<int>("electronics_noise", 1);
    config_.setDefault<double>("gain", 100000.0);
    config_.setDefault<double>("gain_smearing", 10000.0);
    config_.setDefault<int>("threshold", 70000);
    config_.setDefault<int>("threshold_smearing", 10000);

    config_.setDefault<int>("adc_resolution", 0);
    config_.setDefault<int>("adc_smearing", 300);
    config_.setDefault<double>("adc_offset", 0);
    config_.setDefault<double>("adc_slope", 10);

    config_.setDefault<bool>("output_plots", false);
    config_.setDefault<int>("output_plots_scale", 30);
    config_.setDefault<int>("output_plots_bins", 100);
}

void ScintillatorDigitizerModule::init() {
    auto type = geo_manager_->getDetectorType();
    (void)type;
    /* std::map<std::string, std::string> type = geo_manager_->getDetectorType();
    LOG(WARNING) << "Detector_type1" << getDetector()->getType();
    //LOG(WARNING) << "Detector_type2" << type[getDetector()->getType()];
    if(type[getDetector()->getType()] != "scintillator") {
        //throw ModuleError("Detector is no scintillator, use DefaultDigitizer");
    }*/
    // Conversion to ADC units requested:
    if(config_.get<int>("adc_resolution") > 31) {
        throw InvalidValueError(config_, "adc_resolution", "precision higher than 31bit is not possible");
    }
    if(config_.get<int>("adc_resolution") > 0) {
        LOG(INFO) << "Converting charge to ADC units, ADC resolution: " << config_.get<int>("adc_resolution")
                  << "bit, max. value " << ((1 << config_.get<int>("adc_resolution")) - 1);
    }

    if(config_.get<bool>("output_plots")) {
        LOG(TRACE) << "Creating output plots";

        // Plot axis are in kilo electrons - convert from framework units!
        int maximum = static_cast<int>(Units::convert(config_.get<int>("output_plots_scale"), "ke"));
        auto nbins = config_.get<int>("output_plots_bins");

        // Create histograms if needed
        h_pxq = new TH1D("pixelcharge", "raw pixel charge;pixel charge [ke];pixels", nbins, 0, maximum);
        h_pxq_noise = new TH1D("pixelcharge_noise", "pixel charge w/ el. noise;pixel charge [ke];pixels", nbins, 0, maximum);
        h_gain = new TH1D("gain", "applied gain; gain factor;events", 40, -20, 20);
        h_pxq_gain =
            new TH1D("pixelcharge_gain", "pixel charge w/ gain applied;pixel charge [ke];pixels", nbins, 0, maximum);
        h_thr = new TH1D("threshold", "applied threshold; threshold [ke];events", maximum, 0, maximum / 10);
        h_pxq_thr =
            new TH1D("pixelcharge_threshold", "pixel charge above threshold;pixel charge [ke];pixels", nbins, 0, maximum);
        h_pxq_adc_smear = new TH1D(
            "pixelcharge_adc_smeared", "pixel charge after ADC smearing;pixel charge [ke];pixels", nbins, 0, maximum);

        // Create final pixel charge plot with different axis, depending on whether ADC simulation is enabled or not
        if(config_.get<int>("adc_resolution") > 0) {
            int adcbins = ((1 << config_.get<int>("adc_resolution")) - 1);
            h_pxq_adc = new TH1D("pixelcharge_adc", "pixel charge after ADC;pixel charge [ADC];pixels", adcbins, 0, adcbins);
            h_calibration = new TH2D("charge_adc_calibration",
                                     "calibration curve of pixel charge to ADC units;pixel charge [ke];pixel charge [ADC]",
                                     nbins,
                                     0,
                                     maximum,
                                     adcbins,
                                     0,
                                     adcbins);
        } else {
            h_pxq_adc = new TH1D("pixelcharge_adc", "final pixel charge;pixel charge [ke];pixels", nbins, 0, maximum);
        }
    }
}

void ScintillatorDigitizerModule::run(unsigned int) {
    // Loop through all pixels with charges
    std::vector<ScintHit> hits;
    for(auto& scint_hits : deposits_message_->getData()) {
        auto sensor = getDetector()->getName();
        // auto pixel_index = pixel.getIndex();
        auto hit = static_cast<double>(scint_hits.getCharge());

        LOG(INFO) << "Received sensor " << sensor << ", hits " << hit;
        if(config_.get<bool>("output_plots")) {
            h_pxq->Fill(hit /*/ 1e3*/);
        }

        // Add electronics noise from Gaussian:
        std::normal_distribution<double> el_noise(0, config_.get<unsigned int>("electronics_noise"));
        hit += el_noise(random_generator_);

        LOG(DEBUG) << "Charge with noise: " << hit;
        if(config_.get<bool>("output_plots")) {
            h_pxq_noise->Fill(hit /*/ 1e3*/);
        }

        // Smear the gain factor, Gaussian distribution around "gain" with width "gain_smearing"
        std::normal_distribution<double> gain_smearing(config_.get<double>("gain"), config_.get<double>("gain_smearing"));
        double gain = gain_smearing(random_generator_);
        if(config_.get<bool>("output_plots")) {
            h_gain->Fill(gain);
        }

        // Apply the gain to the charge:
        hit *= gain;
        LOG(DEBUG) << "Hit after amplifier (gain): " << hit;
        if(config_.get<bool>("output_plots")) {
            h_pxq_gain->Fill(hit /*/ 1e3*/);
        }

        // Smear the threshold, Gaussian distribution around "threshold" with width "threshold_smearing"
        std::normal_distribution<double> thr_smearing(config_.get<unsigned int>("threshold"),
                                                      config_.get<unsigned int>("threshold_smearing"));
        double threshold = thr_smearing(random_generator_);
        if(config_.get<bool>("output_plots")) {
            h_thr->Fill(threshold /*/ 1e3*/);
        }

        // Discard charges below threshold:
        if(hit < threshold) {
            LOG(DEBUG) << "Below smeared threshold: " << hit << " < " << threshold;
            continue;
        }

        LOG(DEBUG) << "Passed threshold: " << hit << " > " << threshold;
        if(config_.get<bool>("output_plots")) {
            h_pxq_thr->Fill(hit /*/ 1e3*/);
        }

        // Simulate ADC if resolution set to more than 0bit
        if(config_.get<int>("adc_resolution") > 0) {
            // temporarily store old charge for histogramming:
            auto original_hit = hit;

            // Add ADC smearing:
            std::normal_distribution<double> adc_smearing(0, config_.get<unsigned int>("adc_smearing"));
            hit += adc_smearing(random_generator_);
            if(config_.get<bool>("output_plots")) {
                h_pxq_adc_smear->Fill(hit /*/ 1e3*/);
            }
            LOG(DEBUG) << "Smeared for simulating limited ADC sensitivity: " << hit;

            // Convert to ADC units and precision:
            hit = static_cast<double>(std::max(
                std::min(static_cast<int>((config_.get<double>("adc_offset") + hit) / config_.get<double>("adc_slope")),
                         (1 << config_.get<int>("adc_resolution")) - 1),
                0));
            LOG(DEBUG) << "Hits converted to ADC units: " << hit;

            if(config_.get<bool>("output_plots")) {
                h_calibration->Fill(original_hit /*/ 1e3*/, hit);
                h_pxq_adc->Fill(hit);
            }
        } else {
            // Fill the final pixel charge
            if(config_.get<bool>("output_plots")) {
                h_pxq_adc->Fill(hit /*/ 1e3*/);
            }
        }

        // Add the hit to the hitmap
        hits.emplace_back(sensor, 0, hit, &scint_hits);
    }

    // Output summary and update statistics
    LOG(INFO) << "Digitized " << hits.size() << " pixel hits";
    total_hits_ += hits.size();

    if(!hits.empty()) {
        // Create and dispatch hit message
        auto hits_message = std::make_shared<ScintHitMessage>(std::move(hits), getDetector());
        messenger_->dispatchMessage(this, hits_message);
    }
}

void ScintillatorDigitizerModule::finalize() {
    if(config_.get<bool>("output_plots")) {
        // Write histograms
        LOG(TRACE) << "Writing output plots to file";
        h_pxq->Write();
        h_pxq_noise->Write();
        h_gain->Write();
        h_pxq_gain->Write();
        h_thr->Write();
        h_pxq_thr->Write();
        h_pxq_adc->Write();

        if(config_.get<int>("adc_resolution") > 0) {
            h_pxq_adc_smear->Write();
            h_calibration->Write();
        }
    }

    LOG(INFO) << "Digitized " << total_hits_ << " pixel hits in total";
}
