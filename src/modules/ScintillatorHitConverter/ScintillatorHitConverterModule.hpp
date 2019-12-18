/**
 * @file
 * @brief Definition of default digitization module
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_SCINTILLATOR_HIT_CONVERTER_MODULE_H
#define ALLPIX_SCINTILLATOR_HIT_CONVERTER_MODULE_H

#include <map>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include <TH1D.h>

#include "core/config/Configuration.hpp"
#include "core/geometry/DetectorModel.hpp"
#include "core/messenger/Messenger.hpp"
#include "core/module/Module.hpp"

#include "objects/DepositedCharge.hpp"
#include "objects/Pixel.hpp"
#include "objects/PixelCharge.hpp"
#include "objects/ScintillatorHit.hpp"

namespace allpix {
    /**
     * @ingroup Modules
     * @brief Module to project created electrons onto the sensor surface including diffusion
     *
     * The electrons from the deposition message are projected onto the sensor surface as a simple propagation method.
     * Diffusion is added by approximating the drift time and drawing a random number from a 2D gaussian distribution of the
     * calculated width.
     */
    class ScintillatorHitConverterModule : public Module {
    public:
        /**
         * @brief Constructor for this detector-specific module
         * @param config Configuration object for this module as retrieved from the steering file
         * @param messenger Pointer to the messenger object to allow binding to messages on the bus
         * @param detector Pointer to the detector for this module instance
         */
        ScintillatorHitConverterModule(Configuration& config, Messenger* messenger, std::shared_ptr<Detector> detector);

        /**
         * @brief Initialize - create plots if needed
         */
        void init() override;

        /**
         * @brief Projection of the electrons to the surface
         */
        void run(unsigned int) override;

        /**
         * @brief Write plots if needed
         */
        void finalize() override;

    private:
        // std::list<Configuration> quantum_efficiency_;
        Messenger* messenger_;
        std::shared_ptr<const Detector> detector_;

        std::shared_ptr<DetectorModel> model_;

        std::vector<double> efficiency_{};
        std::vector<double> wavelength_{};

        // Output plot for drift time
        bool output_plots_;
        TH1D* wavelength_before_;
        TH1D* wavelength_after_;
        TH1D* photo_electrons_;
        // Random generator for diffusion calculation
        std::mt19937_64 random_generator_;

        // Deposits for the bound detector in this event
        size_t total_propagated_{};
        size_t total_received_{};
        std::shared_ptr<ScintillatorHitMessage> scint_hit_message_;
    };
} // namespace allpix

#endif /* ALLPIX_SCINTILLATOR_HIT_CONVERTER_MODULE_H */
