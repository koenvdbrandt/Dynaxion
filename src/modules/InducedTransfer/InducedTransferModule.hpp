/**
 * @file
 * @brief Definition of InducedTransfer module
 * @copyright Copyright (c) 2019-2020 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <string>

#include "core/config/Configuration.hpp"
#include "core/geometry/DetectorModel.hpp"
#include "core/messenger/Messenger.hpp"
#include "core/module/Module.hpp"

#include "objects/PropagatedCharge.hpp"

namespace allpix {
    /**
     * @ingroup Modules
     * @brief Module to calculate the total induced charge from propagated charge carriers and to assign them to pixels
     * @note This module supports parallelization
     *
     * This module calculates the total induced charge by a charge carrier via the Shockley-Ramo theorem by comparing the
     * weighting potential at the initial and final position. The initial position is retrieved via the DepositedCharge
     * object in the history. The total induced charge is calculated per pixel and published as PixelCharge object.
     *
     * This module requires a weighting potential and only works properly of both electrons and holes are present among the
     * propagated charge carriers.
     */
    class InducedTransferModule : public Module {
    public:
        /**
         * @brief Constructor for the InducedTransfer module
         * @param config Configuration object as retrieved from the steering file
         * @param messenger Pointer to the messenger object to allow binding to messages on the bus
         * @param detector Pointer to the detector for this module instance
         */
        InducedTransferModule(Configuration& config, Messenger* messenger, const std::shared_ptr<Detector>& detector);

        /**
         * @brief Initial check for the presence of a weighting potential
         */
        void init() override;

        /**
         * @brief Calculation of the individual total induced charge and combination for all pixels
         */
        void run(unsigned int) override;

    private:
        Messenger* messenger_;
        std::shared_ptr<Detector> detector_;
        std::shared_ptr<DetectorModel> model_;

        // Message containing the propagated charges
        std::shared_ptr<PropagatedChargeMessage> propagated_message_;

        // Induction matrix size in number of pixels along x and y
        ROOT::Math::DisplacementVector2D<ROOT::Math::Cartesian2D<int>> matrix_;
    };
} // namespace allpix
