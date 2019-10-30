/**
 * @file
 * @brief Definition of deposited charge object
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_SCINTILLATOR_HIT_H
#define ALLPIX_SCINTILLATOR_HIT_H

#include <TRef.h>

#include "MCParticle.hpp"
#include "SensorCharge.hpp"

namespace allpix {
    /**
     * @ingroup Objects
     * @brief Charge deposit in sensor of detector
     */
    class ScintillatorHit : public SensorCharge {

    public:
        /**
         * @brief Construct a charge deposit
         * @param local_position Local position of the deposit in the sensor
         * @param global_position Global position of the propagated set of charges in the sensor
         * @param type Type of the carrier
         * @param charge Total charge of the deposit
         * @param event_time Time of deposition after event start
         * @param mc_particle Optional pointer to related MC particle
         */
        ScintillatorHit(ROOT::Math::XYZPoint local_position,
                        ROOT::Math::XYZPoint global_position,
                        CarrierType type,
                        unsigned int charge,
                        double event_time,
                        double charge_deposit,
                        const MCParticle* mc_particle = nullptr);

        /**
         * @brief Get related Monte-Carlo particle
         * @return Pointer to possible Monte-Carlo particle
         */
        const MCParticle* getMCParticle() const;

        /**
         * @brief Set the Monte-Carlo particle
         * @param mc_particle The Monte-Carlo particle
         * @warning Special method because MCParticle is only known after deposit creation, should not be replaced later.
         */
        void setMCParticle(const MCParticle* mc_particle);

        double getChargeDeposit() const;

        /**
         * @brief Print an ASCII representation of DepositedCharge to the given stream
         * @param out Stream to print to
         */
        void print(std::ostream& out) const override;

        /**
         * @brief ROOT class definition
         */
        ClassDefOverride(ScintillatorHit, 2);
        /**
         * @brief Default constructor for ROOT I/O
         */
        ScintillatorHit() = default;

    private:
        TRef mc_particle_;
        double charge_deposit_;
    };

    /**
     * @brief Typedef for message carrying deposits
     */
    using ScintillatorHitMessage = Message<ScintillatorHit>;
} // namespace allpix

#endif /* ALLPIX_SCINTILLATOR_HIT_H */
