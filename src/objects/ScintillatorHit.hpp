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

#include <Math/Point3D.h>
#include <TRef.h>

#include "MCParticle.hpp"

#include "Object.hpp"

namespace allpix {
    /**
     * @ingroup Objects
     * @brief Charge deposit in sensor of detector
     */
    class ScintillatorHit : public Object {

    public:
        /**
         * @brief Construct a charge deposit
         * @param local_position Local position of the deposit in the sensor
         * @param global_position Global position of the propagated set of charges in the sensor
         * @param type Type of the carrier
         * @param energy Total energy of the deposit
         * @param wavelength Total wavelength of the deposit
         * @param event_time Time of deposition after event start
         * @param mc_particle Optional pointer to related MC particle
         */
        ScintillatorHit(ROOT::Math::XYZPoint local_position,
                        ROOT::Math::XYZPoint global_position,
                        CarrierType type,
                        double energy,
                        double wavelength,
                        double emission_time,
                        double detection_time,
                        const MCParticle* mc_particle = nullptr);
        /**
         * @brief Get local position of the set of charges in the sensor
         * @return Local position of charges
         */
        ROOT::Math::XYZPoint getLocalPosition() const;

        /**
         * @brief Get the global position of the set of charges in the sensor
         */
        ROOT::Math::XYZPoint getGlobalPosition() const;

        /**
         * @brief Get the type of charge carrier
         * @return Type of charge carrier
         */
        CarrierType getType() const;
        /**
         * @brief Get energy of the optical photon
         * @return Total charge stored
         */
        double getEnergy() const;
        /**
         * @brief Get wavelength of the optical photon
         * @return Total charge stored
         */
        double getWavelength() const;
        /**
         * @brief Get time of emmision after start of event
         * @return Time from start event
         */
        double getEmissionTime() const;
        /**
         * @brief Get time of detection after start of event
         * @return Time from start event
         */
        double getDetectionTime() const;
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
        /**
         * @brief Print an ASCII representation of DepositedCharge to the given stream
         * @param out Stream to print to
         */
        void print(std::ostream& out) const override;

        /**
         * @brief ROOT class definition
         */
        ClassDefOverride(ScintillatorHit, 1);
        /**
         * @brief Default constructor for ROOT I/O
         */
        ScintillatorHit() = default;

    private:
        ROOT::Math::XYZPoint local_position_;
        ROOT::Math::XYZPoint global_position_;

        CarrierType type_{};
        double energy_{};
        double wavelength_{};
        double emission_time_{};
        double detection_time_{};
        TRef mc_particle_;
    };

    /**
     * @brief Typedef for message carrying deposits
     */
    using ScintillatorHitMessage = Message<ScintillatorHit>;
} // namespace allpix

#endif /* ALLPIX_SCINTILLATOR_HIT_H */
