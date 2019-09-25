/**
 * @file
 * @brief Definition of object with digitized pixel hit
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_SCINTILLATOR_HIT_H
#define ALLPIX_SCINTILLATOR_HIT_H

#include <Math/DisplacementVector2D.h>

#include <TRef.h>

#include "DepositedCharge.hpp"
#include "MCParticle.hpp"
#include "Object.hpp"

namespace allpix {
    /**
     * @ingroup Objects
     * @brief Pixel triggered in an event after digitization
     */
    class ScintHit : public Object {
    public:
        /**
         * @brief Construct a digitized pixel hit
         * @param pixel Object holding the information of the pixel
         * @param time Timing of the occurence of the hit
         * @param signal Signal data produced by the digitizer
         * @param pixel_charge Optional pointer to the related pixel charge
         */
        ScintHit(std::string sensor, double time, double signal, const DepositedCharge* scint_hit = nullptr);

        /**
         * @brief Get the pixel hit
         * @return Pixel indices in the grid
         */
        const std::string& getSensor() const;

        /**
         * @brief Get the timing data of the hit
         * @return Timestamp
         */
        double getTime() const { return time_; }
        /**
         * @brief Get the signal data for the hit
         * @return Digitized signal
         */
        double getSignal() const { return signal_; }

        /**
         * @brief Get related pixel charge
         * @return Possible related pixel charge
         */
        const DepositedCharge* getScintHits() const;
        /**
         * @brief Get the Monte-Carlo particles resulting in this pixel hit
         * @return List of all related Monte-Carlo particles
         */
        const MCParticle* getMCParticles() const;

        /**
         * @brief Print an ASCII representation of PixelHit to the given stream
         * @param out Stream to print to
         */
        void print(std::ostream& out) const override;

        /**
         * @brief ROOT class definition
         */
        ClassDefOverride(ScintHit, 4);
        /**
         * @brief Default constructor for ROOT I/O
         */
        ScintHit() = default;

    private:
        std::string sensor_;
        double time_{};
        double signal_{};

        TRef scint_hits_;
        TRef mc_particle_{nullptr};
    };

    /**
     * @brief Typedef for message carrying pixel hits
     */
    using ScintHitMessage = Message<ScintHit>;
} // namespace allpix

#endif /* ALLPIX_SCINTILLATOR_HIT_H */
