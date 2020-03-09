/**
 * @file
 * @brief Defines a user hook for Geant4 to assign custom track information via TrackInfoG4 objects. This includes custom
 * (unique) track ids.
 *
 * @copyright Copyright (c) 2018-2020 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef SetTrackInfoUserHookG4_H
#define SetTrackInfoUserHookG4_H 1

#include "G4Track.hh"
#include "G4UserTrackingAction.hh"
#include "core/config/Configuration.hpp"

#include "TrackInfoManager.hpp"

namespace allpix {
    /**
     * @brief Assigns every G4Track a TrackInfoG4 which carries various inforamtion, including the custom track id
     */
    class SetTrackInfoUserHookG4 : public G4UserTrackingAction {
    public:
        /**
         * @brief Constructor taking a TrackInfoManager*
         * @param track_info_mgr_ptr Pointer to TrackInfoManager which must be used to create the TrackInfoG4 instances
         */
        explicit SetTrackInfoUserHookG4(TrackInfoManager* track_info_mgr_ptr, double decay_cutoff_time)
            : track_info_mgr_ptr_(track_info_mgr_ptr), decay_cutoff_time_(decay_cutoff_time){};

        /**
         * @brief Default destructor
         */
        ~SetTrackInfoUserHookG4() override = default;

        /**
         * @brief Called for every G4Track at beginning
         * @param aTrack The pointer to the G4Track for which this routine is called
         */
        void PreUserTrackingAction(const G4Track* aTrack) override;

        /**
         * @brief Called for every G4Track at end
         * @param aTrack The pointer to the G4Track for which this routine is called
         */
        void PostUserTrackingAction(const G4Track* aTrack) override;

    private:
        // Raw ptr to track info manager to create instances of TrackInfoG4
        TrackInfoManager* track_info_mgr_ptr_;
        double decay_cutoff_time_;
    };

} // namespace allpix
#endif /* SetTrackInfoUserHookG4_H */
