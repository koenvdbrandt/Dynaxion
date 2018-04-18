/**
 * @file
 * @brief Implementation of the TrackInfoManager class
 * @copyright Copyright (c) 2018 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "TrackInfoManager.hpp"

using namespace allpix;

TrackInfoManager::TrackInfoManager() : counter_(1) {}

std::unique_ptr<TrackInfoG4> TrackInfoManager::makeTrackInfo(const G4Track* const track) {
    auto i = counter_++;
    auto G4ParentID = track->GetParentID();
    auto parent_track_id = G4ParentID == 0 ? G4ParentID : g4_to_custom_id_.at(G4ParentID);
    g4_to_custom_id_[track->GetTrackID()] = i;
    return std::make_unique<TrackInfoG4>(i, parent_track_id, track);
}

void TrackInfoManager::setTrackInfoToBeStored(int track_id) {
    auto element = std::find(to_store_track_ids_.begin(), to_store_track_ids_.end(), track_id);
    // If track id is not present we add it, otherwise skip as we only need each track once
    if(element == to_store_track_ids_.end()) {
        to_store_track_ids_.emplace_back(track_id);
    }
}

void TrackInfoManager::storeTrackInfo(std::unique_ptr<MCTrack> the_track) {
    auto track_id = the_track->getTrackID();
    auto element = std::find(to_store_track_ids_.begin(), to_store_track_ids_.end(), track_id);
    if(element != to_store_track_ids_.end()) {
        stored_tracks_.push_back(*std::move(the_track));
        to_store_track_ids_.erase(element);
    }
}

void TrackInfoManager::resetTrackInfoManager() {
    counter_ = 1;
    stored_tracks_.clear();
    to_store_track_ids_.clear();
    g4_to_custom_id_.clear();
}

void TrackInfoManager::dispatchMessage(Module* module, Messenger* messenger) {
    setAllTrackParents();
    auto mc_track_message = std::make_shared<MCTrackMessage>(std::move(stored_tracks_));
    messenger->dispatchMessage(module, mc_track_message);
}

MCTrack const* TrackInfoManager::findMCTrack(int track_id) const {
    auto element = std::find_if(stored_tracks_.begin(), stored_tracks_.end(), [&track_id](MCTrack const& track) {
        return (track.getTrackID() == track_id);
    });
    return (element != stored_tracks_.end()) ? &(*element) : nullptr;
}

void TrackInfoManager::setAllTrackParents() {
    for(auto& track : stored_tracks_) {
        track.setParent(findMCTrack(track.getParentTrackID()));
    }
}
