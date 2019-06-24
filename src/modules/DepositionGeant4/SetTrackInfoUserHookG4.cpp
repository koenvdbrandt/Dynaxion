#include "SetTrackInfoUserHookG4.hpp"
#include "TrackInfoG4.hpp"

using namespace allpix;

void SetTrackInfoUserHookG4::PreUserTrackingAction(const G4Track* aTrack) {
    auto theTrack = const_cast<G4Track*>(aTrack); // NOLINT
    auto particle = aTrack->GetDefinition();

    // Unstable particles which are not the primary particle should be killed to stop the decay chain:
    if(!particle->GetPDGStable() && aTrack->GetTrackID() > 1) {
	if(particle->GetParticleName() != "neutron"){
        theTrack->SetTrackStatus(fStopAndKill);
    	}
    }
    if(aTrack->GetUserInformation() == nullptr) {
        auto trackInfo = track_info_mgr_ptr_->makeTrackInfo(aTrack);
        // Release ownership of the TrackInfoG4 instance
        theTrack->SetUserInformation(trackInfo.release());
    }
}

void SetTrackInfoUserHookG4::PostUserTrackingAction(const G4Track* aTrack) {
    auto userInfo = dynamic_cast<TrackInfoG4*>(aTrack->GetUserInformation());
    userInfo->finalizeInfo(aTrack);
    // Regain ownership of the TrackInfoG4, and remove it from the G4Track
    auto userInfoOwningPtr = std::unique_ptr<TrackInfoG4>(userInfo);
    auto theTrack = const_cast<G4Track*>(aTrack); // NOLINT
    theTrack->SetUserInformation(nullptr);
    track_info_mgr_ptr_->storeTrackInfo(std::move(userInfoOwningPtr));
}
