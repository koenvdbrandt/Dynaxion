#include "SetTrackInfoUserHookG4.hpp"
#include "TrackInfoG4.hpp"

using namespace allpix;

void SetTrackInfoUserHookG4::PreUserTrackingAction(const G4Track* aTrack) {
    auto theTrack = const_cast<G4Track*>(aTrack); // NOLINT
    auto particle = aTrack->GetDefinition();
    /* Unstable particles which are not the primary particle should be killed to stop the decay chain:
        && particle->GetPDGEncoding() != 1000581380 && particle->GetPDGEncoding() != 1000581390 
        && particle->GetPDGEncoding() != 1000350820 && particle->GetPDGEncoding() != 1000571399 
        && particle->GetPDGEncoding() != 1000340810 && particle->GetPDGEncoding() != 1000350800 
        && particle->GetPDGEncoding() != 1000350780 && particle->GetPDGEncoding() != 1000330760
        && particle->GetPDGEncoding() != 1000581410 && particle->GetPDGEncoding() != 1000581360
        && particle->GetPDGEncoding() != 1000340790 && particle->GetPDGEncoding() != 1000360809
        && particle->GetPDGEncoding() != 1000360829 && particle->GetPDGEncoding() != 1000581430
        && particle->GetPDGEncoding() != 1000591419 && particle->GetPDGEncoding() != 1000340769
        && particle->GetPDGEncoding() != 1000330780 && particle->GetPDGEncoding() != 1000571400
        && particle->GetPDGEncoding() != 1000581350 && particle->GetPDGEncoding() != 1000571400
        && particle->GetPDGEncoding() != 1000581420 && particle->GetPDGEncoding() != 1000340789 */   
    if(!particle->GetPDGStable() && aTrack->GetTrackID() > 1 && particle->GetPDGEncoding() != 2112){
        theTrack->SetTrackStatus(fStopAndKill);
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
