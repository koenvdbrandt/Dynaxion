/**
 * @file
 * @brief Implementation of [ParticleDistributions] module
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "ParticleDistributionsModule.hpp"

#include <string>
#include <utility>

#include "core/utils/log.h"

using namespace allpix;

ParticleDistributionsModule::ParticleDistributionsModule(Configuration& config,
                                                         Messenger* messenger,
                                                         GeometryManager* geo_manager)
    : Module(config), messenger_(messenger), geo_manager_(geo_manager) {

    // ... Implement ... (Typically bounds the required messages and optionally sets configuration defaults)
    // Input required by this module
    messenger->bindSingle(this, &ParticleDistributionsModule::message_, MsgFlags::REQUIRED);
}

void ParticleDistributionsModule::init() {

    std::vector<std::shared_ptr<Detector>> detectors = geo_manager_->getDetectors();
    for(auto& detector : detectors) {
        // Get the detector name
        std::string detectorName = detector->getName();
        LOG(DEBUG) << "Detector with name " << detectorName;
    }
    energy_distribution_ = new TH1F("energy_distribution", "energy_distribution", 1000, 0, 15);
    zx_distribution_ = new TH2F("zx_distribution", "zx_distribution", 100, -1, 1, 100, -1, 1);
    zy_distribution_ = new TH2F("zy_distribution", "zy_distribution", 100, -1, 1, 100, -1, 1);
    xyz_distribution_ = new TH3F("xyz_distribution", "xyz_distribution", 100, -1, 1, 100, -1, 1, 100, -1, 1);
    xyz_energy_distribution_ =
        new TH3F("xyz_energy_distribution", "xyz_energy_distribution", 100, -12., 12, 100, -12, 12, 100, -12, 12);

    config_.setDefault<bool>("store_particles", false);
    store_particles_ = config_.get<bool>("store_particles");
    neutrons_ = 0;
    (void) neutrons_;
    simple_tree_ = new TTree("neutrons", "neutrons");
    simple_tree_->Branch("initial_energy", &initial_energy_);
    simple_tree_->Branch("final_energy", &final_energy_);
    simple_tree_->Branch("particle_id", &particle_id_);
    simple_tree_->Branch("start_position_x", &start_position_x_);
    simple_tree_->Branch("start_position_y", &start_position_y_);
    simple_tree_->Branch("start_position_z", &start_position_z_);
    simple_tree_->Branch("end_position_x", &end_position_x_);
    simple_tree_->Branch("end_position_y", &end_position_y_);
    simple_tree_->Branch("end_position_z", &end_position_z_);
    simple_tree_->Branch("initial_momentum_x", &initial_momentum_x_);
    simple_tree_->Branch("initial_momentum_y", &initial_momentum_y_);
    simple_tree_->Branch("initial_momentum_z", &initial_momentum_z_);
    simple_tree_->Branch("final_momentum_x", &final_momentum_x_);
    simple_tree_->Branch("final_momentum_y", &final_momentum_y_);
    simple_tree_->Branch("final_momentum_z", &final_momentum_z_);    
    //simple_tree_->Branch("#neutrons", &neutrons_);
}

void ParticleDistributionsModule::run(unsigned int) {

    // Make a store for desired MC tracks
    std::vector<MCTrack> saved_tracks;
    std::vector<MCTrack> neutron_track;
    std::vector<MCTrack> alpha_track;
    std::vector<MCTrack> alpha_track2;
    std::vector<MCTrack> proton_track;
    std::vector<MCTrack> deuteron_track;
    std::vector<MCTrack> random;
    for(auto& particle : message_->getData()) {
        if(particle.getParticleID() == 2112) {
            neutron_track.insert(neutron_track.end() , particle); 
	}
        else if(particle.getParticleID() == 2212) {
            proton_track.insert(proton_track.end() , particle);
	}
        else if(particle.getParticleID() == 1000020040) {
            alpha_track.insert(alpha_track.end() , particle);
	}
	else if(particle.getParticleID() == 1000060120){
            deuteron_track.insert(deuteron_track.end() , particle);
	}
	else if(particle.getParticleID() > 30){
            random.insert(random.end() , particle);
	}
    }
    if(neutron_track.size() != 1) {
            alpha_track.clear();
	}
    
    /*else if(proton_track.size() == 1) {
            neutron_track.clear();
	}*/
    else if(alpha_track.size() != 2) {
            alpha_track.clear();
	}
    

    for(auto& neutron : neutron_track) 
	{
	
	//else {
	//    std::cout << neutron_track.size() <<std::endl;      	
    	//}
        ROOT::Math::XYZVector initial_momentum = neutron.getInitialMomentum();
        ROOT::Math::XYZVector final_momentum = neutron.getFinalMomentum();
        double magnitude = sqrt(initial_momentum.Mag2());
        //double magnitude_final = sqrt(final_momentum.Mag2());
        double energy = neutron.getKineticEnergyInitial();

        ROOT::Math::XYZVector directionVector, energyWeightedDirection;
        directionVector.SetX(initial_momentum.X() / magnitude);
        directionVector.SetY(initial_momentum.Y() / magnitude);
        directionVector.SetZ(initial_momentum.Z() / magnitude);

        energyWeightedDirection.SetX(energy * initial_momentum.X() / magnitude);
        energyWeightedDirection.SetY(energy * initial_momentum.Y() / magnitude);
        energyWeightedDirection.SetZ(energy * initial_momentum.Z() / magnitude);

        energy_distribution_->Fill(energy);
        zx_distribution_->Fill(directionVector.Z(), directionVector.X());
        zy_distribution_->Fill(directionVector.Z(), directionVector.Y());
        xyz_distribution_->Fill(directionVector.X(), directionVector.Y(), directionVector.Z());
        xyz_energy_distribution_->Fill(
            energyWeightedDirection.X(), energyWeightedDirection.Y(), energyWeightedDirection.Z());

        initial_energy_ = neutron.getKineticEnergyInitial();
        final_energy_ = neutron.getKineticEnergyFinal();

        particle_id_ = neutron.getParticleID();
        start_position_x_ = neutron.getStartPoint().X();
        start_position_y_ = neutron.getStartPoint().Y();
        start_position_z_ = neutron.getStartPoint().Z();
        end_position_x_ = neutron.getEndPoint().X();
        end_position_y_ = neutron.getEndPoint().Y();
        end_position_z_ = neutron.getEndPoint().Z();
        initial_momentum_x_ = initial_momentum.X();
        initial_momentum_y_ = initial_momentum.Y();
        initial_momentum_z_ = initial_momentum.Z();
        final_momentum_x_ = final_momentum.X();
        final_momentum_y_ = final_momentum.Y();
        final_momentum_z_ = final_momentum.Z();
        simple_tree_->Fill();

        if(store_particles_) {
            saved_tracks.push_back(neutron);
        }
    

    }

    // Dispatch message of pixel charges
    if(store_particles_) {
        auto mcparticle_message = std::make_shared<MCTrackMessage>(saved_tracks);
        messenger_->dispatchMessage(this, mcparticle_message);
    }

    neutron_track.clear();   
}

void ParticleDistributionsModule::finalize() {
    energy_distribution_->Write();
    zx_distribution_->Write();
    zy_distribution_->Write();
    xyz_distribution_->Write();
    xyz_energy_distribution_->Write();
    simple_tree_->Write();
}
