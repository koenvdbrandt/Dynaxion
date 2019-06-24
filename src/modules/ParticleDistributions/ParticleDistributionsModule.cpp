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
    simple_tree_ = new TTree("protons", "protons");
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
}

void ParticleDistributionsModule::run(unsigned int) {

    // Make a store for desired MC tracks
    std::vector<MCTrack> saved_tracks;
    std::vector<MCTrack> proton_track;
    for(auto& particle : message_->getData()) {
        if(particle.getParticleID() == 2212) {
            proton_track.insert(proton_track.end() , particle);
	}
    }
       

    for(auto& proton : proton_track) 
	{
	
        ROOT::Math::XYZVector initial_momentum = proton.getInitialMomentum();
        ROOT::Math::XYZVector final_momentum = proton.getFinalMomentum();
        double magnitude = sqrt(initial_momentum.Mag2());
        //double magnitude_final = sqrt(final_momentum.Mag2());
        double energy = proton.getKineticEnergyInitial();

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

        initial_energy_ = proton.getKineticEnergyInitial();
        final_energy_ = proton.getKineticEnergyFinal();

        particle_id_ = proton.getParticleID();
        start_position_x_ = proton.getStartPoint().X();
        start_position_y_ = proton.getStartPoint().Y();
        start_position_z_ = proton.getStartPoint().Z();
        end_position_x_ = proton.getEndPoint().X();
        end_position_y_ = proton.getEndPoint().Y();
        end_position_z_ = proton.getEndPoint().Z();
        initial_momentum_x_ = initial_momentum.X();
        initial_momentum_y_ = initial_momentum.Y();
        initial_momentum_z_ = initial_momentum.Z();
        final_momentum_x_ = final_momentum.X();
        final_momentum_y_ = final_momentum.Y();
        final_momentum_z_ = final_momentum.Z();
        simple_tree_->Fill();

        if(store_particles_) {
            saved_tracks.push_back(proton);
        }
    

    }

    // Dispatch message of pixel charges
    if(store_particles_) {
        auto mcparticle_message = std::make_shared<MCTrackMessage>(saved_tracks);
        messenger_->dispatchMessage(this, mcparticle_message);
    }

    //proton_track.clear();   
}

void ParticleDistributionsModule::finalize() {
    energy_distribution_->Write();
    zx_distribution_->Write();
    zy_distribution_->Write();
    xyz_distribution_->Write();
    xyz_energy_distribution_->Write();
    simple_tree_->Write();
}
