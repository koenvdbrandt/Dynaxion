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
}

void ParticleDistributionsModule::run(unsigned int) {

    // Make a store for desired MC tracks
    std::vector<MCTrack> saved_tracks;
    std::vector<MCTrack> neutron_track;
    std::vector<MCTrack> particle_track;

    for(auto& particle : message_->getData()) {
        particle_track.insert(particle_track.end(), particle);
        if(particle.getParticleID() == 2112)
            neutron_track.insert(neutron_track.end(), particle);
    }

    auto counter_neutron = 0;
    auto counter_proton = 0;
    auto counter_alpha = 0;
    auto counter_electron = 0;
    for(auto& vent : particle_track) {
        if(vent.getParticleID() == 2112)
            counter_neutron++;
        if(vent.getParticleID() == 2212)
            counter_proton++;
        if(vent.getParticleID() == 1000020040)
            counter_alpha++;
        if(vent.getParticleID() == 11)
            counter_electron++;
    }

    // if(counter_neutron + counter_proton < 2 ) neutron_track.clear();
    // if(counter_neutron != 1 ) neutron_track.clear();
    // if(counter_alpha != 2 ) neutron_track.clear();
    // if(counter_electron == 0 ) neutron_track.clear();

    // if(particle.getParticleID() < 200 ) {
    //     neutron_track.insert(neutron_track.end(), particle);
    //}

    for(auto& neutron : neutron_track) {

        ROOT::Math::XYZVector initial_momentum = neutron.getInitialMomentum();
        ROOT::Math::XYZVector final_momentum = neutron.getFinalMomentum();
        double magnitude = sqrt(initial_momentum.Mag2());
        // double magnitude_final = sqrt(final_momentum.Mag2());
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
