#include "RCEWriterModule.hpp"

#include <string>
#include <utility>

#include <TBranchElement.h>
#include <TClass.h>
#include <TDirectory.h>

#include "core/utils/log.h"
#include "core/utils/type.h"

#include "objects/Object.hpp"
#include "objects/objects.h"

using namespace allpix;

RCEWriterModule::RCEWriterModule(Configuration config, Messenger* messenger, GeometryManager* geo_mgr)
    : Module(config), config_(std::move(config)), geo_mgr_(geo_mgr) {
    // Bind to all messages
    //  messenger->bindMulti(this, &RCEWriterModule::pixel_charge_messages_, MsgFlags::REQUIRED);
    messenger->bindMulti(this, &RCEWriterModule::pixel_hit_messages_);
}
RCEWriterModule::~RCEWriterModule() = default;

void RCEWriterModule::init() {
    // Create output file
    std::string file_name = getOutputPath(config_.get<std::string>("file_name", "data") + ".root", true);
    output_file_ = std::make_unique<TFile>(file_name.c_str(), "RECREATE");
    output_file_->cd();

    // Initialize the events tree
    event_tree_ = std::make_unique<TTree>("Event", "");
    event_tree_->Branch("TimeStamp", &timestamp_);
    event_tree_->Branch("FrameNumber", &frame_number_);
    event_tree_->Branch("TriggerOffset", &trigger_offset_);
    event_tree_->Branch("TriggerInfo", &trigger_info_);
    event_tree_->Branch("TriggerTime", &trigger_time_);
    event_tree_->Branch("Invalid", &invalid_);

    // Get the detector names
    for(const auto& detector : geo_mgr_->getDetectors()) {
        detector_names_.push_back(detector->getName());
    }
    // Sort the detector names
    std::sort(detector_names_.begin(), detector_names_.end());
    // For each detector name, initialze an instance of SensorData
    for(const auto& detector_name : detector_names_) {
        auto& sensor = sensors_[detector_name];
        // Create directories for each detector
        TDirectory* detector = output_file_->mkdir(detector_name.c_str());
        detector->cd();

        // Initialize the struct for each detector
        sensor.tree = std::make_unique<TTree>("Hits", "");
        LOG(TRACE) << "Detector name is: " << detector_name;
        // initialze tree branches for each instance of the sensorData
        sensor.tree->Branch("NHits", &sensor.nhits_);
        sensor.tree->Branch("PixX", &sensor.pix_x_, "PixX[NHits]/I");
        sensor.tree->Branch("PixY", &sensor.pix_y_, "PixY[NHits]/I");
        sensor.tree->Branch("Value", &sensor.value_, "Value[NHits]/I");
        sensor.tree->Branch("Timing", &sensor.timing_, "Timing[NHits]/I");
        sensor.tree->Branch("HitInCluster", &sensor.hit_in_cluster_, "HitInCluster[NHits]/I");
    }
}

void RCEWriterModule::run(unsigned int event_id) {
    // fill per-event data
    timestamp_ = 0;
    frame_number_ = event_id;
    trigger_offset_ = 0;
    trigger_info_ = 0;
    trigger_time_ = 0;
    invalid_ = false;

    LOG(TRACE) << "Writing new objects to the Events tree";
    event_tree_->Fill();

    // Loop over all the detectors
    for(const auto& detector_name : detector_names_) {
        // reset nhits
        auto& sensor = sensors_[detector_name];
        sensor.nhits_ = 0;
    }

    // Loop over the pixel hit messages
    for(const auto& hit_msg : pixel_hit_messages_) {

        std::string detector_name = hit_msg->getDetector()->getName();
        LOG(TRACE) << "Detector Name: " << detector_name;
        auto& sensor = sensors_[detector_name];

        // Loop over all the hits
        for(const auto& hit : hit_msg->getData()) {
            int i = sensor.nhits_;

            // Get the pixel charge
            sensor.nhits_ += 1;
            sensor.pix_x_[i] = hit.getPixel().x();                  // NOLINT
            sensor.pix_y_[i] = hit.getPixel().y();                  // NOLINT
            sensor.value_[i] = static_cast<Int_t>(hit.getSignal()); // NOLINT
            // Set the  Timing and HitInCluster for each sesnor_tree (= 0 for now)
            sensor.timing_[i] = 0;         // NOLINT
            sensor.hit_in_cluster_[i] = 0; // NOLINT

            LOG(TRACE) << "X: " << hit.getPixel().x();
            LOG(TRACE) << "Y: " << hit.getPixel().y();
            LOG(TRACE) << "Signal: " << hit.getSignal();
        }
    }

    // Loop over all the detectors
    for(const auto& detector_name : detector_names_) {
        LOG(TRACE) << "Writing new objects to the Sensor Tree for " << detector_name;
        sensors_[detector_name].tree->Fill();
    }
}

void RCEWriterModule::finalize() {
    LOG(TRACE) << "Writing objects to file";
    output_file_->Write();
}
