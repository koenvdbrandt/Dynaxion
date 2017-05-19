#include "ElectricFieldInputLinearModule.hpp"

#include <fstream>
#include <limits>
#include <memory>
#include <new>
#include <stdexcept>
#include <string>
#include <utility>

#include <Math/Vector3D.h>
#include <TFile.h>
#include <TH2D.h>

#include "core/config/exceptions.h"
#include "core/geometry/PixelDetectorModel.hpp"
#include "core/utils/log.h"
#include "core/utils/unit.h"

// use the allpix namespace within this file
using namespace allpix;

// constructor to load the module
ElectricFieldInputLinearModule::ElectricFieldInputLinearModule(Configuration config,
                                                               Messenger*,
                                                               std::shared_ptr<Detector> detector)
    : Module(detector), config_(std::move(config)), detector_(std::move(detector)) {}

// init method that reads the electric field from the file
void ElectricFieldInputLinearModule::init() {
    // get field
    LOG(INFO) << "Setting electric field from bias voltage";

    // compute the electric field
    LOG(DEBUG) << getDetector()->getModel()->getSensorSize().z();
    auto field_z = config_.get<double>("voltage") / getDetector()->getModel()->getSensorSize().z();

    auto field = std::make_shared<std::vector<double>>(3);
    (*field)[0] = 0;
    (*field)[1] = 0;
    (*field)[2] = -field_z;

    // set detector field
    detector_->setElectricField(field, {{1, 1, 1}});
}
