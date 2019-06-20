/**
 * @file
 * @brief Implementation of passive material
 *
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

#include <Math/Rotation3D.h>
#include <Math/Translation3D.h>

#include "PassiveMaterial.hpp"
#include "core/module/exceptions.h"

using namespace allpix;

/**
 * @throws InvalidModuleActionException If the passive material model pointer is a null pointer
 *
 * Creates a passive material without any electric field in the sensor.
 */
PassiveMaterial::PassiveMaterial(std::string name,
                   ROOT::Math::XYZPoint position,
                   const ROOT::Math::Rotation3D& orientation)
    :name_(std::move(name)), position_(std::move(position)), orientation_(orientation) {}

std::string PassiveMaterial::getName() const {
    return name_;
}

ROOT::Math::XYZPoint PassiveMaterial::getPosition() const {
    return position_;
}
ROOT::Math::Rotation3D PassiveMaterial::getOrientation() const {
    return orientation_;
}

