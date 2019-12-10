/**
 * @file
 * @brief Parameters of a box passive material model
 *
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_PASSIVE_MATERIAL_BOX_H
#define ALLPIX_PASSIVE_MATERIAL_BOX_H

#include <string>
#include <utility>

#include <Math/Cartesian2D.h>
#include <Math/DisplacementVector2D.h>
#include <Math/Point3D.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>

#include <G4Box.hh>
#include <G4SubtractionSolid.hh>
#include <G4VSolid.hh>
#include "PassiveMaterialModel.hpp"

namespace allpix {

    /**
     * @ingroup PassiveMaterialModel
     * @brief Model of a rectangular box.
     */
    class BoxModel : public PassiveMaterialModel {
    public:
        /**
         * @brief Constructs the box passive material model
         * @param config Configuration with description of the model
         */
        explicit BoxModel(Configuration& config) : PassiveMaterialModel(config), config_(config) {

            // Set the box specifications
            setOuterSize(config.get<ROOT::Math::XYZVector>("size"));
            setInnerSize(config.get<ROOT::Math::XYVector>("inner_size", ROOT::Math::XYVector()));
            auto thickness = config.get<double>("thickness", 0);
            if(thickness != 0) {
                if(inner_size_ != ROOT::Math::XYVector()) {
                    throw InvalidValueError(config_, "thickness", "cannot have both 'thickness' and 'inner_size'");
                }

                setInnerSize({outer_size_.x() - thickness, outer_size_.y() - thickness});
            }

            std::string name = config_.getName();
            // Limit the values that can be given
            if(inner_size_.x() >= outer_size_.x() || inner_size_.y() >= outer_size_.y()) {
                throw InvalidValueError(config_, "inner_size_", "inner_size_ cannot be larger than the outer_size_");
            }

            // Create the G4VSolids which make the Box
            auto outer_volume =
                new G4Box(name + "_outer_volume", outer_size_.x() / 2, outer_size_.y() / 2, outer_size_.z() / 2);
            if(inner_size_ == ROOT::Math::XYVector()) {
                solid_ = outer_volume;
            } else {
                auto inner_volume =
                    new G4Box(name + "_inner_volume", inner_size_.x() / 2, inner_size_.y() / 2, 1.1 * outer_size_.z() / 2);

                solid_ = new G4SubtractionSolid(name + "_volume", outer_volume, inner_volume);
            }
            // Get the maximum of the size parameters
            max_size_ = std::max(outer_size_.x(), std::max(outer_size_.y(), outer_size_.z()));
        }
        /**
         * @brief Set the XYZ-value of the outer size of the box
         * @param val Outer size of the box
         */
        void setOuterSize(ROOT::Math::XYZVector val) { outer_size_ = std::move(val); }
        /**
         * @brief Set the XYZ-value of the outer size of the box
         * @param val Outer size of the box
         */
        void setInnerSize(ROOT::Math::XYVector val) { inner_size_ = std::move(val); }

        // Set the override functions of PassiveMaterialModel
        G4VSolid* getSolid() override { return solid_; }
        double getMaxSize() override { return max_size_; }

    private:
        Configuration& config_;

        // G4VSolid returnables
        G4VSolid* solid_;

        double max_size_;

        // G4VSolid specifications
        ROOT::Math::XYZVector outer_size_;
        ROOT::Math::XYVector inner_size_;
    };
} // namespace allpix

#endif /* ALLPIX_PASSIVE_MATERIAL_BOX_H */
