/**
 * @file
 * @brief Parameters of a hybrid pixel detector model
 *
 * @copyright Copyright (c) 2017-2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_SCINTILLATOR_H
#define ALLPIX_SCINTILLATOR_H

#include <string>
#include <utility>

#include <Math/Cartesian2D.h>
#include <Math/DisplacementVector2D.h>
#include <Math/Point3D.h>
#include <Math/Vector2D.h>
#include <Math/Vector3D.h>

#include "DetectorModel.hpp"

namespace allpix {

    /**
     * @ingroup DetectorModels
     * @brief Model of a scintillator. This a model where a scintillating material creates optical photons that will get
     * measured by the sensor
     */
    class ScintillatorModel : public DetectorModel {
    public:
        explicit ScintillatorModel(std::string type, const ConfigReader& reader) : DetectorModel(std::move(type), reader) {
            auto config = reader.getHeaderConfiguration();
            setScintShape(config.get<std::string>("scintillator_shape", "box"));
            if(scint_shape_ == "box") {
                setScintSize(config.get<ROOT::Math::XYZVector>("scintillator_size"));
            }
            if(scint_shape_ == "cylinder") {
                auto scint_size = config.get<ROOT::Math::XYVector>("scintillator_size");
                setScintSize({scint_size.x(), scint_size.x(), scint_size.y()});
            }
            setScintMaterial(config.get<std::string>("scintillator_material", "vacuum"));
            setHousingThickness(config.get<double>("housing_thickness"));
            setHousingReflectivity(config.get<double>("housing_reflectivity", 0));
            setHousingMaterial(config.get<std::string>("housing_material", "vacuum"));
            // Set defaults for unused DetectorModel values
            setNPixels({1, 1});
            setPixelSize({sensor_size_.x(), sensor_size_.y()});
            setSensorThickness(sensor_size_.z());
        }

        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setScintShape(std::string val) { scint_shape_ = std::move(val); }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        std::string getScintShape() const { return scint_shape_; }
        /**
         * @brief Set the size of the scintillator
         * @param val Size of scintillator
         */
        void setScintSize(ROOT::Math::XYZVector val) { scint_size_ = std::move(val); }
        /**
         * @brief Set the material of the scintillator
         * @param val Material of scintillator
         */
        void setScintMaterial(std::string val) { scint_material_ = std::move(val); }
        /**
         * @brief the size of the scintillator
         * @return Size of scintillator
         */
        ROOT::Math::XYZVector getScintSize() const { return scint_size_; }
        /**
         * @brief Get the material of the scintillator
         * @return Material of scintillator
         */
        std::string getScintMaterial() const { return scint_material_; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val Housing thickness
         */
        void setHousingThickness(double val) { housing_thickness_ = val; }
        /**
         * @brief Set the reflectivity of the housing of the scintillator
         * @param val Housing reflectivity
         */
        void setHousingReflectivity(double val) { housing_reflectivity_ = val; }
        /**
         * @brief Set the material of the housing of the scintillator
         * @param val Hhousing material
         */
        void setHousingMaterial(std::string val) { housing_material_ = std::move(val); }

        /**
         * @brief Get the thickness of the housing of the scintillator
         * @return Housing thickness
         */
        double getHousingThickness() const { return housing_thickness_; }
        /**
         * @brief Set the reflectivity of the housing of the scintillator
         * @return Housing reflectivity
         */
        double getHousingReflectivity() const { return housing_reflectivity_; }
        /**
         * @brief Get the material of the housing of the scintillator
         * @return Housing material
         */
        std::string getHousingMaterial() const { return housing_material_; }

        /**
         * @brief Get size of the Detector
         * @return Size of the Detector, which is the housing size
         */
        ROOT::Math::XYZVector getSize() const override {
            return ROOT::Math::XYZVector(getScintSize().x() + 2 * getHousingThickness(),
                                         getScintSize().y() + 2 * getHousingThickness(),
                                         getScintSize().z() + getSensorSize().z() + 2 * getHousingThickness());
        }
        /**
         * @brief Get local coordinate of the geometric center of the model
         * @note This returns the center of the geometry model, which is the housing center
         */
        ROOT::Math::XYZPoint getGeometricalCenter() const override {
            return ROOT::Math::XYZPoint(
                getCenter().x(), getCenter().y(), getSize().z() / 2 - getHousingThickness() - getSensorSize().z() / 2);
        }

    private:
        // Scintillator stuff
        std::string scint_shape_;
        ROOT::Math::XYZVector scint_size_;
        std::string scint_material_;

        // Housing stuff
        double housing_thickness_{};
        double housing_reflectivity_{};
        std::string housing_material_;
    };
} // namespace allpix

#endif /* ALLPIX_SCINTILLATOR_H */
