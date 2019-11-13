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
            // Excess around the chip from the pixel grid
            setScintShape(config.get<std::string>("scintillator_shape", "box"));
            if(scint_shape_ == "box") {
                setScintSize(config.get<ROOT::Math::XYZVector>("scintillator_size"));
            }
            if(scint_shape_ == "cylinder") {
                setScintRadius(config.get<double>("scintillator_radius"));
                setScintHeight(config.get<double>("scintillator_height"));
                setScintSize({2 * scint_radius_, 2 * scint_radius_, scint_height_});
            }
            setScintMaterial(config.get<std::string>("scintillator_material", "vacuum"));
            setHousingShape(config.get<std::string>("housing_shape", scint_shape_));
            setHousingThickness(config.get<double>("housing_thickness"));
            setHousingReflectivity(config.get<double>("housing_reflectivity", 0));
            setHousingMaterial(config.get<std::string>("housing_material", "vacuum"));
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
        * @brief Set the thickness of the housing of the scintillator
        * @param val housing thickness
        */
        void setScintHeight(double val) { scint_height_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        double getScintHeight() const { return scint_height_; }
        /**
        * @brief Set the thickness of the housing of the scintillator
        * @param val housing thickness
        */
        void setScintRadius(double val) { scint_radius_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        double getScintRadius() const { return scint_radius_; }

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
        * @brief Set the photo mulitplier type which is used by the scintillator
        * @param val Photo multiplier type
        */
        void setHousingShape(std::string val) { housing_shape_ = std::move(val); }

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
         * @brief Set the photo mulitplier type which is used by the scintillator
         * @return Photo multiplier type
         */
        std::string getHousingShape() const { return housing_shape_; }

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
        double scint_radius_;
        double scint_height_;
        std::string scint_material_;

        // Housing stuff
        std::string housing_shape_;
        double housing_thickness_{};
        double housing_reflectivity_{};
        std::string housing_material_;
    };
} // namespace allpix

#endif /* ALLPIX_SCINTILLATOR_H */
