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
         * @return Size of the Detector
         *
         * Calculated from \ref DetectorModel::getGridSize "pixel grid size", chip excess and chip thickness
         */
        ROOT::Math::XYZVector getSize() const override {
            ROOT::Math::XYZVector max(std::numeric_limits<double>::lowest(),
                                      std::numeric_limits<double>::lowest(),
                                      std::numeric_limits<double>::lowest());
            ROOT::Math::XYZVector min(
                std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());

            ROOT::Math::XYZPoint HousingCenter = {getSensorCenter().x(),
                                                  getSensorCenter().y(),
                                                  getSensorCenter().z() + getSensorSize().z() / 2 + getScintSize().z() / 2};
            ROOT::Math::XYZVector HousingSize = {getScintSize().x(), getScintSize().y(), getScintSize().z()};

            std::array<ROOT::Math::XYZPoint, 2> centers = {{getSensorCenter(), HousingCenter}};
            std::array<ROOT::Math::XYZVector, 2> sizes = {{getSensorSize(), HousingSize}};

            for(size_t i = 0; i < 2; ++i) {
                max.SetX(std::max(max.x(), (centers.at(i) + sizes.at(i) / 2.0).x() + getHousingThickness()));
                max.SetY(std::max(max.y(), (centers.at(i) + sizes.at(i) / 2.0).y() + getHousingThickness()));
                max.SetZ(std::max(max.z(), (centers.at(i) + sizes.at(i) / 2.0).z() + getHousingThickness()));
                min.SetX(std::min(min.x(), (centers.at(i) - sizes.at(i) / 2.0).x() + getHousingThickness()));
                min.SetY(std::min(min.y(), (centers.at(i) - sizes.at(i) / 2.0).y() + getHousingThickness()));
                min.SetZ(std::min(min.z(), (centers.at(i) - sizes.at(i) / 2.0).z() + getHousingThickness()));
            }

            for(auto& support_layer : getSupportLayers()) {
                auto size = support_layer.getSize();
                auto center = support_layer.getCenter();
                max.SetX(std::max(max.x(), (center + size / 2.0).x()));
                max.SetY(std::max(max.y(), (center + size / 2.0).y()));
                max.SetZ(std::max(max.z(), (center + size / 2.0).z()));
                min.SetX(std::min(min.x(), (center - size / 2.0).x()));
                min.SetY(std::min(min.y(), (center - size / 2.0).y()));
                min.SetZ(std::min(min.z(), (center - size / 2.0).z()));
            }

            ROOT::Math::XYZVector size;
            size.SetX(2 * std::max(max.x() - getCenter().x(), getCenter().x() - min.x()));
            size.SetY(2 * std::max(max.y() - getCenter().y(), getCenter().y() - min.y()));
            size.SetZ((max.z() - getCenter().z()) +
                      (getCenter().z() - min.z())); // max.z() is positive (chip side) and min.z() is negative (sensor side)

            return size;
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
