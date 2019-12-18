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
            setScintMaterial(config.get<std::string>("scintillator_material"));
            setHousingThickness(config.get<double>("housing_thickness"));
            setHousingReflectivity(config.get<double>("housing_reflectivity", 0));
            setHousingMaterial(config.get<std::string>("housing_material", "vacuum"));
            // FIXME Need something to require scintillator size
            auto sensor_size = config.get<ROOT::Math::XYZVector>("sensor_size");
            (void)sensor_size;
            // Set defaults for unused DetectorModel values

            setNPixels({1, 1});
            setPixelSize({sensor_size_.x(), sensor_size_.y()});
            setSensorThickness(sensor_size_.z());
            // Set Optical Surface properties
            setHousingSurfaceModel(config.get<int>("housing_surface_model", 0));
            setHousingSurfaceType(config.get<int>("housing_surface_type", 0));
            setHousingSurfaceFinish(config.get<int>("housing_surface_finish", 0));
            setHousingSurfaceValue(config.get<double>("housing_surface_value", 1));
            setHousingSurfaceRIndex(config.get<double>("housing_surface_refractive_index", 0));
            setPhotocathodeSurfaceModel(config.get<int>("photocathode_surface_model", 0));
            setPhotocathodeSurfaceType(config.get<int>("photocathode_surface_type", 0));
            setPhotocathodeSurfaceFinish(config.get<int>("photocathode_surface_finish", 0));
            setPhotocathodeSurfaceValue(config.get<double>("photocathode_surface_value", 1));
            setPhotocathodeSurfaceRIndex(config.get<double>("photocathode_surface_refractive_index", 0));
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

        void setHousingSurfaceModel(int val) { housing_surface_model_ = val; }
        void setHousingSurfaceType(int val) { housing_surface_type_ = val; }
        void setHousingSurfaceFinish(int val) { housing_surface_finish_ = val; }
        void setHousingSurfaceValue(double val) { housing_surface_value_ = val; }
        void setHousingSurfaceRIndex(double val) { housing_surface_refractive_index_ = val; }
        void setPhotocathodeSurfaceModel(int val) { photocathode_surface_model_ = val; }
        void setPhotocathodeSurfaceType(int val) { photocathode_surface_type_ = val; }
        void setPhotocathodeSurfaceFinish(int val) { photocathode_surface_finish_ = val; }
        void setPhotocathodeSurfaceValue(double val) { photocathode_surface_value_ = val; }
        void setPhotocathodeSurfaceRIndex(double val) { photocathode_surface_refractive_index_ = val; }

        int getHousingSurfaceModel() const { return housing_surface_model_; }
        int getHousingSurfaceType() const { return housing_surface_type_; }
        int getHousingSurfaceFinish() const { return housing_surface_finish_; }
        double getHousingSurfaceValue() const { return housing_surface_value_; }
        double getHousingSurfaceRIndex() const { return housing_surface_refractive_index_; }
        int getPhotocathodeSurfaceModel() const { return photocathode_surface_model_; }
        int getPhotocathodeSurfaceType() const { return photocathode_surface_type_; }
        int getPhotocathodeSurfaceFinish() const { return photocathode_surface_finish_; }
        double getPhotocathodeSurfaceValue() const { return photocathode_surface_value_; }
        double getPhotocathodeSurfaceRIndex() const { return photocathode_surface_refractive_index_; }

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

        // Optical surface stuff
        int housing_surface_model_{};
        int housing_surface_type_{};
        int housing_surface_finish_{};
        double housing_surface_value_{};
        double housing_surface_refractive_index_{};
        int photocathode_surface_model_{};
        int photocathode_surface_type_{};
        int photocathode_surface_finish_{};
        double photocathode_surface_value_{};
        double photocathode_surface_refractive_index_{};
    };
} // namespace allpix

#endif /* ALLPIX_SCINTILLATOR_H */
