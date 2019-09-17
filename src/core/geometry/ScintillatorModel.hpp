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
            setScintSize(config.get<ROOT::Math::XYZVector>("scintillator_size", {0., 0., 0.}));
            setScintMaterial(config.get<std::string>("scintillator_material", "vacuum"));
            setHousingThickness(config.get<double>("housing_thickness", 0));
            setHousingReflectivity(config.get<double>("housing_reflectivity", 0));
            setHousingMaterial(config.get<std::string>("housing_material", "vacuum"));
            setPMType(config.get<std::string>("PM_type", ""));
            if(PM_type_ == "PMT") {
                setPMTOuterRadius(config.get<double>("PMT_Outer_Radius", 0));
                setPMTHeight(config.get<double>("PMT_Height", 0));
                setPMTMaterial(config.get<std::string>("PMT_Material", "glass"));
                setPhotoCathMaterial(config.get<std::string>("PhotoCath_Material", "aluminum"));
                setNx(config.get<int>("Nx", 0));
                setNy(config.get<int>("Ny", 0));
                setNz(config.get<int>("Nz", 0));
            }
        }

        /**
         * @brief Set the size of the scintillator
         * @param val Size of scintillator
         */
        void setScintSize(ROOT::Math::XYZVector val) { scint_size_ = val; }
        /**
         * @brief Set the material of the scintillator
         * @param val Material of scintillator
         */
        void setScintMaterial(std::string val) { scint_material_ = val; }
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
        void setHousingMaterial(std::string val) { housing_material_ = val; }
        /**
        * @brief Set the photo mulitplier type which is used by the scintillator
        * @param val Photo multiplier type
        */
        void setPMType(std::string val) { PM_type_ = val; }

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
        std::string getPMType() const { return PM_type_; }

        /**
         * @brief Set outer radius of the photo multiplier tube
         * @param val Photo multiplier tube outer radius
         */
        void setPMTOuterRadius(double val) { PMT_outer_radius_ = val; }
        /**
         * @brief Get outer radius of the photo multiplier tube
         * @return Photo multiplier tube outer radius
         */
        double getPMTOuterRadius() const { return PMT_outer_radius_; }
        /**
         * @brief Set height of the photo multiplier tube
         * @param val photo multiplier tube height
         */
        void setPMTHeight(double val) { PMT_height_ = val; }
        /**
         * @brief Get height of the photo multiplier tube
         * @return Photo multiplier tube height
         */
        double getPMTHeight() const { return PMT_height_; }
        /**
         * @brief Set material of the photo multiplier tube
         * @param val Photo multiplier tube material
         */
        void setPMTMaterial(std::string val) { PMT_material_ = val; }
        /**
         * @brief Get material of the photo multiplier tube
         * @return Photo multiplier tube material
         */
        std::string getPMTMaterial() const { return PMT_material_; }
        /**
         * @brief Set the material of the photo cathode of the photo multiplier tube
         * @param val Photo cathode material
         */
        void setPhotoCathMaterial(std::string val) { photo_cath_material_ = val; }
        /**
         * @brief Get the material of the photo cathode of the photo multiplier tube
         * @return Photo cathode material
         */
        std::string getPhotoCathMaterial() const { return photo_cath_material_; }
        /**
         * @brief Set the number of PMT's in the x-direction
         * @param val Nx
         */
        void setNx(int val) { Nx_ = val; }
        /**
         * @brief Set the number of PMT's in the x-direction
         * @return Nx
         */
        int getNx() const { return Nx_; }
        /**
         * @brief Set the number of PMT's in the y-direction
         * @param Ny
         */
        void setNy(int val) { Ny_ = val; }
        /**
         * @brief Get the number of PMT's in the y-direction
         * @return Ny
         */
        int getNy() const { return Ny_; }
        /**
         * @brief Set the number of PMT's in the z-direction
         * @param Nz
         */
        void setNz(int val) { Nz_ = val; }
        /**
         * @brief Get the number of PMT's in the z-direction
         * @return Nz
         */
        int getNz() const { return Nz_; }

    private:
        ROOT::Math::XYZVector scint_size_;
        std::string scint_material_;
        double housing_thickness_{};
        double housing_reflectivity_{};
        std::string housing_material_;
        std::string PM_type_;
        double PMT_outer_radius_{};
        double PMT_height_{};
        std::string PMT_material_;
        std::string photo_cath_material_;
        int Nx_{};
        int Ny_{};
        int Nz_{};
    };
} // namespace allpix

#endif /* ALLPIX_SCINTILLATOR_H */
