/**
 * @file
 * @brief Parameters of a hybrid pixel detector model
 *
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
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
     * @brief Model of a hybrid pixel detector. This a model where the sensor is bump-bonded to the chip
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
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setScintSize(ROOT::Math::XYZVector val) { scint_size_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setScintMaterial(std::string val) { scint_material_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        ROOT::Math::XYZVector getScintSize() const { return scint_size_; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        std::string getScintMaterial() const { return scint_material_; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setHousingThickness(double val) { housing_thickness_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setHousingReflectivity(double val) { housing_reflectivity_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setHousingMaterial(std::string val) { housing_material_ = val; }
        /**
        * @brief Set the thickness of the housing of the scintillator
        * @param val housing thickness
        */
        void setPMType(std::string val) { PM_type_ = val; }

        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        double getHousingThickness() const { return housing_thickness_; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        double getHousingReflectivity() const { return housing_reflectivity_; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        std::string getHousingMaterial() const { return housing_material_; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        std::string getPMType() const { return PM_type_; }

        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setPMTOuterRadius(double val) { PMT_outer_radius_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        double getPMTOuterRadius() const { return PMT_outer_radius_; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setPMTHeight(double val) { PMT_height_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        double getPMTHeight() const { return PMT_height_; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setPMTMaterial(std::string val) { PMT_material_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        std::string getPMTMaterial() const { return PMT_material_; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setPhotoCathMaterial(std::string val) { photo_cath_material_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        std::string getPhotoCathMaterial() const { return photo_cath_material_; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setNx(int val) {
            Nx_ = val;
        } /**
* @brief Set the thickness of the housing of the scintillator
* @param val housing thickness
*/
        int getNx() const {
            return Nx_;
        } /**
* @brief Set the thickness of the housing of the scintillator
* @param val housing thickness
*/
        void setNy(int val) {
            Ny_ = val;
        } /**
* @brief Set the thickness of the housing of the scintillator
* @param val housing thickness
*/
        int getNy() const {
            return Ny_;
        } /**
* @brief Set the thickness of the housing of the scintillator
* @param val housing thickness
*/
        void setNz(int val) {
            Nz_ = val;
        } /**
* @brief Set the thickness of the housing of the scintillator
* @param val housing thickness
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
