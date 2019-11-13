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
            setScintShape(config.get<std::string>("scintillator_shape", "square"));
            if(scint_shape_ == "square")
                setScintSize(config.get<ROOT::Math::XYZVector>("scintillator_size", {0., 0., 0.}));
            if(scint_shape_ == "cylinder") {
                setScintRadius(config.get<double>("scintillator_radius", 0));
                setScintHeight(config.get<double>("scintillator_height", 0));
                setScintSize({scint_radius_, scint_radius_, scint_height_});
            }
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
        void setScintShape(std::string val) { scint_shape_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        std::string getScintShape() const {
            return scint_shape_;
        } /**
* @brief Set the thickness of the housing of the scintillator
* @param val housing thickness
*/
        void setScintHeight(double val) { scint_height_ = val; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        double getScintHeight() const {
            return scint_height_;
        } /**
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
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        void setScintSize(ROOT::Math::XYZVector val) { scint_size_ = val; }
        /**
 * @brief Set the thickness of the housing of the scintillator
 * @param val housing thickness
 */
        /////    void setScintRad(double val) { scint_rad_ = val; }
        /**
 * @brief Set the thickness of the housing of the scintillator
 * @param val housing thickness
 */
        //  ///     void setScintHeight(double val) { scint_height_ = val; }
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
        ///       double getScintHeight() const { return scint_height_; }
        /**
         * @brief Set the thickness of the housing of the scintillator
         * @param val housing thickness
         */
        ///      double getScintRad() const { return scint_rad_; }
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
            size.SetZ(2 * std::max(max.z() - getCenter().z(), getCenter().z() - min.z()));

            return size;
        }

    private:
        std::string scint_shape_;
        ROOT::Math::XYZVector scint_size_;
        double scint_radius_;
        double scint_height_;

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
