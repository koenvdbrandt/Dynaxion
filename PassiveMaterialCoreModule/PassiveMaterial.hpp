/**
 * @file
 * @brief Base of passive material implementation
 *
 * @copyright Copyright (c) 2017 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_PASSIVE_MATERIAL_H
#define ALLPIX_PASSIVE_MATERIAL_H

#include <array>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <typeindex>
#include <vector>

#include <Math/Point3D.h>
#include <Math/Rotation3D.h>
#include <Math/Transform3D.h>

#include "PassiveMaterial.hpp"

#include "objects/Pixel.hpp"

namespace allpix {

    /**
     * @brief Instantiation of a passive material model in the world
     *
     * Contains the passive material in the world with several unique properties (like the electric field).
     */
    class PassiveMaterial {
        friend class GeometryManager;

    public:
        /**
         * @brief Constructs a passive material in the geometry
         * @param name Unique name of the passive material
         * @param position Position in the world frame
         * @param orientation Rotation matrix representing the orientation
         */
        PassiveMaterial(std::string name,
                 ROOT::Math::XYZPoint position,
                 const ROOT::Math::Rotation3D& orientation);

        /**
         * @brief Get name of the passive material
         * @return PassiveMaterial name
         */
        std::string getName() const;

        /**
         * @brief Get position in the world
         * @return Global position in Cartesian coordinates
         */
        ROOT::Math::XYZPoint getPosition() const;
        /**
         * @brief Get orientation in the world
         * @return Rotation matrix representing the orientation
         */
        ROOT::Math::Rotation3D getOrientation() const;

        /**
         * @brief Fetch an external object linked to this passive material
         * @param name Name of the external object
         * @return External object or null pointer if it does not exists
         */
        template <typename T> std::shared_ptr<T> getExternalObject(const std::string& name);
        /**
         * @brief Sets an external object linked to this passive material
         * @param name Name of the external object
         * @param model External object of arbitrary type
         */
        template <typename T> void setExternalObject(const std::string& name, std::shared_ptr<T> model);

    private:

        std::string name_;

        ROOT::Math::XYZPoint position_;
        ROOT::Math::Rotation3D orientation_;

        std::map<std::type_index, std::map<std::string, std::shared_ptr<void>>> external_objects_;
    };

    /**
     * If the returned object is not a null pointer it is guaranteed to be of the correct type
     */
    template <typename T> std::shared_ptr<T> PassiveMaterial::getExternalObject(const std::string& name) {
        return std::static_pointer_cast<T>(external_objects_[typeid(T)][name]);
    }
    /**
     * Stores external representations of objects in this passive material that need to be shared between modules.
     */
    template <typename T> void PassiveMaterial::setExternalObject(const std::string& name, std::shared_ptr<T> model) {
        external_objects_[typeid(T)][name] = std::static_pointer_cast<void>(model);
    }
} // namespace allpix

#endif /* ALLPIX_PASSIVE_MATERIAL_H */
