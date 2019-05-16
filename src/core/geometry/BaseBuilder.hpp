/**
 * @file
 * @brief Base of geometry building
 *
 * @copyright Copyright (c) 2019 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#ifndef ALLPIX_BASE_BUILDER_H
#define ALLPIX_BASE_BUILDER_H

#include <array>
#include <functional>
#include <map>
#include <memory>
#include <string>
#include <tuple>
#include <typeindex>
#include <utility>
#include <vector>

#include <Math/Point3D.h>
#include <Math/Rotation3D.h>
#include <Math/Transform3D.h>
#include "GeometryBuilder.hpp"

namespace allpix {
    template<typename WorldVolume, typename Materials>
    class BaseBuilder: public GeometryBuilder{

    public:
        /**
             * @brief Empty Build function with 2 inputs
             */

        virtual void build(WorldVolume* world_log, std::map<std::string, Materials*> materials_);

        /*
             * @brief Essential virtual destructor
         */    

        //virtual ~BaseBuilder() = default;

    private:
    };

} // namespace allpix

#endif /* ALLPIX_BASE_BUILDER_H */