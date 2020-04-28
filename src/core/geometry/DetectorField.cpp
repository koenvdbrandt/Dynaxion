/**
 * @file
 * @brief Implementation of detector fields
 *
 * @copyright Copyright (c) 2019-2020 CERN and the Allpix Squared authors.
 * This software is distributed under the terms of the MIT License, copied verbatim in the file "LICENSE.md".
 * In applying this license, CERN does not waive the privileges and immunities granted to it by virtue of its status as an
 * Intergovernmental Organization or submit itself to any jurisdiction.
 */

#include "DetectorField.hpp"

namespace allpix {

    /*
     * Vector field template specialization of helper function for field flipping
     */
    template <> void flip_vector_components<ROOT::Math::XYZVector>(ROOT::Math::XYZVector& vec, bool x, bool y) {
        if(x) {
            vec.SetX(-vec.x());
        }
        if(y) {
            vec.SetY(-vec.y());
        }
    }

    /*
     * Scalar field template specialization of helper function for field flipping
     * Here, no inversion of the field components is required
     */
    template <> void flip_vector_components<double>(double&, bool, bool) {}
} // namespace allpix
