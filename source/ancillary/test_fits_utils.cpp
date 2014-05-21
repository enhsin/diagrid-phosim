///
/// @package phosim
/// @file test_fits_utils.cpp
/// @brief Unit tests for fits_utils classes.
///
/// @brief Created by:
/// @author J. Chiang
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <cassert>
#include <cmath>

#include <sstream>
#include <iostream>

#include "raytrace/basic_types.h"
#include "raytrace/rng_mwc.h"

#include "ancillary/fits_utils.h"

#define ASSERT_EQUALS(X, Y) assert(fabs( (X - Y)/Y ) < 1e-5)

void test_FitsImage() {
    long nx(100);
    long ny(100);
    std::vector<float> output_image;
    for (long i(0); i < nx; i++) {
        for (long j(0); j < nx; j++) {
            output_image.push_back(RandomNumbers::RngFloat());
        }
    }
    std::string image_file("test_image.fits");
    fits_utils::FitsImage<float> output(output_image, nx, ny, "!" + image_file);
    output.write();
    output.close();

    long naxes[2];
    std::vector<float> input_image;
    fits_utils::FitsImage<float> input(input_image, naxes, image_file);

    assert(naxes[0] == nx);
    assert(naxes[1] == ny);
    for (long i(0); i < nx*ny; i++) {
        assert(output_image[i] == input_image[i]);
    }
}

int main(int iargc, char * argv[]) {
    test_FitsImage();
    return 0;
}
