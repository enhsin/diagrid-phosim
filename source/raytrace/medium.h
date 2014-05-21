///
/// @package phosim
/// @file medium.h
/// @brief header file for medium class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

class Medium {

 public:
    long index_refraction_number[MAX_SURF];
    double *index_refraction[MAX_SURF];
    double *index_refraction_wavelength[MAX_SURF];

};
