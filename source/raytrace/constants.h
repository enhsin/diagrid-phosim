///
/// @package phosim
/// @file vis_bi_ccd.cpp
/// @brief silicon raytrace
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

// MATH

#ifdef M_PI
    #undef M_PI
#endif
#define M_PI (3.141592653589793238462643)

#define ARCSEC 0.000004841369
#define DEGREE 0.017453293

#define HALFSQ5M1 0.618034
#define HALF3MSQ5 0.381966

// PHYSICS

#define EPSILON_0 8.85419e-14
#define E_CHARGE 1.6022e-19      // Coloumbs
#define K_BOLTZMANN 1.38066e-23  // Joules K^-1
#define EPSILON_SI 11.7     // Relative permittivity for Si
#define E_RADIUS 2.81794e-13     // cm

// OTHER

#define RADIUS_EARTH 6371.0 // km
