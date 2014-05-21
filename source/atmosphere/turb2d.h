///
/// @package phosim
/// @file turb2d.h
/// @brief header file for turb2d
///
/// @brief Created by
/// @author J. Chiang
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef atmosphere_turb2d_h
#define atmosphere_turb2d_h

#include <string>

namespace atmosphere {

void turb2d(long seed, double see5, double outerx, double outers,
            double zenith, double wavelength, const std::string & name,
            long N_size=1024);

} // namespace atmosphere

#endif // atmosphere_turb2d_h
