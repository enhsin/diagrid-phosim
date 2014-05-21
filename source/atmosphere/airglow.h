///
/// @package phosim
/// @file airglow.h
/// @brief header file for airglow
///
/// @brief Created by
/// @author J. Chiang (SLAC)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef atmosphere_airglow_h
#define atmosphere_airglow_h

#include <string>
#include <vector>

namespace atmosphere {

    void airglow(long seed, const std::string & name, long screenSize=1024);

} // namespace atmosphere

#endif // atmosphere_airglow_h
