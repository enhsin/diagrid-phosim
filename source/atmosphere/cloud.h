///
/// @package phosim
/// @file cloud.h
/// @brief header file for cloud
///
/// @brief Created by
/// @author J. Chiang
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef atmosphere_cloud_h
#define atmosphere_cloud_h

#include <string>

namespace atmosphere {

   void cloud(long seed, double cloheight, double pixsz,
              const std::string & name, long N_size = 1024);

} // namespace atmosphere

#endif // atmosphere_cloud_h
