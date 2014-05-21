///
/// @package phosim
/// @file obstruction.h
/// @brief header file for obstruction class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

class Obstruction {

 public:
    double spider_par[MAX_SURF];
    double spider_center[MAX_SURF];
    double spider_width[MAX_SURF];
    double spider_depth[MAX_SURF];
    double spider_angle[MAX_SURF];
    double spider_reference[MAX_SURF];
    double spider_height[MAX_SURF];
    int spider_type[MAX_SURF];
    int nspid;
    int pupil;

};
