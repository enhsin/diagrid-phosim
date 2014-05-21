///
/// @package phosim
/// @file screen.h
/// @brief header file for screen class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

class Screen {

 public:
    double large_sizeperpixel;
    double coarse_sizeperpixel;
    double medium_sizeperpixel;
    double fine_sizeperpixel;
    double *hffunc;
    double *hffunc_n;
    float *seex_coarse,*seey_coarse;
    float *seex_large, *seey_large;
    float *seex_medium, *seey_medium;
    float *phase_large, *phase_coarse, *phase_medium, *phase_fine;
    float *phaseh_medium, *phaseh_fine;
    float *cloud[MAX_LAYER];
    float *see_norm, *phase_norm;
    float secondKickSize;
    double *phasescreen;
    double *focalscreen;
    double *tfocalscreen;
    fftw_complex *inscreen;
    fftw_complex *outscreen;
    double wavelengthfactor_nom;
    double *jitterwind;

};
