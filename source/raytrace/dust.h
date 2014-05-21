///
/// @package phosim
/// @file dust.h
/// @brief dust header file
///
/// @brief Created by
/// @author James Pizagno (UW)
///
/// @brief Modified by:
/// @author John Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

class Dust {

 public:

    double *kGrid;
    double *aGrid;
    double *bGrid;
    double *wavelengthGrid;

    void setup ();
    double calzetti (double wavelength, double A_v, double R_v);
    double calzettiSetup(double l);
    double ccm (double wavelength, double A_v, double R_v);
    void ccmSetup(double l, double *a, double *b);

};
