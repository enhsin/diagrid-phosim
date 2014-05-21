///
/// @package phosim
/// @file silicon.h
/// @brief silicon header file
///
/// @brief Created by:
/// @author Andy Rasmussen (SLAC)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///


class Silicon {

 public:

    double *meanFreePath;
    float *sigma, *sigmaX, *sigmaY, *fsigma, *gsigma, *hsigma;
    float *gammaX, *gammaY;
    float *deltaX, *deltaY;
    float *nbulkmap, *deadLayer;
    double *indexRefraction;
    double *temperatureGrid;
    double *rho;
    double *wavelengthGrid;
    double *thicknessGrid;
    double *dopantGrid;
    double stopMomentPerPixel;
    double channelDepth;
    long numWavelength;
    long numTemperature;
    long numThickness;
    long numDopant;

    double abs_coeff(double lambda, double T);
    double siIndexRefraction(double lambda);
    void setup(double ccdtemp, double N_bulk, double N_f, double N_b, double s_f, double s_b, double t_si, double overdep_bias, std::string instrdir, long nampx, long nampy, double pixsize, long seedchip);
    double dope_profile(double z,double N_bulk, double N_f, double N_b, double s_f, double s_b, double t_si);
    double mu_Si (double E, double T, int polarity);


};
