///
/// @package phosim
/// @file contamination.h
/// @brief header file for contamination class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <vector>

class Contamination {

 public:

    float* transmission;
    long* surfacelistmap;
    std::vector<std::vector<double> > surfacelistx;
    std::vector<std::vector<double> > surfacelisty;
    std::vector<std::vector<double> > surfacelists;
    std::vector<std::vector<int> > surfacelistt;
    long* surfacenumber;
    float* chiptransmission;
    std::vector<double> chiplistx;
    std::vector<double> chiplisty;
    std::vector<double> chiplists;
    long* chiplistmap;
    long chipnumber;
    double *henyey_greenstein;
    double *henyey_greenstein_mu;
    double *henyey_greenstein_w;
    double *henyey_greenstein_mu_w;
    long elements;


    void setup(double *innerRadius, double *outerRadius, long totalSurface, long points, double pixsize, long nx, long ny, float qevariation, long seedchip);

};
