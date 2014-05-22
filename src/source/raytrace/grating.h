/// @brief Grating Class
///
/// @brief Created by:
/// @author Glenn Sembroski (Purdue)
///
/// @brief Modified by:
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef GRATING_H
#define GRATING_H

#include <stdlib.h>
#include <math.h>
#include <vector>
#include "constants.h"
#include <iostream>

//const int    kMaxNumIntervals = 5400;
const int    kMaxNumIntervals = 54000; // Number of intervals to approximate
// the integration i.e. MaxNumIntervals
// = 1800 (accurate to 0.1 degree)
//const double kStepRad       = (double)M_PI/(double)kMaxNumIntervals;
//const double kStepRad       = (double)M_PI/(double)kMaxNumIntervals/2.0;
const double kStepRad       = (double)M_PI/(double)kMaxNumIntervals;
const double kBlazeAngleDeg = 7.09;
//const double kBlazeAngleDeg = -6.41;
//const double kBlazeAngleDeg = 17.45;
//const int    kNumSlits      = 1000;
const int kNumSlits = 300;
//const int kNumSlits = 1200;
//const double kDeltaNM       = 3000;   //Spacing of grating lines?
//const double kDeltaNM = 833.3;
const double kDeltaNM = 3333.33;


class Grating
// ***********************************************************************
// Class to determine diffraction from a diffraction grating
// ***********************************************************************
{
public:
    Grating();
    Grating(double BlazeAngleDeg, int NumSlits, double DeltaNM);
    ~Grating();

    void diffract(double  vxIn,  double  vyIn,  double  vzIn,
                  double vxGratingNormal, double vyGratingNormal,
                  double vzGratingNormal,
                  double& vxOut, double& vyOut, double& vzOut,
                  double wavelengthNM);

    void setAngleBlazeRad(double angRadians){fAngleBlazeRad=angRadians;return;};
    void setNumSlits(int N){fNumSlits=N;return;};
    void setDeltaNM(double deltaNM){fDeltaNM=deltaNM;return;};

private:
    double fAngleBlazeRad;
    double fDeltaNM;
    int    fNumSlits;
    double fVxGratingNormal;
    double fVyGratingNormal;
    double fVzGratingNormal;
    std::vector < double > fIntegral;

    //Methods
    double calculateFunction(double angleInRad, double angleOutRad,
                             double wavelengthNM);
    void setGratingNormal(double vxGratingNormal, double vyGratingNormal,
                          double vzGratingNormal)
    {fVxGratingNormal= vxGratingNormal; fVyGratingNormal= vyGratingNormal;
        fVzGratingNormal= vzGratingNormal; return;};

    void   makeTable(double angleInRad,double wavelengthNM);
    int    binarySearch(double goal);
    double calculateAngle(double angleInRad,double wavelengthNM);
};

#endif
