///
/// @package phosim
/// @file surface.h
/// @brief header file for surface class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef Surface_H
#define Surface_H

#include <math.h>

enum SurfaceTypes {MIRROR=0,LENS=1,FILTER=2,DETECTOR=3,GRATING=4};

class Surface {

 public:
    double* radiusCurvature;
    double* conic;
    double* height;
    double* outerRadius;
    double* innerRadius;
    double* three;
    double* four;
    double* five;
    double* six;
    double* seven;
    double* eight;
    double* nine;
    double* ten;
    double* centerx;
    double* centery;
    double* rmax;

    double* normal;
    double* radius;
    double* profile;

    int* surfacecoating;
    int* surfacetype;
    int* surfacepert;
    int* surfacemed;

    void setup(long surfaceTotal, long points);
    void asphere(long surfaceIndex, long points);

};

#endif
