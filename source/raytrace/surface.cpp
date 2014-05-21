///
/// @package phosim
/// @file surface.cpp
/// @brief surface class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "surface.h"

void Surface::setup (long surfaceTotal, long points) {

    radiusCurvature = new double[surfaceTotal]();
    conic = new double[surfaceTotal]();
    height = new double[surfaceTotal]();
    outerRadius = new double[surfaceTotal]();
    innerRadius = new double[surfaceTotal]();
    three = new double[surfaceTotal]();
    four = new double[surfaceTotal]();
    five = new double[surfaceTotal]();
    six = new double[surfaceTotal]();
    seven = new double[surfaceTotal]();
    eight = new double[surfaceTotal]();
    nine = new double[surfaceTotal]();
    ten = new double[surfaceTotal]();
    centerx = new double[surfaceTotal]();
    centery = new double[surfaceTotal]();
    rmax = new double[surfaceTotal]();

    surfacemed = new int[surfaceTotal]();
    surfacecoating = new int[surfaceTotal]();
    surfacetype = new int[surfaceTotal]();
    surfacepert = new int[surfaceTotal]();

    profile = new double[surfaceTotal*points]();
    radius = new double[surfaceTotal*points]();
    normal = new double[surfaceTotal*points]();

}

void Surface::asphere (long surfaceIndex, long points) {

    double radiusofcurv = -radiusCurvature[surfaceIndex];
    double third = three[surfaceIndex]*1e3;
    double fourth = four[surfaceIndex]*1e3;
    double fifth = five[surfaceIndex]*1e3;
    double sixth = six[surfaceIndex]*1e3;
    double seventh = seven[surfaceIndex]*1e3;
    double eighth = eight[surfaceIndex]*1e3;
    double ninth = nine[surfaceIndex]*1e3;
    double tenth = ten[surfaceIndex]*1e3;

    for (long i = 0; i < points; i++) {

        radius[points*surfaceIndex + i] = innerRadius[surfaceIndex] + static_cast<double>(i)*
            (outerRadius[surfaceIndex] - innerRadius[surfaceIndex])/(static_cast<double>(points) - 1);

        if (radiusofcurv != 0) {

            profile[points*surfaceIndex + i] = height[surfaceIndex] -
                (pow(radius[points*surfaceIndex + i], 2.0)/radiusofcurv/
                (1.0 + sqrt(1.0 - (conic[surfaceIndex] + 1.0)*pow(radius[points*surfaceIndex + i]/radiusofcurv, 2.0))) +
                 third*pow(radius[points*surfaceIndex + i], 3.0) +
                 fourth*pow(radius[points*surfaceIndex + i], 4.0) +
                 fifth*pow(radius[points*surfaceIndex + i], 5.0) +
                 sixth*pow(radius[points*surfaceIndex + i], 6.0) +
                 seventh*pow(radius[points*surfaceIndex + i], 7.0) +
                 eighth*pow(radius[points*surfaceIndex + i], 8.0) +
                 ninth*pow(radius[points*surfaceIndex + i], 9.0) +
                 tenth*pow(radius[points*surfaceIndex + i], 10.0));

            normal[points*surfaceIndex + i]=-(radius[points*surfaceIndex + i]/radiusofcurv/
                                            sqrt(1.0 - (conic[surfaceIndex] + 1.0)*
                                                 pow(radius[points*surfaceIndex + i]/radiusofcurv, 2.0)) +
                                            third*pow(radius[points*surfaceIndex + i], 2.0)*3.0 +
                                            fourth*pow(radius[points*surfaceIndex + i], 3.0)*4.0 +
                                            fifth*pow(radius[points*surfaceIndex + i], 4.0)*5.0 +
                                            sixth*pow(radius[points*surfaceIndex + i], 5.0)*6.0 +
                                            seventh*pow(radius[points*surfaceIndex + i], 6.0)*7.0 +
                                            eighth*pow(radius[points*surfaceIndex + i], 7.0)*8.0 +
                                            ninth*pow(radius[points*surfaceIndex + i], 8.0)*9.0 +
                                            tenth*pow(radius[points*surfaceIndex + i], 9.0)*10.0);

        } else {

            profile[points*surfaceIndex + i] = height[surfaceIndex];
            normal[points*surfaceIndex + i] = 0.0;

        }
    }

}
