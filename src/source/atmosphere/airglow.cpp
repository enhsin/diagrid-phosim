///
/// @package phosim
/// @file airglow.cpp
/// @brief Functions to calculate airglow screens.
///
/// @brief Created by
/// @author En-Hsin Peng (Purdue)
///
/// @brief Modified by
/// @author J. Chiang (SLAC)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>

#include <sstream>
#include <stdexcept>
#include <vector>

#include "raytrace/basic_types.h"
#include "raytrace/helpers.h"
#include "raytrace/rng_mwc.h"

#include "ancillary/fftw_utils.h"
#include "ancillary/fits_utils.h"

#include "atmosphere/airglow.h"

using namespace RandomNumbers;

namespace atmosphere {

    void airglow(long seed, const std::string & name, long screenSize) {
        long nx = screenSize;
        long ny = screenSize;
        long N = nx*ny;

        std::vector<double> in1(N, 0);
        std::vector<double> in(N, 0);

        static int ns(3);
        std::vector<double> kmax(ns, 0);
        std::vector<double> kmin(ns, 0);
        std::vector<double> pixs(ns, 0);
        pixs[0] = 15.0/3600;  //in degrees
        kmin[0] = 4.0/screenSize/pixs[0];
        kmax[0] = 1.0/pixs[0];
        for (int ii=1; ii < ns; ii++) {
            pixs[ii] = pixs[ii-1]*16;
            kmax[ii] = kmin[ii-1];
            kmin[ii] = kmin[ii-1]/16;
        }

        fftw_utils::ComplexArray out(nx, ny);
        fftw_utils::ComplexArray incom(nx, ny);

        /* initialize random seed */
        RngSetSeed32(seed);
        RngUnwind(10000);
        double idxa=1.3, idxb=3.0;
        double n0a=6.6e-3;
        double n0b=1e-5;

        for (int ii=0; ii<ns; ii++) {
            std::vector<double> xgrid(screenSize, 0);
            for (int ix=0; ix<nx; ix++) {
                xgrid[ix]=ix*pixs[ii];
                double kx=(ix<(nx/2)) ? (double)ix: (double)(nx-ix);
                for(int iy=0; iy<ny; iy++) {
                    double ky=(iy<(ny/2)) ? (double)iy: (double)(ny-iy);
                    long i = iy + ny*ix;
                    double dkk=sqrt(kx*kx+ky*ky)/screenSize/pixs[ii];
                    out[i][0]=0.0;
                    out[i][1]=0.0;
                    if (dkk>=kmax[ii] || dkk<kmin[ii]) {
                        out[i][0]=0.0;
                    } else {
                        double dkkp=(n0a*pow(dkk,-idxa)+n0b*pow(dkk,-idxb))/pixs[ii];
                        out[i][0] = dkkp*random_gaussian();
                    }
                }
            }
            out[0][0]=0.0;
            out.inverse_fft(incom);

            for (long i=0; i < N; i++) {
                in1[i] = incom[i][1]/screenSize;
            }

            double xpos = RngDouble()*(nx*pixs[ii] - nx*pixs[0]);
            double ypos = RngDouble()*(ny*pixs[ii] - ny*pixs[0]);
            for(int ix=0; ix<nx; ix++) {
                for(int iy=0; iy<ny; iy++) {
                    long i = iy + ny*ix;
                    double dx, dy;
                    long idx0 = find_linear(const_cast<double *>(&xgrid[0]), screenSize,
                                            xpos + ix*pixs[0], &dx);
                    long idy0 = find_linear(const_cast<double *>(&xgrid[0]), screenSize,
                                            ypos + iy*pixs[0], &dy);
                    in[i] += interpolate_bilinear(&in1[0], screenSize, idx0, dx, idy0, dy);
                }
            }
            //       std::ostringstream outfile;
            //       outfile << "!" << name << "_" << ii << ".fits";
            //       fits_utils::write_fits_image(in1, nx, ny, outfile.str());
        }

        double avg = 0.0;
        double rms = 0.0;
        double imax=0., imin=0.;
        for(long i=0; i<N; i++) { 
            avg += in[i];
            imax=(in[i]>imax)? in[i]:imax;
            imin=(in[i]<imin)? in[i]:imin;
        }
        avg /= N;
        for(long i=0; i<N; i++) {
            rms += (in[i]-avg)*(in[i]-avg);
        }
        rms /= N;
        rms = sqrt(rms);
        //printf("%g %g %g %g\n",avg,imax,imin,rms);

        for(long i=0; i<N; i++) {
            in[i]-=avg;
        }


        std::ostringstream outfile;
        outfile << "!" << name << ".fits.gz";
        fits_utils::write_fits_image(in, nx, ny, outfile.str());
    }


} // namespace atmosphere
