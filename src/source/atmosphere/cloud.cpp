///
/// @package phosim
/// @file cloud.cpp
/// @brief makes cloud screens
///
/// @brief Created by
/// @author J. Garrett Jernigan (SLAC)
///
/// @brief Modified by
/// @author John R. Peterson (Purdue)
/// @author En-Hsin Peng (Purdue)
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

#include "raytrace/basic_types.h"
#include "raytrace/rng_mwc.h"
#include "ancillary/fftw_utils.h"
#include "ancillary/fits_utils.h"
#include "atmosphere/cloud.h"

using namespace RandomNumbers;

namespace atmosphere {

    void cloud(long seed, double cloheight, double pixsz, const std::string & name, long N_size) {

        /*
        // Angular correlation scale degrees.
        static double angcor = 2.0;
        // Spatial correlation scale pixels.
        double spacor = (1000.0*100.0*cloheight/pixsz)*tan(3.1415926*angcor/180.0);
        */

        // correlation scale meters.
        double phycor = 200.0;
        // Spatial correlation scale pixels
        double spacor = phycor*100.0/pixsz;

        // spatial correlation freq inverse pixels
        double freqcor = static_cast<double>(N_size)/spacor;

        long nx = N_size;
        long ny = N_size;
        int N = nx*ny;

        fftw_utils::ComplexArray out(nx, ny);
        fftw_utils::ComplexArray incom(nx, ny);

        // initialize random seed
        RngSetSeed32(seed);
        RngUnwind(10000);

        double dkksqr = 0.5*sqrt(static_cast<double>(nx*nx + ny*ny));

        // build a random cloud screen
        for (int ix=0; ix < nx; ix++) {
            double kx;
            if (ix < (nx/2)) {
                kx = static_cast<double>(ix);
            } else {
                kx = static_cast<double>(nx - ix);
            }
            for (int iy=0; iy < ny; iy++) {
                int i = iy + ny*ix;
                double ky;
                if (iy < (ny/2)) {
                    ky = static_cast<double>(iy);
                } else {
                    ky = static_cast<double>(ny - iy);
                }

                // random Gaussian
                double dkk = std::sqrt(kx*kx + ky*ky);
                double dkkp = 1.0/std::exp(dkk/freqcor);

                if (i == 0) {
                    dkkp = 0.0;
                }

                // cutout the corner of the Fourier plane - no x,y effect
                if (dkk > dkksqr) {
                    dkkp = 0.0;
                }

                out[i][0] = dkkp*random_gaussian();
                out[i][1] = dkkp*random_gaussian();
                int idkk = static_cast<int>(dkk);
                if (idkk == 0) {
                    continue;
                }
            }
        }
        // finish building random cloud screen

        // inverse FFT
        out.inverse_fft(incom);

        // normalize FFT
        std::vector<double> in(N, 0);
        double avg = 0.0;
        for (int i=0; i < N; i++) {
            in[i] = incom[i][0]/static_cast<double>(N);
            avg += in[i];
        }
        avg /= static_cast<double>(N);

        double rms = 0.0;
        for (int i=0; i < N; i++) {
            rms += (in[i] - avg)*(in[i] - avg);
        }
        rms /= static_cast<double>(N);
        rms = std::sqrt(rms);

        double xmax = in[0];
        double xmin = in[0];
        for (int i=0; i < N; i++) {
            if (in[i] > xmax) {
                xmax = in[i];
            }
            if (in[i] < xmin) {
                xmin = in[i];
            }
        }

        // output fits file cloud.fits
        // output scalar density as fits file
        // convert double results to a float

        double denavg = 0.0;
        double denrms = 0.0;
        for (int ix=0; ix < nx; ix++) {
            for (int iy=0; iy < ny; iy++) {
                denavg += in[iy+ny*ix];
                denrms += in[iy+ny*ix]*in[iy+ny*ix];
            }
        }
        denavg /= static_cast<double>(nx*ny);
        denrms = std::sqrt(denrms/static_cast<double>(nx*ny));

        std::vector<float> infl(nx*ny, 0);
        for (int ix=0; ix < nx; ix++) {
            for (int iy=0; iy < ny; iy++) {
                infl[iy+ny*ix] = in[iy+ny*ix]/denrms;
            }
        }

        std::ostringstream outfile;
        outfile << "!" << name << ".fits.gz";
        fits_utils::write_fits_image(infl, nx, ny, outfile.str());

    }

} // namespace atmosphere
