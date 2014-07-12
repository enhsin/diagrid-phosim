///
/// @package phosim
/// @file silicon.cpp
/// @brief silicon class
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

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string>
#include <fstream>
#include <vector>

#include <fitsio.h>
#include <fitsio2.h>

#include "silicon.h"
#include "constants.h"
#include "parameters.h"
#include "raytrace/basic_types.h"
#include "rng_mwc.h"
#include "ancillary/readtext.h"

using readtext::readText;
using namespace RandomNumbers;

double Silicon::abs_coeff(double lambda,double T) {

    // follows Rajkanan et al. Solid-state electronics 22, pp793-795.

    double Eg_T[2],Eg_0[]={1.1557f,2.5f}; // eV
    double Egd_T=0.0,Egd_0=3.2f;          // eV
    double Ep[]={1.827e-2f, 5.773e-2f}; // eV
    double C[]={5.5f,4.0f};            // no dimension
    double A[]={3.231e2f,7.237e3f};    // cm-1 eV-2
    double Ad = 1.052e6f;             // cm-1 eV-2
    double k_boltzmann=8.617e-5f;     // eV K-1
    double eV;                              // eV
    double beta = 7.021e-4f;                 // eV K-1
    double gamma=1108;                      // K
    double alpha=0.0;
    double delta_e0,delta_e1[2][2][2]; // eV

    Egd_T   = Egd_0   - (beta*T*T/(T + gamma));
    Eg_T[0] = Eg_0[0] - (beta*T*T/(T + gamma));
    Eg_T[1] = Eg_0[1] - (beta*T*T/(T + gamma));

    if (lambda>3100) {
        eV=12398/lambda;
        delta_e0=eV - Egd_T;
        for (int j=0;j<2;j++) {
            for (int i=0;i<2;i++) {
                delta_e1[i][j][0] = eV - Eg_T[j] + Ep[i];
                delta_e1[i][j][1] = eV - Eg_T[j] - Ep[i];
            }
        }

        alpha=Ad*sqrt((delta_e0>0)?delta_e0:0);
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                alpha += C[i]*A[j]*
                    (((delta_e1[i][j][0]>0)?pow(delta_e1[i][j][0],2)/(exp(Ep[i]/(k_boltzmann*T)) - 1):0) +
                     ((delta_e1[i][j][1]>0)?pow(delta_e1[i][j][1],2)/(1 - exp(-Ep[i]/(k_boltzmann*T))):0));
            }
        }
        // alpha is in cm-1.
    } else {
        // wavelength range is between 413 & 3100 Angstroms - not well characterized. use a constant.
        alpha=1.35e6;
    }

    return(alpha);
}

double Silicon::siIndexRefraction (double lambda) {
    double energy=12398/lambda;
    // data points were digitized from Philipp & Taft (1960) and fit between 0 and 3.4 eV
    return(double)(3.364967 + 0.2119184*energy + 2.78878*exp((energy - 3.30168)/0.397862));
}

double Silicon::dope_profile(double z,double N_bulk, double N_f, double N_b, double s_f, double s_b, double t_si) {
    return(N_bulk + N_b*exp(-(t_si - z)/s_b) + N_f*exp(-(z)/s_f));
}


double Silicon::mu_Si (double E, double T, int polarity) {

    double vm, Ec, beta;
    // Jacobini et al. (1977) equation (9)
    if (polarity == -1) {
        vm = 1.53e9 * pow(T, -0.87); // cm/s
        Ec = 1.01 * pow(T, 1.55); // V/cm
        beta = 2.57e-2 * pow(T, 0.66); // index
    } else {
        vm = 1.62e8 * pow(T, -0.52); // cm/s
        Ec = 1.24 * pow(T, 1.68); // V/cm
        beta = 0.46 * pow(T, 0.17); // index
    }
    return((vm/Ec)/pow(1 + pow(fabs(E)/Ec,beta),1/beta));

}


void Silicon::setup (double ccdtemp, double N_bulk, double N_f, double N_b,
                     double s_f, double s_b, double t_si, double overdep_bias,
                     std::string instrdir,
                     long nampx, long nampy,double pixsize, long seedchip) {


    s_f    =  s_f/1e4f;
    s_b    =  s_b/1e4f;
    t_si   =  t_si/1e4f;


    numWavelength = 1024;
    numTemperature = 16;
    numDopant = 16;
    numThickness = SILICON_STEPS;

    wavelengthGrid = (double*)calloc(numWavelength,sizeof(double));
    if (wavelengthGrid == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    temperatureGrid = (double*)calloc(numTemperature,sizeof(double));
    if (temperatureGrid == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    rho = (double*)calloc(numTemperature,sizeof(double));
    if (rho == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    dopantGrid = (double*)calloc(numDopant,sizeof(double));
    if (dopantGrid == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    thicknessGrid = (double*)calloc(numThickness,sizeof(double));
    if (thicknessGrid == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    meanFreePath = (double*)calloc(numWavelength*numTemperature,sizeof(double));
    if (meanFreePath == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    indexRefraction = (double*)calloc(numWavelength,sizeof(double));
    if (indexRefraction == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    sigma = (float*)calloc(numTemperature*numThickness*numDopant,sizeof(float));
    if (sigma == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    hsigma = (float*)calloc(numTemperature*numThickness*numDopant,sizeof(float));
    if (hsigma == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    fsigma = (float*)calloc(numThickness*numDopant,sizeof(float));
    if (fsigma == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    gsigma = (float*)calloc(numThickness*numDopant,sizeof(float));
    if (gsigma == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    sigmaX = (float*)calloc(nampx*nampy,sizeof(float));
    if (sigmaX == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    sigmaY = (float*)calloc(nampx*nampy,sizeof(float));
    if (sigmaY == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    gammaX = (float*)calloc(nampx*nampy,sizeof(float));
    if (gammaX == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    gammaY = (float*)calloc(nampx*nampy,sizeof(float));
    if (gammaY == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    deltaX = (float*)calloc(nampx*nampy,sizeof(float));
    if (deltaX == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    deltaY = (float*)calloc(nampx*nampy,sizeof(float));
    if (deltaY == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    nbulkmap = (float*)calloc(nampx*nampy,sizeof(float));
    if (nbulkmap == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}
    deadLayer = (float*)calloc(nampx*nampy,sizeof(float));
    if (deadLayer == NULL) {fprintf(stderr,"Allocation error\n"); exit(1);}

    double period;
    double amplitude;
    double ratio;
    double surfaceCharge;
    double depth;
    double width;
    double height;
    double xoverlap;
    double yoverlap;
    double pixelVarX;
    double pixelVarY;
    int polarity;
    std::vector<double> xinit, yinit, angle;
    readText pars(instrdir + "/silicon.txt");
    for (size_t t(0); t < pars.getSize(); t++) {
        std::string line(pars[t]);
        readText::get(line, "treeRingPeriod", period);
        readText::get(line, "treeRingAmplitude", amplitude);
        readText::get(line, "treeRingRatio", ratio);
        readText::get(line, "edgeSurfaceCharge", surfaceCharge);
        readText::get(line, "deadLayerDepth", depth);
        readText::get(line, "deadLayerWidth", width);
        readText::get(line, "deadLayerHeight", height);
        readText::get(line, "deadLayerXoverlap", xoverlap);
        readText::get(line, "deadLayerYoverlap", yoverlap);
        readText::get(line, "pixelVarX", pixelVarX);
        readText::get(line, "pixelVarY", pixelVarY);
        readText::get(line, "deadLayerXinit", xinit);
        readText::get(line, "deadLayerYinit", yinit);
        readText::get(line, "deadLayerAngle", angle);
        readText::get(line, "channelDepth", channelDepth);
        readText::get(line, "stopMomentPerPixel", stopMomentPerPixel);
        readText::get(line, "siliconType", polarity);
    }


    // tree ring map
    RngSetSeed32_reseed(1000 + seedchip);
    double x0 = RngDouble_reseed()*(ratio*nampx) - (ratio - 1)/2*nampx;
    double y0 = RngDouble_reseed()*(ratio*nampy) - (ratio - 1)/2*nampy;
    for (long i = 0;i<nampx;i++) {
        for (long j = 0;j<nampy;j++) {
            nbulkmap[nampx*j + i] = 1 + amplitude*sin(2*M_PI*sqrt(pow(i - x0,2) + pow(j - y0,2))/period);
        }
    }

    // lateral deflection map
    double scale = fabs(surfaceCharge)/N_bulk;
    if (scale == 0.0) scale = 1e-4*pixsize;
    for (long i = 0;i<nampx;i++) {
        for (long j = 0;j<nampy;j++) {
            deltaX[nampx*j + i] = surfaceCharge/EPSILON_0/EPSILON_SI*exp(-(i*pixsize*1e-4)/scale)*E_CHARGE-
                surfaceCharge/EPSILON_0/EPSILON_SI*exp(-((nampx - 1 - i)*pixsize*1e-4)/scale)*E_CHARGE;
            deltaY[nampx*j + i] = surfaceCharge/EPSILON_0/EPSILON_SI*exp(-(j*pixsize*1e-4)/scale)*E_CHARGE-
                surfaceCharge/EPSILON_0/EPSILON_SI*exp(-((nampy - 1 - j)*pixsize*1e-4)/scale)*E_CHARGE;
        }
    }

    // add tree rings
    for (long i = 1;i<(nampx - 1);i++) {
        for (long j = 1;j<(nampy - 1);j++) {
            sigmaX[nampx*j + i] = (nbulkmap[nampx*j + i + 1] - nbulkmap[nampx*j + i - 1])*
                N_bulk*(2.0*t_si)*(period*pixsize*1e-4)/(2.0*M_PI)/(2.0*M_PI)*E_CHARGE/
                EPSILON_SI/EPSILON_0/2.0/(pixsize*1e-4)/2.0;
            sigmaY[nampx*j + i] = (nbulkmap[nampx*(j + 1) + i] - nbulkmap[nampx*(j - 1) + i])*
                N_bulk*(2.0*t_si)*(period*pixsize*1e-4)/(2.0*M_PI)/(2.0*M_PI)*E_CHARGE/
                EPSILON_SI/EPSILON_0/2.0/(pixsize*1e-4)/2.0;
        }
    }

    // int status;
    // long naxesa[2];
    // fitsfile *faptr;
    // status = 0;
    // fits_create_file(&faptr,"!defl.fits",&status);
    // naxesa[0] = nampx; naxesa[1] = nampy;
    // fits_create_img(faptr,FLOAT_IMG,2,naxesa,&status);
    // fits_write_img(faptr,TFLOAT,1,nampx*nampy,deltaX,&status);
    // fits_close_file(faptr,&status);

    // nonuniform pixels
    for (long i = 0;i<nampx;i++) {
        for (long j = 0;j<nampy;j++) {
            gammaX[nampx*j+i] = random_gaussian_reseed()*pixelVarX*pixsize;
            gammaY[nampx*j + i] = random_gaussian_reseed()*pixelVarY*pixsize;
        }
    }

    // dead layer
    for (size_t k = 0; k < xinit.size(); k++) {
        for (long i = 0; i < nampx; i++) {
            for (long j = 0; j < nampy; j++) {
                double ip = i*cos(angle[k]*M_PI/180.) - j*sin(angle[k]*M_PI/180.);
                double jp = i*sin(angle[k]*M_PI/180.)+j*cos(angle[k]*M_PI/180.);
                double di = fmod(ip*pixsize - xinit[k], width - xoverlap);
                double dj = fmod(jp*pixsize - yinit[k], height - yoverlap);
                while (di<0) di+=width - xoverlap;
                while (dj<0) dj+=height - yoverlap;
                if (di<width && dj<height) {
                    deadLayer[nampx*j+i]+=depth;
                    if (di<xoverlap) deadLayer[nampx*j+i]+=depth;
                    if (dj<yoverlap) deadLayer[nampx*j+i]+=depth;
                    if (di<xoverlap && dj<yoverlap) deadLayer[nampx*j+i]+=depth;
                }
            }
        }
    }

    for (long i = 0;i<numTemperature;i++) {
        for (long j = 0;j<numWavelength;j++) {
            temperatureGrid[i] = ccdtemp+((double)i)/((double)numTemperature)*2.0 - 1.0;
            rho[i] = t_si*((double)(i + 1))/((double)(numTemperature))/5.0;
            wavelengthGrid[j] = 0.3 + ((double)j)/((double)(numWavelength - 1))*0.9;
            meanFreePath[i*numWavelength + j] = 10.0/abs_coeff(wavelengthGrid[j]*10000.0,temperatureGrid[i]);
        }
    }
    for (long i = 0;i<numDopant;i++) dopantGrid[i] = ((double)i)/((double)numDopant)*2*amplitude*N_bulk+
                                       (1.0 - amplitude)*N_bulk;

    for (long i = 0;i<numWavelength;i++) {
        indexRefraction[i] = siIndexRefraction(wavelengthGrid[i]*10000.0);
    }

    double z[SILICON_STEPS],E[SILICON_STEPS],v[SILICON_STEPS],
        tcol[SILICON_STEPS],tp[SILICON_STEPS], gp[SILICON_STEPS], hp[SILICON_STEPS];

    double min_E, dz;
    int i, min_E_index;
    long n = SILICON_STEPS;
    double factor = 3.0*E_CHARGE/EPSILON_0/EPSILON_SI/(4*M_PI)*channelDepth;

    for (long l = 0;l<numDopant;l++) {
        for (long k = 0;k<numTemperature;k++) {
            dz = (-t_si)/(double)(n - 1);
            for (i = 0;i<n;i++) {
                z[i] = t_si + i*dz;
                thicknessGrid[i] = t_si + i*dz;
                if (i == 0) {
                    E[0] = overdep_bias/t_si;
                    v[0] = 0;
                } else {
                    // sample the doping profile to evaluate E[i]
                    int n_sample = 100;
                    double dp = 0;
                    int j;
                    dp = 0;
                    for (j = 0;j<n_sample;j++) {
                        dp +=  dope_profile(t_si+dz*(i + (j)/(double)(n_sample)),dopantGrid[l],N_f,N_b,s_f,s_b,t_si);
                    }
                    dp /= n_sample;
                    E[i] = E[i - 1] + dp*E_CHARGE/(EPSILON_0*EPSILON_SI)*dz;
                    v[i] = v[i - 1]-0.5f*(E[i] + E[i - 1])*dz;
                }
            }
            min_E_index = 0;
            min_E = 0;
            for (i = 0;i<n;i++) {
                // find minimum value for E and count from that index
                if (E[i]<min_E) {
                    min_E_index = i;
                    min_E = E[i];
                }
            }
            i = n - 1;
            tcol[i] = (dz/(min_E*mu_Si(min_E,temperatureGrid[k], polarity)));
            tp[i] = dz/min_E*((float)i)/((float)(n - 1))*(1.0 - ((float)i)/((float)(n - 1)))*4.0;
            gp[i] = dz/min_E;
            hp[i] = dz/min_E*factor*(z[i]*rho[k])/pow(z[i]*z[i] + rho[k]*rho[k],5./2.);
            while (i--) {
                if (i>min_E_index) {
                    tcol[i] = (dz/(min_E*mu_Si(min_E,temperatureGrid[k], polarity)));
                    tp[i] = dz/min_E*((float)i)/((float)(n - 1))*(1.0 - ((float)i)/((float)(n - 1)))*4.0;
                    gp[i] = dz/min_E;
                    hp[i] = dz/min_E*factor*(z[i]*rho[k])/pow(z[i]*z[i] + rho[k]*rho[k],5./2.);
                } else {
                    if (E[i] + E[i + 1]>0) {
                        tcol[i] = 0.1;
                        tp[i] = 0.1*((float)i)/((float)(n - 1))*(1.0 - ((float)i)/((float)(n - 1)))*4.0;
                        gp[i] = 0.1;
                        hp[i] = 0.1*factor*(z[i]*rho[k])/pow(z[i]*z[i] + rho[k]*rho[k],5./2.);
                    } else {
                        tcol[i] = 2*dz/(E[i]*mu_Si(E[i],temperatureGrid[k], polarity) + E[i + 1]*mu_Si(E[i + 1],temperatureGrid[k], polarity));
                        tp[i] = 2*dz/(E[i] + E[i + 1])*((float)i)/((float)(n - 1))*(1.0 - ((float)i)/((float)(n - 1)))*4.0;
                        gp[i] = 2*dz/(E[i] + E[i + 1]);
                        hp[i] = 2*dz/(E[i] + E[i + 1])*factor*(z[i]*rho[k])/pow(z[i]*z[i] + rho[k]*rho[k],5./2.);
                    }
                }
                tcol[i]+=tcol[i+1];
                tp[i]+=tp[i + 1];
                gp[i]+=gp[i + 1];
                hp[i]+=hp[i + 1];
            }
            i = n;
            while (i--) {
                sigma[l*numThickness*numTemperature + k*numThickness + i] = sqrt(2*K_BOLTZMANN*temperatureGrid[k]
                                                                                 /E_CHARGE*mu_Si(0,temperatureGrid[k], polarity) * tcol[i])*1e4f;
                fsigma[l*numThickness + i] = tp[i]*1e4f;
                gsigma[l*numThickness + i] = gp[i]*1e4f;
                hsigma[l*numThickness*numTemperature + k*numThickness + i] = hp[i]*1e4f;
            }

        }
    }

}
