///
/// @package phosim
/// @file telescopesetup.cpp
/// @brief setup for telescope (part of image class)
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <vector>
#include <iostream>
#include <fstream>
#include "ancillary/readtext.h"
#include "fea.h"

using readtext::readText;
using namespace fea;

int Image::telSetup () {

    double runningz;


    long seedchip = 0;
    for (size_t m = 0; m < chipid.size(); m++) {
        seedchip += static_cast<long>((static_cast<int>(chipid.c_str()[m]%10))*pow(10, m));
    }

    chip.nampx = (maxx - minx + 1);
    chip.nampy = (maxy - miny + 1);
    chip.buffer = 400;
    chip.midpoint = pixelsx/2;
    if (saturation == 1) {

        satupmap = static_cast<int*>(calloc(chip.nampx*chip.nampy, sizeof(int)));
        satdownmap = static_cast<int*>(calloc(chip.nampx*chip.nampy, sizeof(int)));
        for (long i = minx; i <= maxx; i++) {
            for (long j = miny; j <= maxy; j++) {
                *(satdownmap + chip.nampx*(j - miny) + (i - minx)) = i - 1;
                *(satupmap + chip.nampx*(j - miny) + (i - minx)) = i + 1;
            }
        }

    }
    chip.focal_plane = static_cast<float*>(calloc(chip.nampx*chip.nampy, sizeof(float)));
    if (opdfile) {
        opd = static_cast<double*>(calloc(SCREEN_SIZE/4*SCREEN_SIZE/4, sizeof(double)));
        opdcount = static_cast<double*>(calloc(SCREEN_SIZE/4*SCREEN_SIZE/4, sizeof(double)));
    }

    // OPTICS AND COATINGS
    fprintf(stdout, "Building Optics.\n");
    std::ostringstream opticsFile;
    opticsFile << instrdir << "/optics_" << filter << ".txt";
    readText opticsPars(opticsFile.str());
    size_t totSurf = opticsPars.getSize();

    surface.setup(totSurf, SURFACE_POINTS);
    coating.setup(totSurf);

    nsurf = 0;
    npertsurf = 0;
    runningz = 0.0;
    nmirror = 0;
    for (size_t t(0); t < totSurf; t++){
        std::istringstream iss(opticsPars[t]);
        std::string surfaceName, surfacetype, coatingFile, mediumFile;
        iss >> surfaceName >> surfacetype;

        if (surfacetype == "mirror") surface.surfacetype[nsurf] = MIRROR;
        else if (surfacetype == "lens") surface.surfacetype[nsurf] = LENS;
        else if (surfacetype == "filter") surface.surfacetype[nsurf] = FILTER;
        else if (surfacetype == "det") surface.surfacetype[nsurf] = DETECTOR;
        else if (surfacetype == "grating") surface.surfacetype[nsurf] = GRATING;
        if (surfacetype == "mirror") nmirror++;

        iss >> surface.radiusCurvature[nsurf];
        double dz;
        iss >> dz;
        runningz += dz;
        surface.height[nsurf] = runningz;
        iss >> surface.outerRadius[nsurf];
        iss >> surface.innerRadius[nsurf];
        iss >> surface.conic[nsurf];
        iss >> surface.three[nsurf];
        iss >> surface.four[nsurf];
        iss >> surface.five[nsurf];
        iss >> surface.six[nsurf];
        iss >> surface.seven[nsurf];
        iss >> surface.eight[nsurf];
        iss >> surface.nine[nsurf];
        iss >> surface.ten[nsurf];
        iss >> coatingFile;
        iss >> mediumFile;
        surface.centerx[nsurf] = 0.0;
        surface.centery[nsurf] = 0.0;
        surface.rmax[nsurf] = surface.outerRadius[nsurf];
        if (surfacetype != "none") {
            surface.asphere(nsurf, SURFACE_POINTS);
            if (surface.surfacetype[nsurf] == MIRROR || surface.surfacetype[nsurf] == DETECTOR || surface.surfacetype[nsurf] == LENS) {
                surface.surfacepert[nsurf] = npertsurf;
                npertsurf++;
            } else {
                surface.surfacepert[nsurf] = 0;
            }

            // COATINGS
            surface.surfacecoating[nsurf] = 0;
            if (coatingFile != "none") {
                readText coatingPars(instrdir + "/" + coatingFile);
                size_t nline = coatingPars.getSize();
                coating.allocate(nsurf, nline);
                long j = 0;
                long i = 0;
                double angle, angle0 = 0.0;
                for (size_t tt(0); tt<nline; tt++) {
                    std::istringstream isst(coatingPars[tt]);
                    isst >> angle;
                    if (i == 0 && j == 0) angle0 = angle;
                    if (angle > angle0) {
                        i++;
                        coating.wavelengthNumber[nsurf] = j;
                        if (j != coating.wavelengthNumber[nsurf] && i > 1) {
                            std::cout << "Error in format of " << coatingFile << std::endl;
                            exit(1);
                        }
                        j = 0;
                    }
                    *(coating.angle[nsurf] + i) = angle;
                    isst >> *(coating.wavelength[nsurf] + j);
                    isst >> *(coating.transmission[nsurf] + tt);
                    isst >> *(coating.reflection[nsurf] + tt);
                    j++;
                    angle0 = angle;
                }
                i++;
                coating.wavelengthNumber[nsurf] = j;
                coating.angleNumber[nsurf] = i;
                coating.wavelength[nsurf] = static_cast<double*>(realloc(coating.wavelength[nsurf], coating.wavelengthNumber[nsurf]*sizeof(double)));
                coating.angle[nsurf] = static_cast<double*>(realloc(coating.angle[nsurf], coating.angleNumber[nsurf]*sizeof(double)));
                surface.surfacecoating[nsurf] = 1;
            }

            // MEDIUM
            surface.surfacemed[nsurf] = 0;
            if (mediumFile == "air") surface.surfacemed[nsurf] = 2;
            if (mediumFile != "vacuum" && mediumFile != "air") {
                readText mediumPars(instrdir + "/" + mediumFile);
                size_t nline = mediumPars.getSize();
                medium.index_refraction_number[nsurf] = static_cast<long>(nline);
                medium.index_refraction[nsurf] = static_cast<double*>(calloc(nline, sizeof(double)));
                medium.index_refraction_wavelength[nsurf] = static_cast<double*>(calloc(nline, sizeof(double)));
                for (size_t tt(0); tt < nline; tt++) {
                    std::istringstream isst(mediumPars[tt]);
                    isst >> *(medium.index_refraction_wavelength[nsurf] + tt);
                    isst >> *(medium.index_refraction[nsurf] + tt);
                }
                surface.surfacemed[nsurf] = 1;
            }
            nsurf++;

        }
    }

    if (telescope_on==1) {
        maxr = surface.outerRadius[0];
        minr = surface.innerRadius[0];
    }

    fprintf(stdout, "Placing Obstructions.\n");
    obstruction.nspid = 0;
    readText spiderPars(instrdir + "/spider.txt");
    for (size_t t(0); t < spiderPars.getSize(); t++) {
        std::istringstream iss(spiderPars[t]);
        std::string spiderType;
        iss >> spiderType;
        iss >> obstruction.spider_height[obstruction.nspid];
        iss >> obstruction.spider_width[obstruction.nspid];
        iss >> obstruction.spider_center[obstruction.nspid];
        double tempf4;
        iss >> tempf4;
        obstruction.spider_depth[obstruction.nspid] = tempf4;
        obstruction.spider_angle[obstruction.nspid] = tempf4;
        obstruction.spider_reference[obstruction.nspid] = tempf4;
        if (spiderType=="crossx") obstruction.spider_type[obstruction.nspid] = 1;
        else if (spiderType=="crossy") obstruction.spider_type[obstruction.nspid] = 2;
        obstruction.nspid++;
    }



    /* PERTURBATIONS */

    fprintf(stdout, "Perturbing Design.\n");
    perturbation.zernike_coeff = static_cast<double*>(calloc(NZERN*npertsurf, sizeof(double)));
    // chuck's old model: 0.82312e-3*pow(zv, -1.2447)

    for (long i = 0; i < npertsurf; i++) {
        for (long j = 0; j < NZERN; j++) {
            *(perturbation.zernike_coeff + i*NZERN + j) = izernike[i][j];
        }
    }


    perturbation.zernikeflag = 0;

    perturbation.eulerPhi.reserve(nsurf + 1);
    perturbation.eulerPsi.reserve(nsurf + 1);
    perturbation.eulerTheta.reserve(nsurf + 1);
    perturbation.decenterX.reserve(nsurf + 1);
    perturbation.decenterY.reserve(nsurf + 1);
    perturbation.defocus.reserve(nsurf + 1);

    perturbation.zernike_r = static_cast<double*>(calloc(NZERN*SURFACE_POINTS, sizeof(double)));
    perturbation.zernike_phi = static_cast<double*>(calloc(NZERN*SURFACE_POINTS, sizeof(double)));
    perturbation.zernike_r_grid = static_cast<double*>(calloc(SURFACE_POINTS, sizeof(double)));
    perturbation.zernike_phi_grid = static_cast<double*>(calloc(SURFACE_POINTS, sizeof(double)));
    perturbation.zernike_normal_r = static_cast<double*>(calloc(NZERN*SURFACE_POINTS, sizeof(double)));
    perturbation.zernike_normal_phi = static_cast<double*>(calloc(NZERN*SURFACE_POINTS, sizeof(double)));


    if (pertType == "zern") {
        zernikes(perturbation.zernike_r, perturbation.zernike_phi, perturbation.zernike_r_grid,
                 perturbation.zernike_phi_grid, perturbation.zernike_normal_r, perturbation.zernike_normal_phi,
                 SURFACE_POINTS, NZERN);
    } else if (pertType == "chebyshev") {
        //will change this later
        zernikes(perturbation.zernike_r, perturbation.zernike_phi, perturbation.zernike_r_grid,
                 perturbation.zernike_phi_grid, perturbation.zernike_normal_r, perturbation.zernike_normal_phi,
                 SURFACE_POINTS, NZERN);
    }
    perturbation.zernike_summed = static_cast<double*>(calloc(npertsurf*SURFACE_POINTS*SURFACE_POINTS, sizeof(double)));
    perturbation.zernike_summed_nr_p = static_cast<double*>(calloc(npertsurf*SURFACE_POINTS*SURFACE_POINTS, sizeof(double)));
    perturbation.zernike_summed_np_r = static_cast<double*>(calloc(npertsurf*SURFACE_POINTS*SURFACE_POINTS, sizeof(double)));

    for (long j = 0; j < SURFACE_POINTS; j++) {
        for (long l = 0; l < SURFACE_POINTS; l++) {
            for (long k = 0; k < npertsurf; k++) {
                *(perturbation.zernike_summed + k*SURFACE_POINTS*SURFACE_POINTS + l*SURFACE_POINTS + j) = 0;
                for (long i = 0; i < NZERN; i++) {
                    *(perturbation.zernike_summed + k*SURFACE_POINTS*SURFACE_POINTS + l*SURFACE_POINTS + j) +=
                        *(perturbation.zernike_coeff + k*NZERN + i)*(*(perturbation.zernike_r + i*SURFACE_POINTS + j))*
                        (*(perturbation.zernike_phi + i*SURFACE_POINTS + l));
                }

            }
        }
    }

    for (long j = 0; j < SURFACE_POINTS; j++) {
        for (long l = 0; l < SURFACE_POINTS; l++) {
            for (long k = 0; k < npertsurf; k++) {
                *(perturbation.zernike_summed_nr_p + k*SURFACE_POINTS*SURFACE_POINTS + l*SURFACE_POINTS + j) = 0;
                for (long i = 0; i < NZERN; i++) {
                    *(perturbation.zernike_summed_nr_p + k*SURFACE_POINTS*SURFACE_POINTS + l*SURFACE_POINTS + j) +=
                        *(perturbation.zernike_coeff + k*NZERN + i)*(*(perturbation.zernike_normal_r + i*SURFACE_POINTS + j))*
                        (*(perturbation.zernike_phi + i*SURFACE_POINTS + l));
                }
            }
        }
    }

    for (long j = 0; j < SURFACE_POINTS; j++) {
        for (long l = 0; l < SURFACE_POINTS; l++) {
            for (long k = 0; k < npertsurf; k++) {
                *(perturbation.zernike_summed_np_r + k*SURFACE_POINTS*SURFACE_POINTS + l*SURFACE_POINTS + j) = 0;
                for (long i = 0; i < NZERN; i++) {
                    *(perturbation.zernike_summed_np_r + k*SURFACE_POINTS*SURFACE_POINTS + l*SURFACE_POINTS + j) +=
                        *(perturbation.zernike_coeff + k*NZERN + i)*(*(perturbation.zernike_r + i*SURFACE_POINTS + j))*
                        (*(perturbation.zernike_normal_phi + i*SURFACE_POINTS + l));
                }
            }
        }
    }


    perturbation.rotationmatrix = static_cast<double*>(calloc((nsurf + 1)*3*3, sizeof(double)));
    perturbation.inverserotationmatrix = static_cast<double*>(calloc((nsurf + 1)*3*3, sizeof(double)));
    for (long i = 0; i< nsurf + 1; i++) {
        perturbation.eulerPhi[i] = body[i][0];
        perturbation.eulerPsi[i] = body[i][1];
        perturbation.eulerTheta[i] = body[i][2];
        perturbation.decenterX[i] = body[i][3];
        perturbation.decenterY[i] = body[i][4];
        perturbation.defocus[i] = body[i][5];
    }
    perturbation.decenterX[nsurf] += centerx*1e-3;
    perturbation.decenterY[nsurf] += centery*1e-3;

    surface.centerx[nsurf] = 0.0;
    surface.centery[nsurf] = 0.0;
    surface.rmax[nsurf] = sqrt(pixelsx*pixelsx + pixelsy*pixelsy)*sqrt(2.0)/2.0*pixsize*1e-3;
    surface.height[nsurf] = runningz;
    for (long i = 0; i < nsurf + 1; i++) {
        *(perturbation.rotationmatrix + 9*i + 0*3 + 0) = cos(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 0*3 + 1) = cos(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 0*3 + 2) = sin(perturbation.eulerPsi[i])*sin(perturbation.eulerTheta[i]);
        *(perturbation.rotationmatrix + 9*i + 1*3 + 0) = -sin(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 1*3 + 1) = -sin(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.rotationmatrix + 9*i + 1*3 + 2) = cos(perturbation.eulerPsi[i])*sin(perturbation.eulerTheta[i]);
        *(perturbation.rotationmatrix + 9*i + 2*3 + 0) = sin(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i]);
        *(perturbation.rotationmatrix + 9*i + 2*3 + 1) = -sin(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i]);
        *(perturbation.rotationmatrix + 9*i + 2*3 + 2) = cos(perturbation.eulerTheta[i]);
        *(perturbation.inverserotationmatrix + 9*i + 0*3 + 0) = cos(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 0*3 + 1) = -sin(perturbation.eulerPsi[i])*cos(perturbation.eulerPhi[i]) -
            cos(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 0*3 + 2) = sin(perturbation.eulerTheta[i])*sin(perturbation.eulerPhi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 1*3 + 0) = cos(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 1*3 + 1) = -sin(perturbation.eulerPsi[i])*sin(perturbation.eulerPhi[i]) +
            cos(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 1*3 + 2) = -sin(perturbation.eulerTheta[i])*cos(perturbation.eulerPhi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 2*3 + 0) = sin(perturbation.eulerTheta[i])*sin(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 2*3 + 1) = sin(perturbation.eulerTheta[i])*cos(perturbation.eulerPsi[i]);
        *(perturbation.inverserotationmatrix + 9*i + 2*3 + 2) = cos(perturbation.eulerTheta[i]);
    }

    /* Update perturbation model from FEA file */
    size_t nNear=4;
    // int degree=1;
    int leafSize=10;
    for (long i=0;i<nsurf;i++) {
        if (feaflag[i]==1) {
            long k=surface.surfacepert[i];
            std::vector<double> ix(SURFACE_POINTS*SURFACE_POINTS), iy(SURFACE_POINTS*SURFACE_POINTS);
            for (long j=0;j<SURFACE_POINTS;j++) {
                for (long l=0;l<SURFACE_POINTS;l++) {
                    ix[l*SURFACE_POINTS+j]=surface.rmax[i]*perturbation.zernike_r_grid[j]*cos(perturbation.zernike_phi_grid[l]);
                    iy[l*SURFACE_POINTS+j]=surface.rmax[i]*perturbation.zernike_r_grid[j]*sin(perturbation.zernike_phi_grid[l]);
                }
            }
            Fea feaTree(feafile[i],leafSize,surface,i);
            feaTree.knnQueryFitDegree1(ix, iy, &perturbation.zernike_summed[k*SURFACE_POINTS*SURFACE_POINTS],
                    &perturbation.zernike_summed_nr_p[k*SURFACE_POINTS*SURFACE_POINTS],
                    &perturbation.zernike_summed_np_r[k*SURFACE_POINTS*SURFACE_POINTS],nNear);
            //feaTree.knnQueryFit(ix, iy, &perturbation.zernike_summed[k*SURFACE_POINTS*SURFACE_POINTS],
            //        &perturbation.zernike_summed_nr_p[k*SURFACE_POINTS*SURFACE_POINTS],
            //        &perturbation.zernike_summed_np_r[k*SURFACE_POINTS*SURFACE_POINTS],nNear, degree);
            feaTree.getTransformation(perturbation,i);

        }
    }


    /// Large Angle Scattering
    perturbation.miescatter = static_cast<double*>(calloc(10000, sizeof(double)));
    for (long j = 0; j < 10000; j++) {
        *(perturbation.miescatter + j) = 1.0/(1.0 + pow(((static_cast<float>(j))/10000.0*0.1*M_PI/180.0)/
                                                        (1.0*M_PI/180.0/3600.0), 3.5))*2.0*M_PI*
            (static_cast<float>(j));
    }
    double tempf1 = 0.0;
    for (long j = 0; j < 10000; j++) {
        tempf1 = tempf1 + *(perturbation.miescatter + j);
    }
    for (long j = 0; j < 10000; j++) {
        *(perturbation.miescatter + j) = *(perturbation.miescatter + j)/tempf1;
    }
    for (long j = 1; j < 10000; j++) {
        *(perturbation.miescatter + j) = *(perturbation.miescatter + j) + *(perturbation.miescatter + j - 1);
    }

    // Tracking

    rotationrate = 15.04*cos(latitude)*cos(azimuth)/cos(M_PI/2 - zenith);

    if (trackingfile == ".") tracking_on = 0;
    if (tracking_on) {

        readText trackingPars(trackingfile);
        size_t j = trackingPars.getSize();
        perturbation.jitterrot = static_cast<double*>(calloc(j, sizeof(double)));
        perturbation.jitterele = static_cast<double*>(calloc(j, sizeof(double)));
        perturbation.jitterazi = static_cast<double*>(calloc(j, sizeof(double)));
        screen.jitterwind = static_cast<double*>(calloc(j, sizeof(double)));
        perturbation.jittertime = static_cast<double*>(calloc(j, sizeof(double)));
        trackinglines = j;

        for (size_t t(0); t < j; t++) {
            std::istringstream iss(trackingPars[t]);
            iss >> perturbation.jittertime[t];
            double f1, f2, f3, f4;
            iss >> f1 >> f2 >> f3 >> f4;
            perturbation.jitterele[t] = f1*elevationjitter;
            perturbation.jitterazi[t] = f2*azimuthjitter;
            perturbation.jitterrot[t] = f3*rotationjitter;
            screen.jitterwind[t] = f4*windjitter;
        }

    }

    /// Setup silicon
    fprintf(stdout, "Electrifying Devices.\n");
    silicon.setup(ccdtemp, nbulk, nf, nb, sf, sb, siliconthickness, overdepbias, instrdir, chip.nampx, chip.nampy, pixsize, seedchip);

    fprintf(stdout, "Contaminating Surfaces.\n");
    contamination.setup(surface.innerRadius, surface.outerRadius, nsurf - 1, SURFACE_POINTS,
                        pixsize, chip.nampx, chip.nampy, qevariation, seedchip);


    fprintf(stdout, "Diffracting.\n");
    // 2nd kick
    {

        Vector position;
        Vector angle;
        double shiftedAngle;
        double r,  phi;
        double rindex;
        long index;
        double radius;
        long atmtempdebug = 0;
        double bestvalue;
        double bestscale = 1.0;
        Vector largeAngle;
        double stdx = 0.0, stdy = 0.0;
        double radx = 0.0, rady = 0.0;
        double ncount = 0.0;
        screen.hffunc = static_cast<double*>(calloc(10000, sizeof(double)));
        screen.hffunc_n = static_cast<double*>(calloc(10000, sizeof(double)));

        atmtempdebug = atmdebug;

        for (long l = 0; l <= atmtempdebug; l++) {
            atmdebug = l;

            for (long i = 0; i < SCREEN_SIZE; i++) {
                for (long j = 0; j < SCREEN_SIZE; j++) {
                    screen.tfocalscreen[i*SCREEN_SIZE + j] = 0;
                }
            }

            prtime = -1.0;

            for (long k = 0; k < 10; k++) {

                wavelength = 0.5;
                wavelengthFactor = pow(wavelength, -0.2)/screen.wavelengthfactor_nom;
                // time = RngDouble()*400.0-200.0;
                time = RngDouble()*exptime;
                shiftedAngle = spiderangle + time*rotationrate*ARCSEC;
                r = sqrt(RngDouble()*(maxr*maxr - minr*minr) + minr*minr);
                phi = RngDouble()*2*M_PI;
                position.x = r*cos(phi);
                position.y = r*sin(phi);
                index = find_linear(const_cast<double *>(&surface.radius[0]), SURFACE_POINTS, r, &rindex);
                position.z = interpolate_linear(const_cast<double *>(&surface.profile[0]), index, rindex);
                xp = position.x;
                yp = position.y;
                angle.x = 0.0;
                angle.y = 0.0;
                angle.z = -1.0;

                atmosphereIntercept(&position, angle, -1, 2);
                for (long layer = 0; layer < natmospherefile; layer++) {
                    atmosphereIntercept(&position, angle, layer, 2);
                    atmosphereRefraction(&angle, layer, 3);
                }
                atmosphereDiffraction(&angle);
                for (long i = 0; i < SCREEN_SIZE; i++) {
                    for (long j = 0; j < SCREEN_SIZE; j++) {
                        screen.tfocalscreen[i*SCREEN_SIZE + j] += screen.focalscreen[i*SCREEN_SIZE + j];
                    }
                }


            }



            for (long i = 0; i < 10000; i++) {
                *(screen.hffunc + i) = 0.0;
            }
            for (long i = 0; i < 10000; i++) {
                *(screen.hffunc_n + i) = 0.0;
            }
            for (long i = 0; i < SCREEN_SIZE; i++) {
                for (long j = 0; j < SCREEN_SIZE; j++) {
                    radius = sqrt(pow((i - (SCREEN_SIZE/2 + 1)), 2) + pow((j - (SCREEN_SIZE/2 + 1)), 2));
                    *(screen.hffunc + ((long)(radius))) += screen.tfocalscreen[i*SCREEN_SIZE + j];
                    *(screen.hffunc_n + ((long)(radius))) += 1;
                }
            }

            for (long j = 0; j < 10000; j++) {
                if (screen.hffunc_n[j] < 1) {
                    *(screen.hffunc + j) = 0;
                } else {
                    *(screen.hffunc + j) = *(screen.hffunc + j)/(screen.hffunc_n[j])*((double)(j));
                }
            }
            double tempf1 = 0.0;
            for (long j = 0; j < 10000; j++) {
                tempf1 = tempf1 + *(screen.hffunc + j);
            }
            for (long j = 0; j < 10000; j++) {
                *(screen.hffunc + j) = *(screen.hffunc + j)/tempf1;
            }
            for (long j = 1; j < 10000; j++) {
                *(screen.hffunc + j) = *(screen.hffunc + j) + *(screen.hffunc + j - 1);
            }

            // FILE *fileptr;
            // fileptr=fopen("temp.txt","w");
            // for (long j = 0 ; j < 10000 ; j++) fprintf(fileptr,"%e\n",*(screen.hffunc+j));
            // fclose(fileptr);
            // FILE *fileptr2;
            // fileptr2=fopen("temp2.txt","w");
            // for (long j = 0 ; j < 10000 ; j++) fprintf(fileptr2,"%e\n",*(perturbation.miescatter+j));
            // fclose(fileptr2);

            if (l == 0) {

                bestvalue = 1e10;
                bestscale = 1.0;
                for (double scales = 0.0; scales <= 1.01; scales += 0.01) {
                    stdx = 0.0;
                    stdy = 0.0;
                    ncount = 0.0;

                    for (long k = 0; k < 10000; k++) {

                        screen.secondKickSize = scales;
                        wavelength = 0.5;
                        wavelengthFactor = pow(wavelength, -0.2)/screen.wavelengthfactor_nom;
                        time = RngDouble()*exptime;
                        shiftedAngle = spiderangle + time*rotationrate*ARCSEC;
                        r = sqrt(RngDouble()*(maxr*maxr - minr*minr) + minr*minr);
                        phi = RngDouble()*2*M_PI;
                        position.x = r*cos(phi);
                        position.y = r*sin(phi);
                        index = find_linear(const_cast<double *>(&surface.radius[0]), SURFACE_POINTS, r, &rindex);
                        position.z = interpolate_linear(const_cast<double *>(&surface.profile[0]), index, rindex);
                        xp = position.x;
                        yp = position.y;
                        angle.x = 0.0;
                        angle.y = 0.0;
                        angle.z = -1.0;
                        largeAngle.x = 0;
                        largeAngle.y = 0;

                        secondKick(&largeAngle);
                        atmosphereIntercept(&position, angle, -1, 1);
                        for (long layer = 0; layer < natmospherefile; layer++) {
                            atmosphereIntercept(&position, angle, layer, 1);
                            atmosphereRefraction(&angle, layer, 1);
                        }
                        angle.x = angle.x + largeAngle.x;
                        angle.y = angle.y + largeAngle.y;

                        radx = angle.x/((totalseeing + 1e-6)/2.35482*ARCSEC*pow(1/cos(zenith), 0.6)*wavelengthFactor);
                        rady = angle.y/((totalseeing + 1e-6)/2.35482*ARCSEC*pow(1/cos(zenith), 0.6)*wavelengthFactor);
                        stdx += radx*radx*exp(-(radx*radx + rady*rady)/2.0/1.0);
                        stdy += rady*rady*exp(-(radx*radx + rady*rady)/2.0/1.0);
                        ncount += exp(-(radx*radx + rady*rady)/2.0/1.0);

                    }

                    stdx /= ncount;
                    stdy /= ncount;
                    screen.secondKickSize = sqrt(stdx + stdy);
                    if (fabs(screen.secondKickSize - 1.0) < bestvalue) {
                        bestvalue = fabs(screen.secondKickSize - 1.0);
                        bestscale = scales;
                    }

                }
                screen.secondKickSize = bestscale;
            }

        }

        atmdebug = atmtempdebug;


    }



    nphot = static_cast<long long>(M_PI*(maxr/10.0*maxr/10.0 - minr/10.0*minr/10.0)*exptime*(totalnorm/(6.626e-27)));


    settings();

    fprintf(stdout, "Number of Sources: %ld\n", nsource);
    fprintf(stdout, "Photons: %6.2e   Flux: %6.2e ergs/cm2/s\n", static_cast<double>(nphot), totalnorm);


    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "Photon Raytrace\n");
    readText versionPars(bindir + "/version");
    fprintf(stdout, "%s\n", versionPars[0].c_str());
    fprintf(stdout, "------------------------------------------------------------------------------------------\n");

    return(0);

}
