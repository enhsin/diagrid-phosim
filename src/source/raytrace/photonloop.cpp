///
/// @package phosim
/// @file photonloop.cpp
/// @brief main photon loop
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include "event.h"
#include "grating.h"
#include "counter.h"

int Image::photonLoop () {

    Vector angle, largeAngle, position, positionDiffraction;
    Vector positionPrevious, anglePrevious, normal;
    double transmission, reflection, distance;

    long long ray;
    long long photonSource;
    long sourceCounter;

    long long detRay;
    long long backSourceOver;
    double sourceSaturation;
    double backSigma = 0.0;
    long preGhost = 0;
    long sourceOver;
    long newSurf = 0;
    long oldSurf;
    long waveIndex;
    long straylightcurrent;

    int miss;
    int initFlag;
    char tempstring[4096];

    EventFile* pEventLogging;
    Clog counterLog;
    Tlog throughputLog;
    Grating* pGrating;

    //    OPTIMIZATION INITIALIZATION
    maxcounter = 0;
    dynamicTransmission = (double**)malloc((natmospherefile*2 + nsurf*2 + 2)*sizeof(double*));
    for (int i = 0; i < (natmospherefile*2 + nsurf*2 + 2); i++) dynamicTransmission[i] = (double*)malloc(901*sizeof(double));
    for (int i = 0; i < (natmospherefile*2 + nsurf*2 + 2); i++)
        for (int j = 0; j < 901; j++) dynamicTransmission[i][j] = -1.0;
    if (checkpointcount  !=  0) readCheckpoint(checkpointcount);

    //    LOGGING INITIALIZATION
    if (eventfile) {
        pEventLogging = new EventFile((int)(nphot*100), outputdir, eventFitsFileName);
    } else {
        pEventLogging = NULL;
    }
    pGrating = new Grating();
    if (throughputfile) initThroughput(&throughputLog, nsurf);
    counterInit(&counterLog);

    //    SOURCE TYPE LOOP
    for (int sourceType = 0; sourceType < 73; sourceType++) {
        sourceCounter = 0;

        if (static_cast<int>(round(sourceType*checkpointtotal/72)) == checkpointcount) {

            //    SOURCE LOOP
            for (long source = 0; source<nsource; source++) {

                if ((sourceType >= 0 && sourceType <= 9 && sources.type[source] == 0 && (source%10) == sourceType) ||
                    (sourceType >= 10 && sourceType <= 19 && sources.type[source] == 1 && (source%10) == (sourceType - 10)) ||
                    (sourceType >= 20 && sourceType <= 29 && sources.type[source] == 2 && (source%10) == (sourceType - 20)) ||
                    (sourceType == 30 && sources.type[source] == 3 && sources.mag[source] > 40.5) ||
                    (sourceType >= 31 && sourceType <= 71 && sources.type[source] == 3 && sources.mag[source] > (70.5 - sourceType)
                     && sources.mag[source] <= (71.5 - sourceType))
                    || (sourceType == 72 && sources.type[source] == 3 && sources.mag[source] <= -0.5)) {

                    //    SETUP PHOTONS
                    photonSource = (long long)(nphot*sources.norm[source]/totalnorm);
                    if (nsource == 1) photonSource = nphot;
                    photonSource = poisson(photonSource);
                    if (telconfig != 0 && sources.type[source] != 0) photonSource = 0;

                    ray = 0;
                    sourceOver = 1;
                    sourceOver_m = 1;
                    sourceSaturation = 1.0;
                    sourceSaturationRadius = 0.0;
                    obstruction.pupil = 0;
                    prtime = -1.0;
                    op = 0.0;
                    if (sources.type[source]<3) {
                        double rate = ((double)counterLog.accepted + 10)/((double)counterLog.totalPhoton + 1000);
                        if (rate < 0.01) rate = 0.01;
                        if (backGamma < 1.0) backSourceOver = (long long)(backAlpha*sqrt(rate*photonSource));
                        else backSourceOver = (long long)(backAlpha/backGamma*sqrt(rate*photonSource));
                        if (backSourceOver < 1) backSourceOver = 1;
                        if (backSourceOver > photonSource) backSourceOver = photonSource;
                    } else {
                        backSourceOver = 1;
                    }
                    if (sources.type[source] >= 3) backSigma = 0.0; else backSigma = backBeta*backRadius/3600.0*platescale/1000.0;
                    if (sources.mag[source] < straylightcut && straylight == 1) straylightcurrent = 1; else straylightcurrent = 0;

                    //   MAIN PHOTON LOOP
                photonloop:
                    while (ray<photonSource) {
                        sourceOver = 1;
                        sourceOver_m = 1;

                        //   Get Wavelength and Time
                        miss = getWavelengthTime(&wavelength, &time, source);
                        // Get Angle
                        getAngle(&angle, time, source);
                        if (eventfile){
                            pEventLogging->logPhoton(angle.x, angle.y, wavelength, 0);
                            pEventLogging->logPhoton(time, 0, 0, 1);
                        }

                        //   Saturation Optimization
                        if (saturation) {
                            if (sourceSaturation > 1.0 && sourceSaturationRadius > 0.0 && sources.type[source] >= 3) {
                                sourceOver = round(sourceSaturation);
                                if (sourceOver > well_depth) sourceOver = well_depth;
                                if (sourceOver < 1) sourceOver = 1;
                                sourceOver_m = round(((sourceSaturation - np)/(1 - np)));
                                if (sourceOver_m < 1) {
                                    sourceOver_m = 1;
                                    sourceOver = 1;
                                }
                            }
                        }
                        if (sources.type[source] < 3 && backGamma > 1.0) {
                            double rate = ((double)counterLog.accepted + 10)/((double)counterLog.totalPhoton + 1000);
                            if (rate < 0.01) rate = 0.01;
                            if (backAlpha <= 0.0) sourceOver = static_cast<long long>(backGamma*sqrt(rate*photonSource));
                            else sourceOver = (long long)(backGamma);
                            sourceOver_m = sourceOver;
                            if (backSourceOver*sourceOver > photonSource) {
                                sourceOver = (long long)(photonSource/backSourceOver);
                                sourceOver_m = (long long)(photonSource/backSourceOver);
                            }
                            if (sourceOver < 1) {
                                sourceOver = 1;
                                sourceOver_m = 1;
                            }
                        }


                        if (miss) {
                            countBad(&counterLog, sourceOver*backSourceOver, &ray);
                            goto photonloop;
                        }
                        waveIndex = ((long)(wavelength*1000 - 300));
                        initFlag = 0;
                        if (throughputfile) addThroughput(&throughputLog, -1, waveIndex, sourceOver*backSourceOver);

                        // Dynamic Transmission Optimization
                        long kmax;
                        if (sources.type[source] < 3 && backGamma > 1.0) kmax = 1; else kmax = sourceOver;
                        for (long k = 0; k < kmax; k++) {

                            if ((k == 0) || (k > 0 && straylightcurrent == 1)) {
                            redoghlp:;
                                newSurf = -1;
                                photon.direction = 1;
                                counter = -1;
                                preGhost = 0;

                                for (int i = -1; i < natmospherefile; i++) {

                                    if (i >= 0) {
                                        if (counter >= (MAX_BOUNCE - 1)) goto maxbounce;
                                        if (transmissionPreCheck(i*2 + 1, waveIndex)) {
                                            if (k > 0) goto redoghlp;
                                            countBad_dt(&counterLog, sourceOver*backSourceOver, &ray);
                                            goto photonloop;
                                        }
                                    }

                                    if (counter >= (MAX_BOUNCE - 1)) goto maxbounce;
                                    if (transmissionPreCheck(i*2 + 2, waveIndex)) {
                                        if (k > 0) goto redoghlp;
                                        countBad_dt(&counterLog, sourceOver*backSourceOver, &ray);
                                        goto photonloop;
                                    }

                                }
                                while (1) {
                                    if (photon.direction == 1) newSurf++; else newSurf--;
                                    if (newSurf == -1) {
                                        if (k > 0) goto redoghlp;
                                        countBad_dt(&counterLog, sourceOver*backSourceOver, &ray);
                                        goto photonloop;
                                    }

                                    if (throughputfile == 1 && photon.direction == 1) {
                                        addThroughput(&throughputLog, newSurf, waveIndex, sourceOver*backSourceOver);
                                    }
                                    if (counter >= (MAX_BOUNCE - 1)) goto maxbounce;
                                    if (transmissionPreCheck(natmospherefile*2 + 1 + 2*newSurf, waveIndex)) {
                                        if (straylightcurrent == 1) {
                                            if (counter>= (MAX_BOUNCE - 1)) goto maxbounce;
                                            if (transmissionPreCheck(natmospherefile*2 + 1 + 2*newSurf + 1, waveIndex)) {
                                                if (k > 0) goto redoghlp;
                                                countBad_dt(&counterLog, sourceOver*backSourceOver, &ray);
                                                goto photonloop;
                                            } else {
                                                photon.direction = -photon.direction;
                                                preGhost++;
                                            }
                                        } else {
                                            if (k > 0) goto redoghlp;
                                            countBad_dt(&counterLog, sourceOver*backSourceOver, &ray);
                                            goto photonloop;
                                        }
                                    }
                                    if (photon.direction == 1 && newSurf == nsurf - 1) {
                                        if (counter>= (MAX_BOUNCE - 1)) goto maxbounce;
                                        if (transmissionPreCheck(natmospherefile*2 + 1 + 2*nsurf, waveIndex)) {
                                            if (k > 0) goto redoghlp;
                                            countBad_dt(&counterLog, sourceOver*backSourceOver, &ray);
                                            goto photonloop;
                                        }
                                        goto maxbounce;
                                    }

                                    if (photon.direction == -1 && newSurf == 0) {
                                        if (k > 0) goto redoghlp;
                                        countBad_dt(&counterLog, sourceOver*backSourceOver, &ray);
                                        goto photonloop;
                                    }
                                }
                            maxbounce:;
                                maxcounter = counter;
                                if (throughputfile)    {
                                    for (long j = 0; j < nsurf; j++) addThroughput(&throughputLog, j,
                                                                                   waveIndex, -sourceOver*backSourceOver);
                                }

                            }


                        redodiff:;
                            vectorInit(&largeAngle);

                            //  Sample Pupil
                            miss = samplePupil(&positionDiffraction, ray);
                            if (miss) {
                                countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                goto photonloop;
                            }
                            miss = samplePupil(&position, ray);
                            if (miss) {
                                countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                goto photonloop;
                            }
                            if (finiteDistance != 0.0) {
                                double finiteDistanceR = sqrt(pow(finiteDistance, 2) + pow(position.x, 2) +
                                                              pow(position.y, 2));
                                angle.x += position.x/finiteDistanceR;
                                angle.y += position.y/finiteDistanceR;
                                angle.z = smallAnglePupilNormalize(angle.x, angle.y);
                            }
                            if (diffraction_on == 5) vectorCopy(position, &positionDiffraction);
                            xp = position.x;
                            yp = position.y;
                            if (initFlag == 0) {
                                shiftedAngle = spiderangle + time*rotationrate*ARCSEC;
                                wavelengthFactor = pow(wavelength, -0.2)/screen.wavelengthfactor_nom;
                                vectorInit(&largeAngle);
                                initFlag = 1;
                            }

                            //  Diffraction
                            if (diffraction_on >= 1) {
                                miss = diffraction(&positionDiffraction, angle, &largeAngle);
                                if (miss) {
                                    if (k > 0) goto redodiff;
                                    countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                    goto photonloop;
                                }
                            }

                            //   Large Angle Scattering
                            largeAngleScattering(&largeAngle);

                            //   Second Kick
                            if (diffraction_on == 1 && sources.type[source] !=  0) secondKick(&largeAngle);

                            // Saturation Optimization Calculation
                            if (sourceSaturationRadius > 0.0) {
                                if (modulus(&largeAngle)/DEGREE*platescale/pixsize
                                    > sourceSaturationRadius || preGhost>= 2) {
                                    sourceOver_m = 1;
                                    sourceSaturation -= ((1.0 - np)/np + 1.0)*0.001;
                                    if (sourceSaturation<1) sourceSaturation = 1;
                                    break;
                                }
                            }

                        }
                        if (sourceSaturationRadius > 0.0) sourceSaturation += 0.001;
                        if (np<= 0.0) sourceSaturation = 1;
                        counter = -1;
                        if (opdfile) op = 0.0;

                        // Get Delta Angle
                        getDeltaAngle(&angle, &position, source);
                        airRefraction = airIndexRefraction();
                        ncurr = 1.0 + airRefraction/1e6;

                        // ATMOSPHERE

                        // Atmospheric Dispersion
                        if (sources.type[source] !=  0) atmosphericDispersion(&angle);

                        // Loop over Atmosphere Layers
                        for (int layer =-1; layer < natmospherefile; layer++) {

                            // Atmosphere Intercept
                            atmosphereIntercept(&position, angle, layer, diffraction_on);

                            if (sources.type[source] !=  0) {
                                if (layer>= 0) {
                                    // Atmosphere Refraction
                                    atmosphereRefraction(&angle, layer, diffraction_on);

                                    // Clouds
                                    transmission = cloudOpacity(layer);
                                    if (transmissionCheck(transmission, 1 + layer*2, waveIndex)) {
                                        countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                        goto photonloop;
                                    }
                                }

                                // Atmosphere Opacity
                                transmission = atmosphereOpacity(angle, layer);
                                if (transmissionCheck(transmission, 2 + layer*2, waveIndex)) {
                                    countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                    goto photonloop;
                                }

                                if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, layer + 100);

                            }
                        }

                        // Atmosphere Diffraction
                        if (diffraction_on == 2 && sources.type[source] !=  0) atmosphereDiffraction(&angle);

                        // Dome Seeing
                        if (domeseeing > 0.0 || toypsf > 0.0) domeSeeing(&angle);

                        // Tracking
                        if (tracking_on) tracking(&angle, time);

                        // Large Angle
                        angle.x += largeAngle.x;
                        angle.y += largeAngle.y;
                        angle.z = smallAnglePupilNormalize(angle.x, angle.y);

                        if (telescope_on == 0) {
                            newSurf = nsurf - 2;
                            if (nmirror % 2 == 0) {
                                propagate(&position, angle, ((surface.height[nsurf-1] + platescale/DEGREE/1000) - position.z)/angle.z);
                                angle.x -= position.x/(platescale/DEGREE/1000);
                                angle.y -= position.y/(platescale/DEGREE/1000);
                                angle.z = -1.0;
                                normalize(&angle);
                            } else {
                                propagate(&position, angle, ((surface.height[nsurf-1] - platescale/DEGREE/1000) - position.z)/angle.z);
                                angle.x -= position.x/(platescale/DEGREE/1000);
                                angle.y -= position.y/(platescale/DEGREE/1000);
                                angle.z = 1.0;
                                normalize(&angle);
                            }
                            if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, 200);
                        } else {
                            newSurf = -1;
                        }

                        // OPTICS AND DETECTOR
                        photon.direction = 1;
                        ghostFlag = 0;
                    surfaceloop: while (1) {
                            oldSurf = newSurf;
                            if (photon.direction == 1) newSurf++; else newSurf--;

                        redostraylight:;

                            // Find intercept
                            if (newSurf>= 0 && newSurf<nsurf) {

                                transform(&position, &angle, newSurf);
                                if (surface.surfacetype[newSurf] == DETECTOR) {
                                    transform(&position, &angle, newSurf + 1);
                                }
                                miss = findSurface(angle, position, &distance, newSurf);

                            } else {
                                miss = 1;
                            }

                            //   Missed surface or wrong direction
                            if (miss == 1 || (distance < 0)) {
                                if (straylightcurrent == 0) {
                                    countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                    goto photonloop;
                                } else {
                                    if (chooseSurface(&newSurf, &oldSurf) == 1) {
                                        countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                        goto photonloop;
                                    } else {
                                        goto redostraylight;
                                    }
                                }
                            }

                            propagate(&position, angle, distance);


                            if (opdfile) {
                                op += distance*ncurr;

                                if (newSurf == 0) {
                                    opdx = position.x;
                                    opdy = position.y;
                                }
                                if (surface.surfacetype[newSurf] == DETECTOR) {
                                    distance = sqrt(pow(xporig - position.x, 2) + pow(yporig - position.y, 2) +
                                                    pow(zporig - position.z, 2));
                                    distance += sqrt(pow(70*1e6 - position.z, 2) - xp*xp - yp*yp) - (70*1e6 - position.z);
                                    op -= ncurr*distance;
                                }
                            }

                            if (throughputfile == 1 && photon.direction == 1) {
                                addThroughput(&throughputLog, newSurf, waveIndex, sourceOver*backSourceOver);
                            }

                            if (photon.direction == -1) ghostFlag = 1;

                            //   DERIVATIVES
                            interceptDerivatives(&normal, position, newSurf);

                            //   CONTAMINATION
                            if (surface.surfacetype[newSurf] !=  DETECTOR && contaminationmode==1) {
                                miss = contaminationSurfaceCheck(position, &angle, newSurf);
                                if (miss) {
                                    countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                    goto photonloop;
                                }
                            }

                            //   SURFACE COATINGS
                            transmission = surfaceCoating(wavelength, angle, normal, newSurf, &reflection);

                            if (transmissionCheck(transmission, natmospherefile*2 + 1 + newSurf*2, waveIndex)) {
                                if (straylightcurrent == 1 && ghost[newSurf] == 0) {
                                    if (transmissionCheck(reflection + transmission, natmospherefile*2 + 1 + newSurf*2 + 1, waveIndex)) {
                                        countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                        goto photonloop;
                                    } else {
                                        photon.direction = -photon.direction;
                                        reflect(&angle, normal);
                                        transformInverse(&position, &angle, newSurf);
                                        if (surface.surfacetype[newSurf] == DETECTOR) {
                                            transformInverse(&position, &angle, newSurf + 1);
                                        }
                                        goto surfaceloop;
                                    }
                                } else {
                                    countBad(&counterLog, sourceOver*backSourceOver, &ray);
                                    goto photonloop;
                                }
                            }


                            //   INTERACTIONS
                            if (surface.surfacetype[newSurf] == MIRROR) {

                                //   MIRROR
                                reflect(&angle, normal);
                                transformInverse(&position, &angle, newSurf);
                                if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);

                            } else if (surface.surfacetype[newSurf] == LENS || surface.surfacetype[newSurf] == FILTER) {

                                //   LENS/FILTER
                                newRefractionIndex(newSurf);
                                refract(&angle, normal, nprev, ncurr);
                                transformInverse(&position, &angle, newSurf);
                                if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);

                            } else if (surface.surfacetype[newSurf] == GRATING) {

                                //   GRATING
                                double wavelengthNm = wavelength*1000.0;
                                Vector angleOut;
                                pGrating->diffract(angle.x, angle.y, angle.z, normal.x, normal.y, normal.z,
                                                   angleOut.x, angleOut.y, angleOut.z, wavelengthNm);
                                vectorCopy(angleOut, &angle);
                                transformInverse(&position, &angle, newSurf);
                                if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);

                            } else if (surface.surfacetype[newSurf] == DETECTOR) {

                                if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, newSurf + 200);
                                vectorCopy(position, &positionPrevious);
                                vectorCopy(angle, &anglePrevious);
                                detRay = 0;

                            detectorloop: while (detRay<backSourceOver) {


                                    position.x = positionPrevious.x + random_gaussian()*backSigma;
                                    position.y = positionPrevious.y + random_gaussian()*backSigma;
                                    position.z = positionPrevious.z;
                                    vectorCopy(anglePrevious, &angle);

                                    //   SILICON
                                    if (detector_on) {

                                        double dh;
                                        miss = getDeltaIntercept(position.x, position.y, &dh, newSurf);

                                        xPos = static_cast<long>(floor(position.x*1000/pixsize + pixelsx/2));
                                        yPos = static_cast<long>(floor(position.y*1000/pixsize + pixelsy/2));
                                        xPosR = position.x*1000/pixsize-floor(position.x*1000/pixsize);
                                        yPosR = position.y*1000/pixsize-floor(position.y*1000/pixsize);

                                        if (xPos >= minx && xPos <= maxx && yPos >= miny && yPos <= maxy && contaminationmode == 1) {
                                            if (RngDouble()>(double)(*(contamination.chiptransmission +
                                                                       chip.nampx*(yPos - miny) + (xPos - minx)))) {
                                                countBad(&counterLog, sourceOver, &ray);
                                                detRay++;
                                                goto detectorloop;
                                            }
                                            if (*(contamination.chiplistmap + chip.nampx*(yPos - miny) + (xPos - minx)) != -1) {
                                                double xx = ((position.x*1000/pixsize) + pixelsx/2)*1e-3*pixsize;
                                                double yy = ((position.y*1000/pixsize) + pixelsy/2)*1e-3*pixsize;
                                                long cc = *(contamination.chiplistmap + chip.nampx*(yPos - miny) + (xPos - minx));
                                                if (sqrt(pow(xx - contamination.chiplistx[cc], 2.0) +
                                                         pow(yy - contamination.chiplisty[cc], 2.0)) <
                                                    contamination.chiplists[cc]) {
                                                    long index;
                                                    find(contamination.henyey_greenstein, contamination.elements, RngDouble(), &index);
                                                    double mu = contamination.henyey_greenstein_mu[index];
                                                    double phi = RngDouble()*2*M_PI;
                                                    shift_mu(&angle, mu, phi);
                                                    if (mu<0) {
                                                        countBad(&counterLog, sourceOver, &ray);
                                                        detRay++;
                                                        goto detectorloop;
                                                    }
                                                }
                                            }
                                        }


                                        miss = siliconPropagate(&angle, &position, wavelength, normal, &interact, &collect, dh, waveIndex);

                                        if (miss == 0) {

                                            vectorCopy(interact, &position);

                                            if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, 300);
                                            vectorCopy(collect, &position);

                                            if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, 301);

                                        } else {

                                            vectorCopy(interact, &position);
                                            if (eventfile) pEventLogging->logPhoton(position.x, position.y, position.z, 302);
                                            countBad(&counterLog, sourceOver, &ray);
                                            detRay++;
                                            goto detectorloop;

                                        }


                                    }


                                    xPos = (long)(floor(position.x*1000/pixsize + pixelsx/2));
                                    yPos = (long)(floor(position.y*1000/pixsize + pixelsy/2));

                                    if (eventfile) {
                                        if (xPos >= minx && xPos <= maxx && yPos >= miny && yPos <= maxy) {
                                            pEventLogging->logPhoton((double)xPos, (double)yPos, 0.0, 303);
                                        }
                                    }

                                    if (centroidfile) {
                                        if (xPos >= minx && xPos <= maxx && yPos >= miny && yPos <= maxy) {
                                            source_xpos[source] += xPos*sourceOver;
                                            source_ypos[source] += yPos*sourceOver;
                                            source_photon[source] += sourceOver;
                                        }
                                    }


                                    if (opdfile) {
                                        if (xPos >= minx && xPos <= maxx && yPos >= miny && yPos <= maxy) {
                                            long xx = opdx/(screen.fine_sizeperpixel*4) + SCREEN_SIZE/8;
                                            long yy = opdy/(screen.fine_sizeperpixel*4) + SCREEN_SIZE/8;
                                            if (xx >= 0 && xx < SCREEN_SIZE/4 && yy >= 0 && yy < SCREEN_SIZE/4) {
                                                *(opd + SCREEN_SIZE/4*xx + yy) += op;
                                                *(opdcount + SCREEN_SIZE/4*xx + yy) += 1;
                                            }
                                        }
                                    }

                                    if (sources.type[source]<3 && backGamma>1.0) {
                                        Vector newPosition;
                                        vectorCopy(position, &newPosition);
                                        long long lmax;
                                        lmax = sourceOver_m;
                                        sourceOver_m = 1;
                                        for (long long l = 0; l < lmax; l++) {
                                            position.x = newPosition.x + random_gaussian()*backSigma/backDelta;
                                            position.y = newPosition.y + random_gaussian()*backSigma/backDelta;
                                            xPos = static_cast<long>(floor(position.x*1000/pixsize + pixelsx/2));
                                            yPos = static_cast<long>(floor(position.y*1000/pixsize + pixelsy/2));
                                            if (xPos >= minx && xPos <= maxx && yPos >= miny && yPos <= maxy) {
                                                if (saturation) {
                                                    saturate(source, &largeAngle);
                                                } else {
                                                    *(chip.focal_plane + chip.nampx*(yPos - miny) + (xPos - minx)) += sourceOver_m;
                                                }
                                                countGood(&counterLog, sourceOver_m, &ray);
                                            } else {
                                                countBad(&counterLog, sourceOver_m, &ray);
                                            }
                                        }
                                        sourceOver_m = lmax;
                                    } else {
                                        if (xPos >= minx && xPos <= maxx && yPos >= miny && yPos <= maxy) {
                                            if (saturation) {
                                                saturate(source, &largeAngle);
                                            } else {
                                                *(chip.focal_plane + chip.nampx*(yPos - miny) + (xPos - minx)) += sourceOver_m;
                                            }
                                        } else {
                                            countBad(&counterLog, sourceOver, &ray);
                                            detRay++;
                                            goto detectorloop;
                                        }
                                    }


                                    if (throughputfile) addThroughput(&throughputLog, nsurf, waveIndex, sourceOver_m);
                                    detRay++;
                                    if ((sources.type[source] < 3 && backGamma > 1.0) || sources.type[source]>= 3
                                        || (sources.type[source] < 3 && backGamma <= 1.0)) countGood(&counterLog, sourceOver_m, &ray);
                                }
                                break;

                            }
                        }
                    }
                    sourceCounter++;
                }
            }
        }

        if (sourceType >= 0  && sourceType < 10) sprintf(tempstring, "Dome Light         ");
        if (sourceType >= 10 && sourceType < 20) sprintf(tempstring, "Dark Sky           ");
        if (sourceType >= 20 && sourceType < 30) sprintf(tempstring, "Moon               ");
        if (sourceType == 30) sprintf(tempstring, "Astrophysical m>40 ");
        if (sourceType >= 31 && sourceType < 70) sprintf(tempstring, "Astrophysical m=%2d ", 71 - sourceType);
        if (sourceType == 72) sprintf(tempstring, "Astrophysical m<0  ");
        if (sourceCounter > 0) counterCheck(&counterLog, sourceCounter, tempstring);
    }

    // COSMIC RAYS
    if (checkpointcount == checkpointtotal) {
        detRay = 0; ray = 0;
        cosmicRays(&detRay);
        if (detRay > 0) {
            for (long i = 0; i < detRay; i++) countGood(&counterLog, 1, &ray);
            sprintf(tempstring, "Cosmic Rays        ");
            counterCheck(&counterLog, detRay, tempstring);
        }
    }

    // OUTPUT DATA
    if (checkpointcount == checkpointtotal) writeImageFile();
    if (opdfile) writeOPD();
    if (checkpointcount !=  checkpointtotal) writeCheckpoint(checkpointcount);
    if (centroidfile) writeCentroidFile(outputdir, outputfilename, source_photon, source_xpos, source_ypos, sources.id, nsource);
    if (throughputfile) writeThroughputFile(outputdir, outputfilename, &throughputLog, nsurf);
    if (eventfile) pEventLogging->eventFileClose();
    return(0);

}
