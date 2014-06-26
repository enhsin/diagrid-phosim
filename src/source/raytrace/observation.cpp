///
/// @package phosim
/// @file observation.cpp
/// @brief observation class for raytrace code
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author Alan Meert (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <iostream>
#include <sstream>
#include <fstream>
#include <zlib.h>
#include <string.h>
#include <math.h>
#include <fitsio.h>
#include <fitsio2.h>

#include "raytrace.h"
#include "observation.h"
#include "constants.h"
#include "helpers.h"
#include "constants.h"
#include "ancillary/readtext.h"

using readtext::readText;

int Observation::addSource (const std::string & object, int sourcetype) {

    char tempstring[4096];
    int nspatialpar, ndustpar, ndustparz;
    long i, j;
    char tempstring1[512];
    char tempstring2[512];
    double tempf1, tempf2;
    long lsedptr;
    FILE *indafile;
    gzFile ingzfile;
    char line[4096];
    int status;
    char *ffptr;
    fitsfile *faptr;
    long naxes[2];
    int nfound;
    int anynull;
    float nullval;
    double cdw, dw;
    long closestw = 0;
    double x, y;
    double nn;
    double mag;
    long badfile, oldsed;
    char *sptr, *sptr2, *sptr3;
    double ra, dec, id, redshift, gamma1, gamma2, kappa, magnification,
        deltara, deltadec;

    nspatialpar = 0;

    std::istringstream iss(object);
    iss >> id >> ra >> dec;
    sources.id.push_back(id);
    sources.ra.push_back(ra);
    sources.dec.push_back(dec);
    sources.ra[nsource] *= M_PI/180.0;
    sources.dec[nsource] *= M_PI/180.0;
    iss >> mag >> sources.sedfilename[nsource] >>redshift >> gamma1 >> gamma2 >> kappa >> deltara >> deltadec;
    sources.redshift.push_back(redshift);
    sources.gamma1.push_back(gamma1);
    sources.gamma2.push_back(gamma2);
    sources.kappa.push_back(kappa);
    sources.deltara.push_back(deltara);
    sources.deltadec.push_back(deltadec);

    sources.type[nsource] = sourcetype;
    magnification = 2.5*log10((1 - sources.kappa[nsource])*(1 - sources.kappa[nsource])-
                              sources.gamma1[nsource]*sources.gamma1[nsource]-
                              sources.gamma2[nsource]*sources.gamma2[nsource]);
    sources.norm[nsource] = pow(10.0, ((mag + magnification + 48.6)/(-2.5)));
    sources.mag[nsource] = mag + magnification;

    oldsed = 0;
    // read SED file
    if (nsource > 0) {
        for (i = 0; i < nsource; i++) {
            if (sources.sedfilename[i] == sources.sedfilename[nsource]) {
                sources.sedptr[nsource] = sources.sedptr[i];
                oldsed = 1;
                goto skipsedread;
            }
        }
    }

    sources.sedptr[nsource] = nsedptr;

    if (sedptr == 0) {
        nreallocs = 0;
        sed_max = 4194304; // 4M elements, 32MB memory
        sed_w = (double*)malloc((long)(sed_max*sizeof(double)));
        sed_c = (double*)malloc((long)(sed_max*sizeof(double)));
    } else {
        if (sedptr > (sed_max - 25000)) {
            ++nreallocs;
            sed_max = 2*sed_max;
            sed_w = (double*)realloc(sed_w, (long)((sed_max)*sizeof(double)));
            sed_c = (double*)realloc(sed_c, (long)((sed_max)*sizeof(double)));
        }
    }

    lsedptr = 0;
    sprintf(tempstring, "%s/%s", seddir.c_str(), sources.sedfilename[nsource].c_str());
    if (flatdir == 1) {
        sptr = strtok_r(tempstring, "/", &sptr2);
        do {
            sptr3 = sptr;
            sptr = strtok_r(NULL, "/", &sptr2);
        } while (sptr != NULL);
        sprintf(tempstring, "%s", sptr3);
    }

    if (strstr(tempstring, ".gz") == NULL) {

        indafile = fopen(tempstring, "r");
        if (indafile == NULL) {
            fprintf(stderr, "Can't find SED file: %s\n", tempstring);
            exit(1);
        }

        closestw = 0;
        cdw = 1e30;
        badfile = 0;
        while (fgets(line, 4096, indafile)) {
            sscanf(line, "%s %s", tempstring1, tempstring2);
            tempf1 = strtod(tempstring1, NULL);
            tempf2 = strtod(tempstring2, NULL);
            *(sed_w + sedptr + lsedptr) = tempf1;
            *(sed_c + sedptr + lsedptr) = tempf2*tempf1*1e-7/(1.98645e-16);
            if (tempstring2[0] == 'n' || tempstring2[0] == 'N') {
                if (badfile == 0) {
                    fprintf(stderr, "Error:   SED file: %s contains a NaN!\n", tempstring);
                }
                *(sed_c + sedptr + lsedptr) = 1e-20;
                badfile = 1;
            }

            dw = fabs(tempf1 - 500.0);

            if (dw < cdw) {
                cdw = dw;
                closestw = lsedptr;
            }
            lsedptr = lsedptr + 1;
            if (lsedptr >= 25000) {
                fprintf(stderr, "Error:  Too many lines in SED file: %s\n", tempstring);
                exit(1);
            }
        }
        fclose(indafile);

    } else {

        ingzfile = gzopen(tempstring, "r");
        if (ingzfile == NULL) {
            fprintf(stderr, "Can't find SED file: %s\n", tempstring);
            exit(1);
        }

        closestw = 0;
        cdw = 1e30;
        badfile = 0;
        while (gzgets(ingzfile, line, 4096)) {
            sscanf(line, "%s %s", tempstring1, tempstring2);
            tempf1 = strtod(tempstring1, NULL);
            tempf2 = strtod(tempstring2, NULL);

            *(sed_w + sedptr + lsedptr) = tempf1;
            *(sed_c + sedptr + lsedptr) = tempf2*tempf1*1e-7/(1.98645e-16);
            if (tempstring2[0] == 'n' || tempstring2[0] == 'N') {
                if (badfile == 0) fprintf(stderr, "Error:   SED file: %s contains a NaN!\n", tempstring);
                *(sed_c + sedptr + lsedptr) = 1e-20;
                badfile = 1;
            }

            dw = fabs(tempf1 - 500.0);

            if (dw < cdw) {
                cdw = dw;
                closestw = lsedptr;
            }
            lsedptr = lsedptr + 1;
            if (lsedptr >= 25000) {
                fprintf(stderr, "Error:  Too many lines in SED file: %s\n", tempstring);
                exit(1);
            }
        }
        gzclose(ingzfile);

    }


    for (j = 0; j < lsedptr; j++) {
        if (j != 0 && j != (lsedptr - 1)) *(sed_c + sedptr + j) = *(sed_c + sedptr + j)*(*(sed_w + sedptr + j + 1) - *(sed_w + sedptr + j - 1))/2.0;
        if (j == 0) *(sed_c + sedptr + j) = *(sed_c + sedptr + j)*(*(sed_w + sedptr + j + 1) - *(sed_w + sedptr + j));
        if (j == (lsedptr - 1)) *(sed_c + sedptr + j) = *(sed_c + sedptr + j)*(*(sed_w + sedptr + j) - *(sed_w + sedptr + j - 1));
    }

    tempf1 = 0;
    for (j = 0; j < lsedptr; j++) {
        tempf1 += *(sed_c + sedptr + j);
    }
    for (j = 0; j < lsedptr; j++) {
        *(sed_c + sedptr + j)=*(sed_c + sedptr + j)/tempf1;
    }

    if (closestw == 0) {
        sed_dwdp[nsedptr] = (*(sed_w + sedptr + closestw + 1) - *(sed_w + sedptr + closestw))/(*(sed_c + sedptr + closestw))/1.0;
    } else {
        sed_dwdp[nsedptr] = (*(sed_w + sedptr + closestw + 1) - *(sed_w + sedptr + closestw - 1))/(*(sed_c + sedptr + closestw))/2.0;
    }
    if (*(sed_c + sedptr + closestw) <= 0.0) {
        printf("Error in SED file; 0 value at 500 nm\n");
        sed_dwdp[nsedptr] = 0.0;
    }

    for (j = 1; j < lsedptr; j++) {
        *(sed_c + sedptr + j) += *(sed_c + sedptr + j - 1);
    }

    if (oldsed == 0) {
        sed_n[nsedptr] = lsedptr;
        sed_ptr[nsedptr] = sedptr;
        sedptr = sedptr + lsedptr;
        nsedptr++;
    }
    if (nsedptr >= 10000) {
        printf("Error:   Too many SED files\n");
        exit(1);
    }

 skipsedread:;
    sources.norm[nsource] = sources.norm[nsource]/(500.0)*(1 + sources.redshift[nsource])*sed_dwdp[sources.sedptr[nsource]];

    iss >> sources.spatialname[nsource];

    if (sources.spatialname[nsource] == "point") {
        sources.spatialtype[nsource] = POINT;
        nspatialpar = 0;
    }
    if (sources.spatialname[nsource].find("fit") != std::string::npos) {
        sources.spatialtype[nsource] = IMAGE;
        nspatialpar = 2;
    }
    if (sources.spatialname[nsource] == "gauss") {
        sources.spatialtype[nsource] = GAUSSIAN;
        nspatialpar = 1;
    }
    if (sources.spatialname[nsource] == "sersic") {
        sources.spatialtype[nsource] = SERSIC;
        nspatialpar = 6;
    }
    if (sources.spatialname[nsource] == "sersic2d") {
        sources.spatialtype[nsource] = SERSIC2D;
        nspatialpar = 4;
    }
    if (sources.spatialname[nsource] == "sersic2D") {
        sources.spatialtype[nsource] = SERSIC2D;
        nspatialpar = 4;
    }
    if (sources.spatialname[nsource] == "movingpoint") {
        sources.spatialtype[nsource] = MOVINGPOINT;
        nspatialpar = 2;
    }
    if (sources.spatialname[nsource] == "pinhole") {
        sources.spatialtype[nsource] = PINHOLE;
        nspatialpar = 4;
    }

    if (sources.spatialtype[nsource] == 1) {
        if (nsource > 0) {
            sources.skysameas[nsource] = -1;
            for (i = 0; i < nsource; i++) {
                if (sources.spatialname[nsource] == sources.spatialname[i]) sources.skysameas[nsource] = i;
            }
        } else {
            sources.skysameas[nsource] = -1;
        }


        /* read image file */

        if (sources.skysameas[nsource] == -1) {
            sprintf(tempstring, "%s/%s+0", imagedir.c_str(), sources.spatialname[nsource].c_str());

            ffptr = tempstring;
            status = 0;

            if (fits_open_file(&faptr, ffptr, READONLY, &status)) printf("Error opening %s\n", ffptr);
            fits_read_keys_lng(faptr, (char*)"NAXIS", 1, 2, naxes, &nfound, &status);
            if ((tempptr[nimage] = (float*)malloc(naxes[0]*naxes[1]*sizeof(float))) == NULL) {
                fprintf(stderr, "Can't allocate image.\n");
                exit(1);
            }
            naxesb[nimage][0] = naxes[1];
            naxesb[nimage][1] = naxes[0];
            fits_read_img(faptr, TFLOAT, 1, naxes[0]*naxes[1], &nullval,
                          tempptr[nimage], &anynull, &status);
            fits_close_file(faptr, &status);


            cumulative[nimage] = 0;
            cumulativex[nimage] = (float*)malloc(naxesb[nimage][0]*sizeof(float));
            for (i = 0; i < naxesb[nimage][0]; i++) {
                cumulativex[nimage][i] = cumulative[nimage];
                for (j = 0; j < naxesb[nimage][1]; j++) {
                    if (*(tempptr[nimage] + i*naxesb[nimage][1] + j) < 0) {
                        *(tempptr[nimage] + i*naxesb[nimage][1] + j) = 0;
                    }
                    cumulative[nimage] += *(tempptr[nimage] + i*naxesb[nimage][1] + j);
                }
            }

            sources.spatialpar[nsource][2] = nimage;
            nimage++;
        }

    }



    if (nspatialpar > 0) {
        for (i = 0; i < nspatialpar; i++) {
            iss >> sources.spatialpar[nsource][i];
        }
    }

    iss >> sources.dustnamez[nsource];

    sources.dusttypez[nsource] = 0;
    ndustparz = 0;
    if (sources.dustnamez[nsource] == "ccm") {
        sources.dusttypez[nsource] = 1;
        ndustparz = 2;
    } else if (sources.dustnamez[nsource] == "calzetti") {
        sources.dusttypez[nsource] = 2;
        ndustparz = 2;
    } else if (sources.dustnamez[nsource] == "CCM") {
        sources.dusttypez[nsource] = 1;
        ndustparz = 2;
    } else if (sources.dustnamez[nsource] == "CALZETTI") {
        sources.dusttypez[nsource] = 2;
        ndustparz = 2;
    }

    if (ndustparz > 0) {
        for (i = 0; i < ndustparz; i++) {
            iss >> sources.dustparz[nsource][i];
        }
    }


    iss >> sources.dustname[nsource];

    sources.dusttype[nsource] = 0; ndustpar = 0;
    if (sources.dustname[nsource] == "ccm") {sources.dusttype[nsource] = 1; ndustpar = 2;}
    else if (sources.dustname[nsource] == "calzetti") {sources.dusttype[nsource] = 2; ndustpar = 2;}
    else if (sources.dustname[nsource] == "CCM") {sources.dusttype[nsource] = 1; ndustpar = 2;}
    else if (sources.dustname[nsource] == "CALZETTI") {sources.dusttype[nsource] = 2; ndustpar = 2;}

    if (ndustpar > 0) {
        for (i = 0; i < ndustpar; i++) {
            iss >> sources.dustpar[nsource][i];
        }
    }

    setup_tangent(pra, pdec, &tpx, &tpy, &tpz);

    tangent(sources.ra[nsource] + sources.deltara[nsource], sources.dec[nsource] + sources.deltadec[nsource], &x, &y, &tpx, &tpy, &tpz);

    sources.vx[nsource] = x*cos(rotationangle) - y*sin(rotationangle);
    sources.vy[nsource] = x*sin(rotationangle) + y*cos(rotationangle);
    sources.vz[nsource] = -1.0;
    nn = sqrt((sources.vx[nsource])*(sources.vx[nsource]) +
              (sources.vy[nsource])*(sources.vy[nsource]) + 1);
    sources.vx[nsource] = sources.vx[nsource]/nn;
    sources.vy[nsource] = sources.vy[nsource]/nn;
    sources.vz[nsource] = sources.vz[nsource]/nn;

    nsource++;
    return(0);
}


int Observation::background () {

    char tempstring[4096];
    double focal_length = platescale*180.0/M_PI;
    double A, F, xp, yp, ra, dec, xv, yv, y_max_distance, currbuffer, source_dist_x, source_dist_y, x_max_distance;
    double dx , dy;
    long ndeci,  nrai, deci, rai;
    long over, i, j;
    double dra, dis, cosdis;

    int nspatialpar, ii, jj;
    char tempstring1[512];
    char tempstring2[512];
    double tempf1,  tempf2;
    long lsedptr;
    FILE *indafile;
    gzFile ingzfile;
    char line[4096];
    std::string dir;
    double cdw, dw;
    long closestw = 0;
    double x, y;
    double nn;
    double mag = 100;
    long badfile, oldsed;
    char *sptr, *sptr2, *sptr3;

    int diffusetype;
    double angular_sep_degrees, angular_sep_radians, moon_magnitude, moon_apparent_magnitude,
        scatter_function, moon_illuminance, lunar_illuminance, darksky_magnitude;

    if (flatdir == 0) dir=datadir + "/sky";
    else if (flatdir == 1) dir = ".";

    airglow = (float*)calloc((airglowScreenSize)*(airglowScreenSize), sizeof(float));
    sprintf(tempstring, "airglowscreen_%ld.fits", obshistid);
    {
        char *ffptr;
        fitsfile *faptr;
        long naxes[2];
        int nfound;
        int anynull;
        float nullval;
        int status;

        ffptr = tempstring;
        status = 0;
        if (fits_open_file(&faptr, ffptr, READONLY, &status)) {printf("Error opening %s\n", ffptr); exit(1);}
        fits_read_keys_lng(faptr, (char*)"NAXIS", 1, 2, naxes, &nfound, &status);
        fits_read_img(faptr, TFLOAT, 1, naxes[0]*naxes[1], &nullval, airglow, &anynull, &status);
        fits_close_file(faptr, &status);
    }

    //CALCULATIONS FOR BACKGROUND
    moon_apparent_magnitude = -12.73 + .026 * (fabs(phaseang * 180 / M_PI)) + (4E-9) * pow(phaseang * 180 / M_PI,  4);
    moon_illuminance = pow(10, -0.4 * (moon_apparent_magnitude + 16.57));
    lunar_illuminance = 26.3311157 - 1.08572918 * log(moon_illuminance);

    //TWILIGHT
    float temp;
    float background_magnitude = zenith_v;
    float background_brightness = 34.08 * exp( 20.72333 - .92104 * background_magnitude );
    float darksky[6];
    float darksky_data[2511][2];
    char darksky_sedfile[4096];
    FILE *darksky_fp;

    sprintf( darksky_sedfile, "%s/darksky_sed.txt",dir.c_str());

    darksky_fp = fopen( darksky_sedfile, "r" );
    for ( i = 0; i < 2511; i++ ){
        fgets( line, 4096, darksky_fp );
        sscanf( line, "%f %f\n", &darksky_data[i][0], &temp );
        darksky_data[i][1] = temp * background_brightness / 6.299537E-18;
    }
    fclose( darksky_fp );
    darksky[0] = darksky_data[39][1];
    darksky[1] = darksky_data[491][1];
    darksky[2] = darksky_data[1019][1];
    darksky[3] = darksky_data[1547][1];
    darksky[4] = darksky_data[1961][1];
    darksky[5] = darksky_data[2250][1];

    float lunar[6];
    float lunar_data[7500][2];
    char lunar_sedfile[4096];
    FILE *lunar_fp;

    sprintf( lunar_sedfile, "%s/lunar_sed.txt",dir.c_str());

    lunar_fp = fopen ( lunar_sedfile, "r" );
    for ( i = 0; i < 1441; i++ ){
        fgets( line, 200, lunar_fp );
        sscanf( line, "%f %f\n", &lunar_data[i][0], &temp );
        lunar_data[i][1] = temp * background_brightness / 3.882815E-16;
    }
    fclose( lunar_fp );

    lunar[0] = lunar_data[600][1];
    lunar[1] = lunar_data[1800][1];
    lunar[2] = lunar_data[3200][1];
    lunar[3] = lunar_data[4600][1];
    lunar[4] = lunar_data[5700][1];
    lunar[5] = lunar_data[6700][1];

    float a[6][3] = {{ 11.78, 1.376, -.039 }, { 11.84, 1.411, -.041 }, { 11.84, 1.518, -.057 }, \
                     { 11.40, 1.567, -.064 }, { 10.93, 1.470, -.062 }, { 10.43, 1.420, -.052 }};
    float color[6] = { 0.67, 1.03, 0, -0.74, -1.90, -2.20 };

    float magnitude = 100.0, brightness = 0.0, angle = solarzen / M_PI * 180.0;
    j = filter;
    float alpha = 0.0, beta = 0.0;

    if ( angle <= 106 ) {
        alpha = 1.0 - ( angle - 95.0 ) / 11.0;
        beta = 1.0 - alpha;
        magnitude = a[j][0] + a[j][1] * ( angle - 95.0 ) + a[j][2] * ( angle - 95.0 ) * ( angle - 95.0 );
    }
    else if ( angle >= 254 ){
        alpha = 1.0 - ( 265 - angle ) / 11.0;
        beta = 1.0 - alpha;
        magnitude = a[j][0] + a[j][1] * ( 265.0 - angle ) + a[j][2] * ( 265.0 - angle ) * ( 265.0 - angle );
    }
    else if ( (angle > 106) && (angle < 254) ){
        alpha = 0.0;
        beta = 1.0;
        magnitude = a[j][0] + a[j][1] * ( 11.0 ) + a[j][2] * ( 11.0 ) * ( 11.0 );
    }

    brightness = 34.08 * exp( 20.72333 - .92104 * ( magnitude - color[j] ) );

    float lunar_like_magnitude, darksky_like_magnitude;
    float lunar_like_brightness, darksky_like_brightness;
    float lunar_like_photon_count, darksky_like_photon_count;

    lunar_like_brightness = 0.5 * ( brightness - 2.0 * alpha * lunar[j] );
    darksky_like_brightness = 0.5 * ( brightness - 2.0 * beta * darksky[j] );

    if ( lunar_like_brightness < 0 ) lunar_like_brightness = 0.0;
    if ( darksky_like_brightness < 0 ) darksky_like_brightness = 0.0;

    lunar_like_magnitude = 26.33111 - 1.08573 * log( lunar_like_brightness ) + color[j];
    darksky_like_magnitude = 26.33111 - 1.08573 * log( darksky_like_brightness ) + color[j];

    lunar_like_photon_count = expf(-.4 * lunar_like_magnitude * 2.30258509);
    if (lunar_like_photon_count < 0) lunar_like_photon_count=0;

    darksky_like_photon_count = expf(-.4 * darksky_like_magnitude * 2.30258509);
    if (darksky_like_photon_count < 0) darksky_like_photon_count = 0;

    if ( (angle > 106) && (angle < 130) ) darksky_like_photon_count *= exp( 1 - 24.0 / fabs( angle - 130.0 ) );
    if ( (angle > 230) && (angle < 254) ) darksky_like_photon_count *= exp( 1 - 24.0 / fabs( angle - 230.0 ) );
    if ( (angle >= 130) && (angle <= 230) ) darksky_like_photon_count = 0;

    if ( (angle > 106) && (angle < 130) ) lunar_like_photon_count *= exp( 1 - 24.0 / fabs( angle - 130.0 ) );
    if ( (angle > 230) && (angle < 254) ) lunar_like_photon_count *= exp( 1 - 24.0 / fabs( angle - 230.0 ) );
    if ( (angle >= 130) && (angle <= 230) ) lunar_like_photon_count = 0;


    for (diffusetype = 0; diffusetype < 3; diffusetype++){
        if ((diffusetype == 0 && domelight < 100) || (diffusetype == 1 && zenith_v < 100) || (diffusetype == 2 && moonalt > 0)) {

            //POPULATE SOURCES
            ndeci = (long)(180*3600/backRadius);
            over = (long)(60/backRadius);
            for (deci = 0; deci < ndeci; deci += over) {
                dec = ((deci + 0.5 - ndeci/2)/(ndeci))*M_PI;
                nrai = ((long)(ndeci*cos(dec)*2));
                for (rai = 0; rai < nrai; rai += over) {
                    ra = 2*M_PI*((rai + 0.5 - nrai/2)/(nrai));

                    A = cos(dec)*cos(ra - pra);
                    F = focal_length/(sin(pdec)*sin(dec)+A*cos(pdec));
                    yp = F*(cos(pdec)*sin(dec)-A*sin(pdec));
                    xp = F*cos(dec)*sin(ra-pra);
                    xv = (xp * cos(-rotationangle) + yp * sin(-rotationangle));
                    yv = (-xp * sin(-rotationangle) + yp * cos(-rotationangle));
                    currbuffer = backBuffer+(3*backBeta*backRadius+60.0)*ARCSEC*focal_length/pixsize;
                    x_max_distance = pixelsx * pixsize/2.0 + currbuffer * pixsize;
                    y_max_distance = pixelsy * pixsize/2.0 + currbuffer * pixsize;
                    dx = xv - centerx - decenterx;
                    dy = yv - centery - decentery;
                    source_dist_x = fabs(cos(chipangle)*dx + sin(chipangle)*dy);
                    source_dist_y = fabs(-sin(chipangle)*dx + cos(chipangle)*dy);

                    dra = fabs(ra-pra);
                    if (dra > M_PI) dra = 2*M_PI-dra;
                    cosdis = sin(dec)*sin(pdec)+cos(dec)*cos(pdec)*cos(dra);
                    if (cosdis>1) cosdis = 1.0;
                    if (cosdis<-1) cosdis = -1.0;
                    dis = acos(cosdis);

                    if ( (source_dist_x <= x_max_distance) && (source_dist_y <= y_max_distance) &&
                         (dis<M_PI/2) ){

                        for (i = 0; i < over; i++) {
                            for (j = 0; j < over; j++) {
                                dec = ((deci + i + 0.5- ndeci/2)/(ndeci))*M_PI;
                                ra = 2*M_PI*((rai + j + 0.5 - nrai/2)/(nrai));

                                A = cos(dec)*cos(ra - pra);
                                F = focal_length/(sin(pdec)*sin(dec) + A*cos(pdec));
                                yp = F*(cos(pdec)*sin(dec) - A*sin(pdec));
                                xp = F*cos(dec)*sin(ra - pra);
                                xv = (xp * cos(-rotationangle) + yp * sin(-rotationangle));
                                yv = (-xp * sin(-rotationangle) + yp * cos(-rotationangle));
                                currbuffer = backBuffer+3*backBeta*backRadius*ARCSEC*focal_length/pixsize;
                                x_max_distance = pixelsx * pixsize/2.0 + currbuffer * pixsize;
                                y_max_distance = pixelsy * pixsize/2.0 + currbuffer * pixsize;
                                dx = xv - centerx - decenterx;
                                dy = yv - centery - decentery;
                                source_dist_x = fabs(cos(chipangle)*dx + sin(chipangle)*dy);
                                source_dist_y = fabs(-sin(chipangle)*dx + cos(chipangle)*dy);

                                dra = fabs(ra - pra);
                                if (dra > M_PI) dra = 2*M_PI - dra;
                                cosdis = sin(dec)*sin(pdec) + cos(dec)*cos(pdec)*cos(dra);
                                if (cosdis>1) cosdis = 1.0;
                                if (cosdis<-1) cosdis = -1.0;
                                dis = acos(cosdis);

                                if ( (source_dist_x <= x_max_distance) && (source_dist_y <= y_max_distance) &&
                                     (dis<M_PI/2)){

                                    dra = fabs(ra-moonra);
                                    if (dra > M_PI) dra = 2*M_PI-dra;
                                    cosdis = sin(dec)*sin(moondec)+cos(dec)*cos(moondec)*cos(dra);
                                    if (cosdis>1) cosdis = 1.0;
                                    if (cosdis<-1) cosdis = -1.0;
                                    angular_sep_radians = acos(cosdis);
                                    angular_sep_degrees  =  angular_sep_radians * 180.0 / M_PI;

                                    scatter_function = 1.08572918*log(2.27E5 * (1.06 + cos(angular_sep_radians) * cos(angular_sep_radians))*
                                                                      pow((0.55/(central_wavelength)), 4.0) +
                                                                      exp((6.15 - angular_sep_degrees / 40.0)*2.30258509));
                                    if ( moonalt < 0 ) moon_magnitude = 10000;
                                    else moon_magnitude = lunar_illuminance - scatter_function;

                                    moon_magnitude = -2.5*log10(lunar_like_photon_count+pow(10.0, -0.4*moon_magnitude));
                                    darksky_magnitude = -2.5*log10(darksky_like_photon_count+pow(10.0, -0.4*zenith_v));


                                    // airglow variation
                                    double airglowv;
                                    long ax0, ax1, ay0, ay1;
                                    double dx, dy;
                                    find_linear_wrap(xv, platescale*15.0/3600, airglowScreenSize, &ax0, &ax1, &dx);
                                    find_linear_wrap(yv, platescale*15.0/3600, airglowScreenSize, &ay0, &ay1, &dy);

                                    airglowv = airglowvariation*((double)interpolate_bilinear_float_wrap(airglow, airglowScreenSize,
                                                                                                         ax0, ax1, dx, ay0, ay1, dy));

                                    if (diffusetype == 0) mag = domelight-2.5*log10(backRadius*backRadius);
                                    if (diffusetype == 1) mag = darksky_magnitude+airglowv-2.5*log10(backRadius*backRadius);
                                    if (diffusetype == 2) mag = moon_magnitude-2.5*log10(backRadius*backRadius);

                                    if (mag<100) {

                                        // fprintf(tempfile,"%lf %lf %d %lf\n",ra,dec,diffusetype,mag);

                                        nspatialpar = 0;


                                        sources.id.push_back(0.0);
                                        sources.ra.push_back(ra);
                                        sources.dec.push_back(dec);
                                        if (diffusetype == 0 && domewave == 0.0) sources.sedfilename[nsource] = dir+"/sed_dome.txt";
                                        if (diffusetype == 0 && domewave!=0.0) sources.sedfilename[nsource]="laser";
                                        if (diffusetype == 1) sources.sedfilename[nsource] = dir+"/darksky_sed.txt";
                                        if (diffusetype == 2) sources.sedfilename[nsource] = dir+"/lunar_sed.txt";
                                        sources.redshift.push_back(0.0);
                                        sources.gamma1.push_back(0.0);
                                        sources.gamma2.push_back(0.0);
                                        sources.kappa.push_back(0.0);
                                        sources.deltara.push_back(0.0);
                                        sources.deltadec.push_back(0.0);
                                        sources.type[nsource] = diffusetype;
                                        sources.norm[nsource] = pow(10.0, ((mag+48.6)/(-2.5)));
                                        sources.mag[nsource] = mag;

                                        double normwave = 500.0;
                                        oldsed = 0;
                                        /* read SED file */
                                        if (nsource > 0) {
                                            for (ii = 0;ii<nsource;ii++) {
                                                if (sources.sedfilename[ii] == sources.sedfilename[nsource]) {
                                                    sources.sedptr[nsource] = sources.sedptr[ii];
                                                    oldsed = 1;
                                                    goto skipsedread;
                                                }
                                            }
                                        }

                                        sources.sedptr[nsource] = nsedptr;

                                        if (sedptr == 0) {
                                            nreallocs = 0;
                                            sed_max = 4194304; // 4M elements, 32MB memory
                                            sed_w = (double*)malloc((long)(sed_max*sizeof(double)));
                                            sed_c = (double*)malloc((long)(sed_max*sizeof(double)));
                                        } else {
                                            if (sedptr > (sed_max-25000)) {
                                                ++nreallocs;
                                                sed_max = 2*sed_max;
                                                // fprintf(stdout,"%ld reallocs of sed_w and sed_c, now with sedptr=%ld and sed_max=%ld\n",
                                                //         nreallocs, sedptr, sed_max);
                                                sed_w = (double*)realloc(sed_w, (long)((sed_max)*sizeof(double)));
                                                sed_c = (double*)realloc(sed_c, (long)((sed_max)*sizeof(double)));
                                            }
                                        }

                                        lsedptr = 0;
                                        sprintf(tempstring, "%s", sources.sedfilename[nsource].c_str());
                                        if (flatdir == 1) {
                                            sptr = strtok_r(tempstring, "/", &sptr2);
                                            do {
                                                sptr3 = sptr;
                                                sptr = strtok_r(NULL, "/", &sptr2);
                                            } while (sptr!=NULL);
                                            sprintf(tempstring, "%s", sptr3);
                                        }

                                        if (strstr(tempstring, "laser")!=NULL) {
                                            normwave = domewave;
                                            closestw = 0; cdw = 1e30; badfile = 0;
                                            tempf1 = 0.0; tempf2 = 0.0;
                                            for (int k = 0; k < 5; k++) {

                                                if (k == 0) {tempf1=300.0; tempf2=0.0;}
                                                if (k == 1) {tempf1 = domewave-1.0; tempf2 = 0.0;}
                                                if (k == 2) {tempf1 = domewave; tempf2 = 1.98645e-16/(1e-7)/tempf1;}
                                                if (k == 3) {tempf1 = domewave+1.0; tempf2 = 0.0;}
                                                if (k == 4) {tempf1 = 1200.0; tempf2 = 0.0;}

                                                *(sed_w + sedptr + lsedptr) = tempf1;
                                                *(sed_c + sedptr + lsedptr) = tempf2*tempf1*1e-7/(1.98645e-16);

                                                dw = fabs(tempf1-domewave);

                                                if (dw < cdw) {
                                                    cdw = dw;
                                                    closestw = lsedptr;
                                                }
                                                lsedptr = lsedptr + 1;
                                            }

                                        } else {

                                            if (strstr(tempstring, ".gz") == NULL) {

                                                indafile = fopen(tempstring, "r");
                                                if (indafile == NULL) {fprintf(stderr, "Can't find SED file: %s\n", tempstring); exit(1);}

                                                closestw = 0; cdw = 1e30; badfile = 0;
                                                while (fgets(line, 4096, indafile)) {
                                                    sscanf(line, "%s %s", tempstring1, tempstring2);
                                                    tempf1 = strtod(tempstring1, NULL);
                                                    tempf2 = strtod(tempstring2, NULL);

                                                    *(sed_w + sedptr + lsedptr) = tempf1;
                                                    *(sed_c + sedptr + lsedptr) = tempf2*tempf1*1e-7/(1.98645e-16);
                                                    if (tempstring2[0] == 'n' || tempstring2[0] == 'N') {
                                                        if (badfile == 0) fprintf(stderr, "Error:   SED file: %s contains a NaN!\n", tempstring);
                                                        *(sed_c + sedptr + lsedptr) = 1e-20;
                                                        badfile = 1;
                                                    }

                                                    dw = fabs(tempf1 - 500.0);

                                                    if (dw < cdw) {
                                                        cdw = dw;
                                                        closestw = lsedptr;
                                                    }
                                                    lsedptr = lsedptr + 1;
                                                    if (lsedptr >= 25000) {
                                                        fprintf(stderr, "Error:  Too many lines in SED file: %s\n", tempstring);
                                                        exit(1);
                                                    }
                                                }
                                                fclose(indafile);

                                            } else {

                                                ingzfile = gzopen(tempstring, "r");
                                                if (ingzfile == NULL) {fprintf(stderr, "Can't find SED file: %s\n", tempstring); exit(1);}

                                                closestw = 0; cdw = 1e30; badfile = 0;
                                                while (gzgets(ingzfile, line, 4096)) {
                                                    sscanf(line, "%s %s", tempstring1, tempstring2);
                                                    tempf1 = strtod(tempstring1, NULL);
                                                    tempf2 = strtod(tempstring2, NULL);

                                                    *(sed_w + sedptr + lsedptr) = tempf1;
                                                    *(sed_c + sedptr + lsedptr) = tempf2*tempf1*1e-7/(1.98645e-16);
                                                    if (tempstring2[0] == 'n' || tempstring2[0] == 'N') {
                                                        if (badfile == 0) fprintf(stderr, "Error:   SED file: %s contains a NaN!\n", tempstring);
                                                        *(sed_c + sedptr + lsedptr) = 1e-20;
                                                        badfile = 1;
                                                    }

                                                    dw = fabs(tempf1 - 500.0);

                                                    if (dw < cdw) {
                                                        cdw = dw;
                                                        closestw = lsedptr;
                                                    }
                                                    lsedptr = lsedptr + 1;
                                                    if (lsedptr >= 25000) {
                                                        fprintf(stderr, "Error:  Too many lines in SED file: %s\n", tempstring);
                                                        exit(1);
                                                    }
                                                }
                                                gzclose(ingzfile);

                                            }

                                        }

                                        for (jj = 0; jj < lsedptr; jj++) {
                                            if (jj !=  0 && jj != (lsedptr-1)) {
                                                *(sed_c+sedptr+jj)=*(sed_c+sedptr+jj)*(*(sed_w+sedptr+jj+1)-*(sed_w+sedptr+jj-1))/2.0;
                                            }
                                            if (jj == 0) *(sed_c+sedptr+jj)=*(sed_c+sedptr+jj)*(*(sed_w+sedptr+jj+1)-*(sed_w+sedptr+jj));
                                            if (jj == (lsedptr-1)) *(sed_c+sedptr+jj)=*(sed_c+sedptr+jj)*(*(sed_w+sedptr+jj)-*(sed_w+sedptr+jj-1));
                                        }

                                        tempf1 = 0;
                                        for (jj = 0; jj < lsedptr; jj++) tempf1+=*(sed_c+sedptr+jj);
                                        for (jj = 0; jj < lsedptr; jj++) *(sed_c+sedptr+jj) = *(sed_c+sedptr+jj)/tempf1;

                                        if (closestw == 0)
                                            sed_dwdp[nsedptr] = (*(sed_w+sedptr+closestw+1)-*(sed_w+sedptr+closestw))/(*(sed_c+sedptr+closestw))/1.0;
                                        else
                                            sed_dwdp[nsedptr] = (*(sed_w+sedptr+closestw+1)-*(sed_w+sedptr+closestw-1))/(*(sed_c+sedptr+closestw))/2.0;
                                        if (*(sed_c+sedptr+closestw) <= 0.0) {
                                            printf("Error in SED file; 0 value at 500 nm\n");
                                            sed_dwdp[nsedptr] = 0.0;
                                        }

                                        for (jj = 1; jj < lsedptr; jj++) *(sed_c+sedptr+jj)+=*(sed_c+sedptr+jj-1);

                                        if (oldsed  ==  0) {
                                            sed_n[nsedptr] = lsedptr;
                                            sed_ptr[nsedptr] = sedptr;
                                            sedptr = sedptr+lsedptr;
                                            nsedptr++;
                                        }
                                        if (nsedptr >= 10000) { printf("Error:   Too many SED files\n"); exit(1); }

                                        /*moved up for catalog fix*/
                                    skipsedread:;
                                        sources.norm[nsource] = sources.norm[nsource]/(normwave)*(1+sources.redshift[nsource])*sed_dwdp[sources.sedptr[nsource]];



                                        sources.spatialname[nsource] = "gauss";
                                        if (sources.spatialname[nsource] == "gauss") {sources.spatialtype[nsource] = 2; nspatialpar = 1;}
                                        sources.spatialpar[nsource][0] = backRadius;
                                        sources.dustnamez[nsource] = "none";
                                        sources.dustname[nsource] = "none";

                                        setup_tangent(pra, pdec, &tpx, &tpy, &tpz);

                                        tangent(sources.ra[nsource]+sources.deltara[nsource], sources.dec[nsource]+sources.deltadec[nsource],
                                                &x, &y, &tpx, &tpy, &tpz);

                                        sources.vx[nsource] = x*cos(rotationangle)-y*sin(rotationangle);
                                        sources.vy[nsource] = x*sin(rotationangle)+y*cos(rotationangle);
                                        sources.vz[nsource] = -1.0;
                                        nn = sqrt((sources.vx[nsource])*(sources.vx[nsource])+
                                                  (sources.vy[nsource])*(sources.vy[nsource])+1);
                                        sources.vx[nsource] = sources.vx[nsource]/nn;
                                        sources.vy[nsource] = sources.vy[nsource]/nn;
                                        sources.vz[nsource] = sources.vz[nsource]/nn;
                                        nsource++;

                                    }

                                }


                            }
                        }
                    }

                }
            }
        }
    }

    return(0);
}


int Observation::parser () {


    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "Photon Raytrace\n");
    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "Installing Universe.\n");

    double altitude,  solaralt;
    int extraCommandFlag=0;
    int allocateZernikes=0;
    nphot = 0;
    totalnorm = 0.0;
    nsedptr = 0;
    sedptr = 0;
    nsource = 0;
    nimage = 0;
    nsurf = 0;
    maxr = 4180;
    minr = 2558;
    exptime = 15.0;
    vistime = 33.0;
    nsnap = 2;
    shuttererror = 0.0;
    timeoffset = 0.0;
    pra = 0.0;
    pdec = 0.0;
    rotationangle = 0.0;
    spiderangle = 0.0;
    zenith = 0.0;
    airmass = 1.0;
    azimuth = 0.0;
    windjitter = 2.5;
    rotationjitter = 1;
    elevationjitter = 0.02;
    azimuthjitter = 0.02;
    impurityvariation = 1;
    fieldanisotropy = 1;
    fringeflag = 1;
    deadlayer = 1;
    chargesharing = 1;
    pixelerror = 1;
    telescope_on = 1;
    coatingmode = 1;
    chargediffusion = 1;
    contaminationmode = 1;
    tracking_on = 1;
    detector_on = 1;
    diffraction_on = 1;
    pressure = 520;
    water_pressure = 8;
    temperature = 5;
    airrefraction = 1;
    raynorm = 1;
    o2norm = 1;
    o3norm = 1;
    h2onorm = 1;
    aerosoltau = 0.02;
    aerosolindex = -1.28;
    ranseed = -1;
    obsseed = 0;
    zernikemode = 1;
    atmospheric_dispersion = 1;
    atmosphericdispcenter = 1;
    outputdir = ".";
    seddir = "../data/SEDs";
    imagedir = "../data/images";
    datadir = "../data";
    instrdir = "../data/lsst";
    bindir = "../data/lsst";
    outputfilename = "focalplane";
    chipid = "R22_S11";
    trackingfile = ".";
    natmospherefile = 0;
    straylight = 1;
    straylightcut = 10.0;
    aperturemode = 0;
    ghostonly = 0;
    saturation = 1;
    eventfile = 0;
    opdfile = 0;
    eventFitsFileName = "output.fits";
    centroidfile = 0;
    throughputfile = 0;
    filter = 0;
    blooming = 1;
    obshistid = 0;
    pairid = 0;
    tai = 0.0;
    domeseeing = 0.1;
    finiteDistance = 0.0;
    transtol = 0.0;
    backAlpha = 0.1;
    backGamma = 2.0;
    backDelta = 80.0;
    screentol = 0.01;
    backBeta = 4.0;
    backRadius = 10.0;
    backBuffer = 100.0;
    np = 0.9;
    date = 1;
    siliconthickness = 100;
    overdepbias = -45.0;
    ccdtemp = 173;
    qevariation = 0.0;
    airglowvariation = 1.0;
    airglowScreenSize = 1024;
    miescatter_scat = 0.135;
    totalseeing = 0.67;
    flatdir = 0;
    atmdebug = 0;
    large_scale = 1.0;
    coarse_scale = 1.0;
    medium_scale = 1.0;
    fine_scale = 1.0;
    large_grid = 1;
    coarse_grid = 1;
    medium_grid = 1;
    fine_grid = 1;
    moonalt  =  - 1.0 * M_PI / 2.0 + 0.00001;
    moondist  =  M_PI - 0.00001;
    phaseang  =  M_PI - 0.00001;
    solarzen  =  M_PI - 0.00001;
    // zenith_v = 21.6437057;
    // + RngDouble() * 0.9;
    zenith_v = 22.08;
    moonra = 0;
    moondec = 0;
    domelight = 1000.0;
    domewave = 0.0;
    raydensity = 0.6;
    scalenumber = 8.0;
    checkpointtotal = 0;
    checkpointcount = 0;

    //update these values when instrdir is set
    //variables in focalplanelayout.txt
    centerx = -1;
    centery = -1;
    pixsize = -1;
    pixelsx = -1;
    pixelsy = -1;
    minx = -1;
    miny = -1;
    maxx = -1;
    maxy = -1;
    //variables in location.txt
    groundlevel = -1;
    xtelloc = -1;
    ytelloc = -1;
    latitude = -1;
    //variables in central_wavelengths.txt
    central_wavelength = -1;
    platescale = -1;
    //variables in silicon.txt
    well_depth = -1;
    nbulk = -1;
    nf = -1;
    nb = -1;
    sf = -1;
    sb = -1;
    //variables in tracking.txt
    windjitter = -1;
    rotationjitter = -1;
    elevationjitter = -1;
    azimuthjitter = -1;

    /* atmosphere parameter arrays */
    atmospherefile.resize(MAX_LAYER);
    cloudfile.resize(MAX_LAYER);
    seefactor.resize(MAX_LAYER, 0);
    wind.resize(MAX_LAYER, 0);
    winddir.resize(MAX_LAYER, 0);
    outerscale.resize(MAX_LAYER, 0);
    height.resize(MAX_LAYER, 0);
    densityfluctuation.resize(MAX_LAYER, 0);
    densitymean.resize(MAX_LAYER, 1.0);
    cloudmean.resize(MAX_LAYER, 0);
    cloudvary.resize(MAX_LAYER, 0);
    dtau = (double*)calloc(MAX_LAYER,sizeof(double));

    /* telescope parameter arrays */
    izernike.resize(MAX_SURF);
    body.resize(MAX_SURF);
    for (int i = 0;i<MAX_SURF;i++) body[i].resize(6, 0);
    ghost.resize(MAX_SURF, 0);
    feaflag.resize(MAX_SURF, 0);
    feafile.resize(MAX_SURF);


    /* sky parameter arrays */
    sources.vx = (double*)calloc(MAX_SOURCE, sizeof(double));
    sources.vy = (double*)calloc(MAX_SOURCE, sizeof(double));
    sources.vz = (double*)calloc(MAX_SOURCE, sizeof(double));
    sources.norm = (double*)calloc(MAX_SOURCE, sizeof(double));
    sources.mag = (double*)calloc(MAX_SOURCE, sizeof(double));
    sources.type = (int*)calloc(MAX_SOURCE, sizeof(int));
    sed_corr = (double*)calloc(MAX_SOURCE, sizeof(double));
    sed_dwdp = (double*)calloc(MAX_SOURCE, sizeof(double));
    sources.spatialtype = (int*)calloc(MAX_SOURCE, sizeof(int));
    sources.dusttype = (int*)calloc(MAX_SOURCE, sizeof(int));
    sources.dusttypez = (int*)calloc(MAX_SOURCE, sizeof(int));
    sources.sedptr = (long*)calloc(MAX_SOURCE, sizeof(long));
    // sources.id = (double*)calloc(MAX_SOURCE, sizeof(double));
    sources.skysameas = (long*)calloc(MAX_SOURCE, sizeof(long));
    sed_n = (long*)calloc(MAX_SOURCE, sizeof(long));
    sed_ptr = (long*)calloc(MAX_SOURCE, sizeof(long));
    source_xpos = (long long*)calloc(MAX_SOURCE, sizeof(long long));
    source_ypos = (long long*)calloc(MAX_SOURCE, sizeof(long long));
    source_photon = (long long*)calloc(MAX_SOURCE, sizeof(long long));
    sources.spatialpar = (double**)calloc(MAX_SOURCE, sizeof(double*));
    for (int i = 0;i<MAX_SOURCE;i++) sources.spatialpar[i] = (double*)calloc(6, sizeof(double));
    sources.dustpar = (double**)calloc(MAX_SOURCE, sizeof(double*));
    for (int i = 0;i<MAX_SOURCE;i++) sources.dustpar[i] = (double*)calloc(2, sizeof(double));
    sources.dustparz = (double**)calloc(MAX_SOURCE, sizeof(double*));
    for (int i = 0;i<MAX_SOURCE;i++) sources.dustparz[i] = (double*)calloc(2, sizeof(double));
    sources.sedfilename.resize(MAX_SOURCE);
    sources.spatialname.resize(MAX_SOURCE);
    sources.dustname.resize(MAX_SOURCE);
    sources.dustnamez.resize(MAX_SOURCE);


    readText pars(std::cin);
    for (size_t t(0); t<pars.getSize(); t++) {
        std::string line(pars[t]);
        std::istringstream iss(line);
        std::string keyName;
        iss >> keyName;

        if (keyName == "object") {
            std::string object;
            std::getline(iss, object);
            addSource(object,  3);
            continue;
        }

        readText::get(line, "outputdir", outputdir);
        readText::get(line, "seddir", seddir);
        readText::get(line, "imagedir", imagedir);
        readText::get(line, "datadir", datadir);
        readText::get(line, "instrdir", instrdir);
        readText::get(line, "bindir", bindir);
        readText::get(line, "telescopemode", telescope_on);
        readText::get(line, "impurityvariation", impurityvariation);
        readText::get(line, "fieldanisotropy", fieldanisotropy);
        readText::get(line, "fringing", fringeflag);
        readText::get(line, "deadlayer", deadlayer);
        readText::get(line, "chargesharing", chargesharing);
        readText::get(line, "pixelerror", pixelerror);
        readText::get(line, "chargediffusion", chargediffusion);
        readText::get(line, "coatingmode", coatingmode);
        readText::get(line, "contaminationmode", contaminationmode);
        readText::get(line, "trackingmode", tracking_on);
        readText::get(line, "detectormode", detector_on);
        readText::get(line, "diffractionmode", diffraction_on);
        readText::get(line, "zernikemode", zernikemode);
        readText::get(line, "straylight", straylight);
        readText::get(line, "straylightcut", straylightcut);// not exposed
        readText::get(line, "ghost", ghost);// no exposed
        readText::get(line, "ghostonly", ghostonly);
        readText::get(line, "aperturemode", aperturemode);
        readText::get(line, "minr", minr);
        readText::get(line, "maxr", maxr);
        readText::get(line, "exptime", exptime);
        if (readText::getKey(line, "nsnap", nsnap)) extraCommandFlag=1;
        readText::get(line, "shuttererror", shuttererror);
        readText::get(line, "timeoffset", timeoffset);
        readText::get(line, "finitedistance", finiteDistance);
        readText::get(line, "transtol", transtol); // not exposed
        readText::get(line, "np", np); // not exposed
        readText::get(line, "backalpha", backAlpha);// not exposed
        readText::get(line, "backgamma", backGamma);// not exposed
        readText::get(line, "backdelta", backDelta);// not exposed
        readText::get(line, "screentol", screentol);// not exposed
        readText::get(line, "backbeta", backBeta);//  not exposed
        readText::get(line, "backradius", backRadius);// not exposed
        readText::get(line, "backbuffer", backBuffer);// not exposed
        readText::get(line, "date", date);// not exposed
        readText::get(line, "flatdir", flatdir);// not exposed
        readText::get(line, "atmdebug", atmdebug);// not exposed
        readText::get(line, "large_grid", large_grid);// not exposed
        readText::get(line, "coarse_grid", coarse_grid);// not exposed
        readText::get(line, "medium_grid", medium_grid);// not exposed
        readText::get(line, "fine_grid", fine_grid);// not exposed
        readText::get(line, "large_scale", large_scale);// not exposed
        readText::get(line, "coarse_scale", coarse_scale);// not exposed
        readText::get(line, "medium_scale", medium_scale);// not exposed
        readText::get(line, "fine_scale", fine_scale);// not exposed
        readText::get(line, "opdfile", opdfile);// not exposed
        readText::get(line, "filter", filter);
        readText::get(line, "saturation", saturation);
        readText::get(line, "blooming", blooming);
        readText::get(line, "eventfile", eventfile);
        readText::get(line, "eventFitsFileName", eventFitsFileName);
        readText::get(line, "centroidfile", centroidfile);
        readText::get(line, "throughputfile", throughputfile);
        readText::get(line, "well_depth", well_depth);
        readText::get(line, "nbulk", nbulk);
        readText::get(line, "nf", nf);
        readText::get(line, "nb", nb);
        readText::get(line, "sf", sf);
        readText::get(line, "sb", sb);
        readText::get(line, "siliconthickness", siliconthickness);
        readText::get(line, "overdepbias", overdepbias);
        readText::get(line, "ccdtemp", ccdtemp);
        readText::get(line, "qevariation", qevariation);
        readText::get(line, "obshistid", obshistid); //
        readText::get(line, "exposureid", pairid);//
        readText::get(line, "tai", tai);//
        readText::get(line, "windjitter", windjitter);
        readText::get(line, "rotationjitter", rotationjitter);
        readText::get(line, "elevationjitter", elevationjitter);
        readText::get(line, "azimuthjitter", azimuthjitter);
        if ( allocateZernikes == 0 && (keyName == "izernike" || keyName == "ichebyshev")) {
            if (keyName == "izernike" ) {
                NZERN=22;
                pertType.assign("zern");
            } else if (keyName == "ichebyshev" ) {
                NZERN=22;
                pertType.assign("chebyshev");
            }
            for (int i = 0;i<MAX_SURF;i++) izernike[i].resize(NZERN, 0);
            allocateZernikes=1;
        }
        readText::get(line, "izernike", izernike);
        readText::get(line, "ichebyshev", izernike);
        readText::get(line, "body", body);
        readText::get(line, "natmospherefile", natmospherefile);
        readText::get(line, "atmospherefile", atmospherefile);
        readText::get(line, "cloudfile", cloudfile);
        readText::get(line, "trackingfile", trackingfile);
        readText::get(line, "chipid", chipid);
        readText::get(line, "seeing", seefactor);
        readText::get(line, "wind", wind);
        readText::get(line, "winddir", winddir);
        readText::get(line, "outerscale", outerscale);
        readText::get(line, "height", height);
        readText::get(line, "densityfluc", densityfluctuation);
        readText::get(line, "densitymean", densitymean);
        readText::get(line, "cloudmean", cloudmean);
        readText::get(line, "cloudvary", cloudvary);
        readText::get(line, "atmosphericdispersion", atmospheric_dispersion);
        readText::get(line, "atmosphericdispcenter", atmosphericdispcenter);
        readText::get(line, "seed", ranseed);
        readText::get(line, "obsseed", obsseed);//
        readText::get(line, "vistime", vistime);//
        readText::get(line, "pressure", pressure);
        readText::get(line, "waterpressure", water_pressure);
        readText::get(line, "temperature", temperature);
        readText::get(line, "airrefraction", airrefraction);
        readText::get(line, "reldensity", raynorm);
        readText::get(line, "relo2", o2norm);
        readText::get(line, "relo3", o3norm);
        readText::get(line, "relh2o", h2onorm);
        readText::get(line, "aerosoltau", aerosoltau);
        readText::get(line, "aerosolindex", aerosolindex);
        readText::get(line, "lascatprob", miescatter_scat);
        readText::get(line, "domeseeing", domeseeing);
        readText::get(line, "airglowvariation", airglowvariation);//
        readText::get(line, "totalseeing", totalseeing);//
        readText::get(line, "zenith_v", zenith_v);//
        readText::get(line, "domelight", domelight);//
        readText::get(line, "telconfig", telconfig);//
        readText::get(line, "checkpointcount", checkpointcount);//
        readText::get(line, "checkpointtotal", checkpointtotal);//
        readText::get(line, "domewave", domewave);//
        readText::get(line, "raydensity", raydensity);//
        readText::get(line, "scalenumber", scalenumber);//
        readText::get(line, "centerx", centerx);
        readText::get(line, "centery", centery);
        readText::get(line, "pixelsize", pixsize);
        readText::get(line, "pixelsx", pixelsx);
        readText::get(line, "pixelsy", pixelsy);
        readText::get(line, "minx", minx);
        readText::get(line, "miny", miny);
        readText::get(line, "maxx", maxx);
        readText::get(line, "maxy", maxy);
        readText::get(line, "wavelength", central_wavelength);//
        readText::get(line, "platescale", platescale);
        readText::get(line, "groundlevel", groundlevel);
        readText::get(line, "xtellocation", xtelloc);
        readText::get(line, "ytellocation", ytelloc);// up to here
        if (readText::getKey(line, "latitude", latitude)) latitude *= M_PI/180.0;
        if (readText::getKey(line, "longitude", longitude)) longitude *= M_PI/180.0;
        if (readText::getKey(line, "pointingra", pra)) pra *= M_PI/180.0;
        if (readText::getKey(line, "pointingdec", pdec)) pdec *= M_PI/180.0;
        if (readText::getKey(line, "rotationangle", rotationangle)) rotationangle *= (-M_PI/180.);
        if (readText::getKey(line, "spiderangle", spiderangle)) spiderangle *= M_PI/180.;
        if (readText::getKey(line, "altitude", altitude)) {
            altitude *= M_PI/180.;
            zenith = M_PI/2.0 - altitude;
        }
        if (readText::getKey(line, "zenith", zenith)) zenith *= M_PI/180.0;
        if (readText::getKey(line, "azimuth", azimuth)) azimuth *= M_PI/180.0;
        if (readText::getKey(line, "moonra", moonra)) moonra *= M_PI/180.0;
        if (readText::getKey(line, "moondec", moondec)) moondec *= M_PI/180.0;
        if (readText::getKey(line, "moonalt", moonalt)) moonalt *= M_PI/180.0;
        if (readText::getKey(line, "moondist", moondist)) moondist *= M_PI/180.0;
        if (readText::getKey(line, "phaseang", phaseang)) phaseang = M_PI-phaseang*M_PI/100.0;
        if (readText::getKey(line, "solaralt", solaralt)) {
            solaralt *= M_PI/180.0;
            solarzen = M_PI/2.0 - solaralt;
        }
        if (readText::getKey(line, "solarzen", solarzen)) solarzen *= M_PI/180.;
        if (keyName == "clearperturbations") {
            for (int i = 0;i<MAX_SURF;i++) for (int j = 0; j < 6; j++) body[i][j] = 0.0;
            for (int i = 0;i<MAX_SURF;i++) for (int j = 0;j < NZERN; j++) izernike[i][j] = 0.0;
        }
        if (keyName == "cleartracking") {
            rotationjitter = 0.0;
            elevationjitter = 0.0;
            azimuthjitter = 0.0;
        }
        if (keyName == "clearclouds") {
            for (int i = 0;i<MAX_LAYER;i++) {
                cloudmean[i] = 0.0;
                cloudvary[i] = 0.0;
            }
        }
        if (keyName == "clearopacity") {
            h2onorm = 0.0;
            raynorm = 0.0;
            o2norm = 0.0;
            o3norm = 0.0;
            aerosoltau = 0.0;
            for (int i = 0;i<MAX_LAYER;i++) {
                cloudmean[i] = 0.0;
                cloudvary[i] = 0.0;
                densitymean[i] = 0.0;
            }
        }
        if (keyName == "clearturbulence") {
            for (int i = 0;i<MAX_LAYER;i++) {
                seefactor[i] = 0.0;
                densityfluctuation[i] = 0.0;
            }
            domeseeing = 0.0;
        }
        if (keyName == "cleardefects") {
            impurityvariation = 0;
            fieldanisotropy = 0;
            deadlayer = 0;
            chargesharing = 0;
            pixelerror = 0;
        }
        if (keyName == "fea") {
            long surfaceIndex;
            iss >> surfaceIndex;
            iss >> feafile[surfaceIndex];
            feaflag[surfaceIndex] = 1;
        }

        if (extraCommandFlag > 1) {
            extraCommandString.push_back(line);
            extraCommandFlag++;
        }
        if (extraCommandFlag == 1) extraCommandFlag++;

    }

    if (telconfig!=2 && telconfig!=3) domelight = 1000.0;

    std::ostringstream outfile;
    unsigned pos = instrdir.rfind("/")+1;
    for (unsigned i = pos; i<instrdir.length(); i++) outfile<<instrdir[i];
    outfile << "_e_"  << obshistid << "_" << chipid << "_E" << std::setfill('0') << std::setw(3) << pairid;
    outputfilename = outfile.str();

    if (flatdir == 1) {
        instrdir  =  ".";
        bindir  =  ".";
        imagedir = ".";
        std::ostringstream tarName;
        tarName << "SEDs_" << obshistid << "_" << chipid << ".tar";
        std::ifstream tarFile(tarName.str().c_str());
        if (tarFile.good()) {
            std::string tarCommand = "tar xf " + tarName.str();
            system(tarCommand.c_str());
        }
    }

    focalplanefile = instrdir + "/focalplanelayout.txt";
    std::istringstream focalplanePars(readText::get(focalplanefile, chipid));
    double centerx_t, centery_t, pixsize_t;
    long pixelsx_t, pixelsy_t;
    double angle1, angle2;
    std::string grouptype;
    focalplanePars>>centerx_t>>centery_t>>pixsize_t>>pixelsx_t>>pixelsy_t>>devtype>>devvalue>>grouptype>>chipangle>>angle1>>angle2>>decenterx>>decentery;
    decenterx *= 1000.0;
    decentery *= 1000.0;
    chipangle *= M_PI/180.0;
    if (centerx == -1) centerx=centerx_t;
    if (centery == -1) centery = centery_t;
    if (pixsize == -1) pixsize = pixsize_t;
    if (pixelsx == -1) pixelsx = pixelsx_t;
    if (pixelsy == -1) pixelsy = pixelsy_t;
    if (minx == -1) minx = 0;
    if (miny == -1) miny = 0;
    if (maxx == -1) maxx = pixelsx-1;
    if (maxy == -1) maxy = pixelsy-1;

    std::istringstream wavelengthPars(readText::get(instrdir + "/central_wavelengths.txt", filter));
    double wavelength_t, platescale_t;
    wavelengthPars >> wavelength_t >> platescale_t;
    if (central_wavelength == -1) central_wavelength = wavelength_t;
    if (platescale == -1) platescale = platescale_t;

    readText locationPars(instrdir + "/location.txt");
    for (size_t t(0); t < locationPars.getSize(); t++) {
        std::string line(locationPars[t]);
        if (groundlevel  ==  -1) readText::get(line, "groundlevel", groundlevel);
        if (xtelloc == -1) readText::get(line, "xtellocation", xtelloc);
        if (ytelloc == -1) readText::get(line, "ytellocation", ytelloc);
        if (latitude == -1) if (readText::getKey(line, "latitude", latitude)) latitude *= M_PI/180.;
    }

    readText siliconPars(instrdir + "/silicon.txt");
    for (size_t t(0); t < siliconPars.getSize(); t++) {
        std::string line(siliconPars[t]);
        if (well_depth == -1) readText::get(line, "wellDepth", well_depth);
        if (nbulk == -1) readText::get(line, "nbulk", nbulk);
        if (nf == -1) readText::get(line, "nf", nf);
        if (nb == -1) readText::get(line, "nb", nb);
        if (sf == -1) readText::get(line, "sf", sf);
        if (sb == -1) readText::get(line, "sb", sb);
    }

    readText trackingPars(instrdir + "/tracking.txt");
    for (size_t t(0); t < trackingPars.getSize(); t++) {
        std::string line(trackingPars[t]);
        if (windjitter == -1) readText::get(line, "windjitter", windjitter);
        if (rotationjitter == -1) readText::get(line, "rotationjitter", rotationjitter);
        if (elevationjitter == -1) readText::get(line, "elevationjitter", elevationjitter);
        if (azimuthjitter == -1) readText::get(line, "azimuthjitter", azimuthjitter);
    }

    return(0);

}


int Observation::settings() {

    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    fprintf(stdout, "Basic Setup \n");
    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "[outputdir] Output directory: %s\n", outputdir.c_str());
    fprintf(stdout, "[outputfilename] Output filename: %s\n", outputfilename.c_str());
    fprintf(stdout, "[seddir] SED directory: %s\n", seddir.c_str());
    fprintf(stdout, "[imagedir] Image directory: %s\n", imagedir.c_str());
    fprintf(stdout, "[centroidfile] Output centroid file (0=no/1=yes):          %ld\n", centroidfile);
    fprintf(stdout, "[throughputfile] Output throughput file (0=no/1=yes):      %ld\n", throughputfile);
    fprintf(stdout, "[eventfile] Output event file (0=no/1=yes):                %ld\n", eventfile);
    fprintf(stdout, "[eventFitsFileName] Output event Fits file name:           %s\n", eventFitsFileName.c_str());
    fprintf(stdout, "[bindir] Binary Directory:                                 %s\n", bindir.c_str());
    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "Module Switches\n");
    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "[telescopemode] Telescope mode (0=off/1=on):               %ld\n", telescope_on);
    fprintf(stdout, "[trackingmode] Tracking mode (0=off/1=on):                 %ld\n", tracking_on);
    fprintf(stdout, "[detectormode] Detector mode (0=off/1=on):                 %ld\n", detector_on);
    fprintf(stdout, "[diffractionmode] Diffraction mode (0=off/1=on):           %ld\n", diffraction_on);
    fprintf(stdout, "[zernikemode] Zernike mode (0=off/1=on):                   %ld\n", zernikemode);
    fprintf(stdout, "[straylight] Straylight mode (0=off/1=on):                 %ld\n", straylight);
    fprintf(stdout, "[aperturemode] Aperture mode (0=normal/1=on):              %ld\n", aperturemode);
    fprintf(stdout, "[ghostonly] Ghost-only mode (0=normal/1=on):               %ld\n", ghostonly);
    fprintf(stdout, "[saturation] Saturation mode (0=off/1=on):                 %ld\n", saturation);
    fprintf(stdout, "[blooming] Blooming mode (0=off/1=on):                     %ld\n", blooming);
    fprintf(stdout, "[atmosphericdispersion] Atmos. Dispersion (0=off/1=on):    %ld\n", atmospheric_dispersion);
    fprintf(stdout, "[atmosphericdispcenter] Atmos. Disp. Ctr. Corr.:           %ld\n", atmosphericdispcenter);
    fprintf(stdout, "[impurityvariation] Impurity Variation (0=off/1=on):       %d\n", impurityvariation);
    fprintf(stdout, "[fieldanisotropy] Field Anisotropy (0=off/1=on):           %d\n", fieldanisotropy);
    fprintf(stdout, "[fringing] Fringing (0=off/1=on):                          %d\n", fringeflag);
    fprintf(stdout, "[deadlayer] Dead Layer (0=off/1=on):                       %d\n", deadlayer);
    fprintf(stdout, "[chargediffusion] Charge Diffusion (0=off/1=on):           %d\n", chargediffusion);
    fprintf(stdout, "[chargesharing] Charge Sharing (0=off/1=on):               %d\n", chargesharing);
    fprintf(stdout, "[pixelerror] Pixel Error (0=off/1=on):                     %d\n", pixelerror);
    fprintf(stdout, "[coatingmode] Coating Mode (0=off/1=on):                   %ld\n", coatingmode);
    fprintf(stdout, "[contaminationmode] Contamination Mode (0=off/1=on):       %ld\n", contaminationmode);
    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "Telescope Operator and Bookkeeping\n");
    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "Pointing RA (degrees):                                     %lf\n", pra/DEGREE);
    fprintf(stdout, "Pointing Dec (degrees):                                    %lf\n", pdec/DEGREE);
    fprintf(stdout, "Rotation Angle (rotSkyPos) (degrees):                      %lf\n", -rotationangle/DEGREE);
    fprintf(stdout, "Angle of Spider (rotTelPos) (degrees):                     %lf\n", spiderangle/DEGREE);
    fprintf(stdout, "Zenith Angle (degrees):                                    %lf\n", zenith/DEGREE);
    fprintf(stdout, "Azimuthal Angle (degrees):                                 %lf\n", azimuth/DEGREE);
    fprintf(stdout, "Filter (number starting with 0):                           %ld\n", filter);
    fprintf(stdout, "Random seed:                                               %ld\n", ranseed);
    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "Instantaneous Instrument and Site Characteristics\n");
    fprintf(stdout, "------------------------------------------------------------------------------------------\n");
    fprintf(stdout, "[instrdir] Instrument & Site Directory:                    %s\n", instrdir.c_str());

    //optics
    fprintf(stdout, "[platescale] Plate Scale:                                  %lf\n", platescale);
    fprintf(stdout, "[minr] Minimum aperture radius:                            %lf\n", minr);
    fprintf(stdout, "[maxr] Maximum aperture radius:                            %lf\n", maxr);
    fprintf(stdout, "[chipid] Chip/Amplifier ID:                                %s\n", chipid.c_str());
    fprintf(stdout, "[centerx] Chip center x (microns):                         %lf\n", centerx);
    fprintf(stdout, "[centery] Chip center y (microns):                         %lf\n", centery);
    fprintf(stdout, "[pixelsx] Chip x pixels:                                   %ld\n", pixelsx);
    fprintf(stdout, "[pixelsy] Chip y pixels:                                   %ld\n", pixelsy);
    fprintf(stdout, "[minx] Minimum x pixel of amplifier:                       %ld\n", minx);
    fprintf(stdout, "[maxx] Maximum x pixel of amplifier:                       %ld\n", maxx);
    fprintf(stdout, "[miny] Minimum y pixel of amplifier:                       %ld\n", miny);
    fprintf(stdout, "[maxy] Maximum y pixel of amplifier:                       %ld\n", maxy);
    fprintf(stdout, "[pixelsize] Pixel Size (microns):                          %lf\n", pixsize);
    fprintf(stdout, "[welldepth] Full well depth:                               %ld\n", well_depth);
    fprintf(stdout, "[nbulk] Bulk doping density:                               %e\n", nbulk);
    fprintf(stdout, "[nf] Front side doping density:                            %e\n", nf);
    fprintf(stdout, "[nb] Back side doping density:                             %e\n", nb);
    fprintf(stdout, "[sf] Front side doping scale:                              %e\n", sf);
    fprintf(stdout, "[sb] Back side doping scale:                               %e\n", sb);
    fprintf(stdout, "[siliconthickness] Silicon Thickness (microns):            %f\n", siliconthickness);
    fprintf(stdout, "[overdepbias] Over depletion bias (volts):                 %f\n", overdepbias);
    fprintf(stdout, "[ccdtemp] CCD temperature (K):                             %f\n", ccdtemp);
    fprintf(stdout, "[qevariation] QE variation:                                %f\n", qevariation);
    fprintf(stdout, "[exptime] Exposure time (s):                               %lf\n", exptime);
    fprintf(stdout, "[nsnap] Number of Snaps:                                   %ld\n", nsnap);
    fprintf(stdout, "[shuttererror] Shutter error (s):                          %lf\n", shuttererror);
    fprintf(stdout, "[timeoffset] Time offset (s):                              %lf\n", timeoffset);
    fprintf(stdout, "[windjitter] Wind Jitter (degrees):                        %lf\n", windjitter);
    fprintf(stdout, "[rotationjitter] Rotation Jitter (arcseconds):             %lf\n", rotationjitter);
    fprintf(stdout, "[elevationjitter] Elevation Jitter (arcseconds):           %lf\n", elevationjitter);
    fprintf(stdout, "[azimuthjitter] Azimuthal Jitter (arcseconds):             %lf\n", azimuthjitter);
    fprintf(stdout, "[izernike optic# zernike#] Zernike amplitude:        \n");
    for (long i = 0; i < npertsurf; i++) {
        for (long j = 0; j < NZERN/2; j++) {
            printf("%lf ", izernike[i][j]);
        } printf("\n");
        printf(" ");
        for (long j = NZERN/2; j < NZERN; j++) {
            printf("%lf ", izernike[i][j]);
        } printf("\n");
    }
    fprintf(stdout, "[body optic# dof#] Body motion of optics:                    \n");
    for (long i = 0; i < nsurf; i++) {
        for (long j = 0; j < 6; j++) {
            printf("%lf ", body[i][j]);
        } printf("\n");
    }
    fprintf(stdout, "[lascatprob] Large angle scattering probability:           %lf\n", miescatter_scat);

    //site
    fprintf(stdout, "[domeseeing] Dome seeing:                                  %lf\n", domeseeing);
    fprintf(stdout, "[groundlevel] Ground level (m):                            %lf\n", groundlevel);
    fprintf(stdout, "[xtellocation] X Telescope location (m):                   %lf\n", xtelloc);
    fprintf(stdout, "[ytellocation] Y Telescope location (m):                   %lf\n", ytelloc);
    fprintf(stdout, "[latitude] Telescope latitude (degrees):                   %lf\n", latitude/DEGREE);
    fprintf(stdout, "[longitude] Telescope longitude (degrees):                 %lf\n", longitude/DEGREE);
    fprintf(stdout, "[pressure] Air pressure (mmHg):                            %lf\n", pressure);
    fprintf(stdout, "[waterpressure] Water vapor pressure (mmHg):               %lf\n", water_pressure);
    fprintf(stdout, "[temperature] Ground Temperature (degrees C):              %lf\n", temperature);
    fprintf(stdout, "[reldensity] Relative density:                             %lf\n", raynorm);
    fprintf(stdout, "[relo2] Relative O2 fraction:                              %lf\n", o2norm);
    fprintf(stdout, "[relh2o] Relative H2O fraction:                            %lf\n", h2onorm);
    fprintf(stdout, "[aerosoltau] Aerosol optical depth:                        %lf\n", aerosoltau);
    fprintf(stdout, "[aerosolindex] Aerosol index:                              %lf\n", aerosolindex);
    fprintf(stdout, "[relo3] Relative O3 fraction:                              %lf\n", o3norm);
    fprintf(stdout, "[natmospherefile] Number of atmosphere layers:             %ld\n", natmospherefile);
    fprintf(stdout, "[seeing layer#] Seeing at 5000 Angstroms (arcsec):\n");
    for (long i = 0; i < natmospherefile; i++) {
        printf("%lf ", seefactor[i]/pow(1/cos(zenith), 0.6));
    }
    printf("\n");
    fprintf(stdout, "[wind layer#] Wind speed (m/s):\n");
    for (long i = 0; i < natmospherefile; i++) {
        printf("%lf ", wind[i]);
    }
    printf("\n");
    fprintf(stdout, "[winddir layer#] Wind direction (degrees):\n");
    for (long i = 0; i < natmospherefile; i++) {
        printf("%lf ", winddir[i]);
    }
    printf("\n");
    fprintf(stdout, "[height layer#] Layer height (km):\n");
    for (long i = 0; i < natmospherefile; i++) {
        printf("%lf ", height[i]);
    }
    printf("\n");
    fprintf(stdout, "[outerscale layer#] Outer scale (m):\n");
    for (long i = 0; i < natmospherefile; i++) {
        printf("%lf ", outerscale[i]);
    }
    printf("\n");
    fprintf(stdout, "[densityfluc layer#] Density fluctuation:\n");
    for (long i = 0; i < natmospherefile; i++) {
        printf("%lf ", densityfluctuation[i]);
    }
    printf("\n");
    fprintf(stdout, "[densitymean layer#] Density mean:\n");
    for (long i = 0; i < natmospherefile; i++) {
        printf("%lf ", densitymean[i]);
    }
    printf("\n");
    fprintf(stdout, "[cloudmean layer#] Mean cloud extinction (km):\n");
    for (long i = 0; i < natmospherefile; i++) {
        printf("%lf ", cloudmean[i]);
    }
    printf("\n");
    fprintf(stdout, "[cloudvary layer#] Variation of cloud extinction (km):\n");
    for (long i = 0; i < natmospherefile; i++) {
        printf("%lf ", cloudvary[i]);
    }
    printf("\n");

    fprintf(stdout, "------------------------------------------------------------------------------------------\n");

    return(1);

}


int Observation::header( fitsfile *faptr) {

    long i; int status;
    char tempstring[4096];
    char tempstring2[4096];
    double tempf1;
    double starttime;
    status=0;

    for (long i = 0; i < (long)extraCommandString.size(); i++) {
        sprintf(tempstring, "PHOV%ld",i);
        sprintf(tempstring2, "Physics Override Command %ld",i);
        fits_write_key(faptr, TSTRING, (char*)tempstring, (char*)extraCommandString[i].c_str(), (char*)tempstring2, &status);
    }
    fits_write_key(faptr, TLONG, (char*)"OBSID", &obshistid, (char*)"Opsim observation ID", &status);
    sprintf(tempstring, "%ld", obshistid);
    fits_write_key(faptr, TSTRING, (char*)"DATASET", tempstring, (char*)"Opsim observation ID", &status);
    starttime=tai + (timeoffset - exptime/2.0)/24./3600.;
    fits_write_key(faptr, TDOUBLE, (char*)"TAI", &starttime, (char*)"International Atomic Time scale", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"MJD-OBS", &starttime, (char*)"Modified Julian date (also TAI)", &status);

    fits_write_key(faptr, TSTRING, (char*)"OUTPDIR", (char*)outputdir.c_str(), (char*)"Output directory", &status);
    fits_write_key(faptr, TSTRING, (char*)"OUTFILE", (char*)outputfilename.c_str(), (char*)"Output filename", &status);
    fits_write_key(faptr, TSTRING, (char*)"SEDDIR", (char*)seddir.c_str(), (char*)"SED directory", &status);
    fits_write_key(faptr, TSTRING, (char*)"IMGDIR", (char*)imagedir.c_str(), (char*)"Image directory", &status);
    fits_write_key(faptr, TLONG, (char*)"TELMODE", &telescope_on, (char*)"Telescope mode (0=off/1=on)", &status);
    fits_write_key(faptr, TLONG, (char*)"TRKMODE", &tracking_on, (char*)"Tracking mode (0=off/1=on)", &status);
    fits_write_key(faptr, TLONG, (char*)"DIFMODE", &diffraction_on, (char*)"Diffraction mode (0=off/1=on)", &status);
    fits_write_key(faptr, TLONG, (char*)"DETMODE", &detector_on, (char*)"Detector mode (0=off/1=on)", &status);
    fits_write_key(faptr, TLONG, (char*)"ZERMODE", &zernikemode, (char*)"Zernike mode (0=off/1=on)", &status);
    fits_write_key(faptr, TLONG, (char*)"STRLGHT", &straylight, (char*)"Straylight mode (0=off/1=on)", &status);
    fits_write_key(faptr, TLONG, (char*)"APRMODE", &aperturemode, (char*)"Aperture mode (0=normal/1=on)", &status);
    fits_write_key(faptr, TLONG, (char*)"GHOMODE", &ghostonly, (char*)"Ghost-only mode (0=normal/1=on)", &status);
    fits_write_key(faptr, TLONG, (char*)"SATMODE", &saturation, (char*)"Saturation mode (0=off/1=on)", &status);
    fits_write_key(faptr, TLONG, (char*)"BLOOMNG", &blooming, (char*)"Blooming mode (0=off/1=on)", &status);
    fits_write_key(faptr, TLONG, (char*)"EVTFILE", &eventfile, (char*)"Output event file (0=no/1=yes)", &status);

    fits_write_key(faptr, TSTRING, (char*)"EVTFITSFILE", (char*)eventFitsFileName.c_str(), (char*)"Event Fits file name", &status);

    fits_write_key(faptr, TLONG, (char*)"THRFILE", &throughputfile, (char*)"Output throughput file (0=no/1=yes)", &status);
    tempf1 = (3600.0*1000.0)/platescale;
    fits_write_key(faptr, TDOUBLE, (char*)"PLTSCAL", &tempf1, (char*)"Approx. Plate scale (arcsec/mm)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"MINR", &minr, (char*)"Minimum aperture radius", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"MAXR", &maxr, (char*)"Maximum aperture radius", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"EXPTIME", &exptime, (char*)"Exposure time", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DARKTIME", &exptime, (char*)"Actual Exposed time", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SHUTERR", &shuttererror, (char*)"Shutter error", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TIMEOFF", &timeoffset, (char*)"Time offset", &status);
    fits_write_key(faptr, TLONG, (char*)"FILTNM", &filter, (char*)"Filter (0=u, 1=g, 2=r, 3=i, 4=z, 5=y)", &status);
    if (filter == 0) {
        sprintf(tempstring, "u");
        fits_write_key(faptr, TSTRING, (char*)"FILTER", &tempstring, (char*)"Filter", &status);
    }
    if (filter == 1) {
        sprintf(tempstring, "g");
        fits_write_key(faptr, TSTRING, (char*)"FILTER", &tempstring, (char*)"Filter", &status);
    }
    if (filter == 2) {
        sprintf(tempstring, "r");
        fits_write_key(faptr, TSTRING, (char*)"FILTER", &tempstring, (char*)"Filter", &status);
    }
    if (filter == 3) {
        sprintf(tempstring, "i");
        fits_write_key(faptr, TSTRING, (char*)"FILTER", &tempstring, (char*)"Filter", &status);
    }
    if (filter == 4) {
        sprintf(tempstring, "z");
        fits_write_key(faptr, TSTRING, (char*)"FILTER", &tempstring, (char*)"Filter", &status);
    }
    if (filter == 5) {
        sprintf(tempstring, "y");
        fits_write_key(faptr, TSTRING, (char*)"FILTER", &tempstring, (char*)"Filter", &status);
    }
    fits_write_key(faptr, TLONG, (char*)"SEED", &ranseed, (char*)"Random seed", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRA", &pra, (char*)"Pointing RA (radians)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PDEC", &pdec, (char*)"Pointing Dec (radians)", &status);
    tempf1 = pra*180/M_PI;
    fits_write_key(faptr, TDOUBLE, (char*)"RA_DEG", &tempf1, (char*)"Pointing RA (decimal degrees)", &status);
    tempf1 = pdec*180/M_PI;
    fits_write_key(faptr, TDOUBLE, (char*)"DEC_DEG", &tempf1, (char*)"Pointing Dec (decimal degrees)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"AIRMASS", &airmass, (char*)"Airmass", &status);


    tempf1 = -rotationangle*180/M_PI;
    fits_write_key(faptr, TDOUBLE, (char*)"ROTANG", &tempf1, (char*)"Rotation Angle (rotSkyPos) (degrees)", &status);
    tempf1 = spiderangle*180/M_PI;
    fits_write_key(faptr, TDOUBLE, (char*)"SPIDANG", &tempf1, (char*)"Angle of Spider (rotTelPos)", &status);
    tempf1 = zenith*180/M_PI;
    fits_write_key(faptr, TDOUBLE, (char*)"ZENITH", &tempf1, (char*)"Zenith Angle (degrees)", &status);
    tempf1 = azimuth*180/M_PI;
    fits_write_key(faptr, TDOUBLE, (char*)"AZIMUTH", &tempf1, (char*)"Azimuthal Angle (degrees)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"ROTJITT", &rotationjitter, (char*)"Rotation Jitter (arcseconds)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"ELEJITT", &elevationjitter, (char*)"Elevation Jitter (arcseconds)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"AZIJITT", &azimuthjitter, (char*)"Azimuthal Jitter (arcseconds)", &status);
    fits_write_key(faptr, TSTRING, (char*)"TRKFILE", (char*)trackingfile.c_str(), (char*)"Tracking File", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"GNDLEVL", &groundlevel, (char*)"Ground level (m)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"XTELLOC", &xtelloc, (char*)"X telescope location (m)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"YTELLOC", &ytelloc, (char*)"Y telescope location (m)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ0", &izernike[0][0], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ1", &izernike[0][1], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ2", &izernike[0][2], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ3", &izernike[0][3], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ4", &izernike[0][4], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ5", &izernike[0][5], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ6", &izernike[0][6], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ7", &izernike[0][7], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ8", &izernike[0][8], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ9", &izernike[0][9], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ10", &izernike[0][10], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ11", &izernike[0][11], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ12", &izernike[0][12], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ13", &izernike[0][13], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ14", &izernike[0][14], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ15", &izernike[0][15], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ16", &izernike[0][16], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ17", &izernike[0][17], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ18", &izernike[0][18], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ19", &izernike[0][19], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ20", &izernike[0][20], (char*)"Primary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRIZ21", &izernike[0][21], (char*)"Primary zernike amplitude", &status);

    fits_write_key(faptr, TDOUBLE, (char*)"SECZ0", &izernike[1][0], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ1", &izernike[1][1], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ2", &izernike[1][2], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ3", &izernike[1][3], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ4", &izernike[1][4], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ5", &izernike[1][5], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ6", &izernike[1][6], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ7", &izernike[1][7], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ8", &izernike[1][8], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ9", &izernike[1][9], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ10", &izernike[1][10], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ11", &izernike[1][11], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ12", &izernike[1][12], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ13", &izernike[1][13], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ14", &izernike[1][14], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ15", &izernike[1][15], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ16", &izernike[1][16], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ17", &izernike[1][17], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ18", &izernike[1][18], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ19", &izernike[1][19], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ20", &izernike[1][20], (char*)"Secondary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"SECZ21", &izernike[1][21], (char*)"Secondary zernike amplitude", &status);

    fits_write_key(faptr, TDOUBLE, (char*)"TERZ0", &izernike[2][0], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ1", &izernike[2][1], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ2", &izernike[2][2], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ3", &izernike[2][3], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ4", &izernike[2][4], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ5", &izernike[2][5], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ6", &izernike[2][6], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ7", &izernike[2][7], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ8", &izernike[2][8], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ9", &izernike[2][9], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ10", &izernike[2][10], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ11", &izernike[2][11], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ12", &izernike[2][12], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ13", &izernike[2][13], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ14", &izernike[2][14], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ15", &izernike[2][15], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ16", &izernike[2][16], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ17", &izernike[2][17], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ18", &izernike[2][18], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ19", &izernike[2][19], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ20", &izernike[2][20], (char*)"Tertiary zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TERZ21", &izernike[2][21], (char*)"Tertiary zernike amplitude", &status);

    fits_write_key(faptr, TDOUBLE, (char*)"DETZ0", &izernike[3][0], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ1", &izernike[3][1], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ2", &izernike[3][2], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ3", &izernike[3][3], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ4", &izernike[3][4], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ5", &izernike[3][5], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ6", &izernike[3][6], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ7", &izernike[3][7], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ8", &izernike[3][8], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ9", &izernike[3][9], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ10", &izernike[3][10], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ11", &izernike[3][11], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ12", &izernike[3][12], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ13", &izernike[3][13], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ14", &izernike[3][14], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ15", &izernike[3][15], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ16", &izernike[3][16], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ17", &izernike[3][17], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ18", &izernike[3][18], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ19", &izernike[3][19], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ20", &izernike[3][20], (char*)"Detector zernike amplitude", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DETZ21", &izernike[3][21], (char*)"Detector zernike amplitude", &status);

    fits_write_key(faptr, TDOUBLE, (char*)"BOD00", &body[0][0], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD10", &body[1][0], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD20", &body[2][0], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD30", &body[3][0], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD40", &body[4][0], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD50", &body[5][0], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD60", &body[6][0], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD70", &body[7][0], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD80", &body[8][0], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD90", &body[9][0], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD100", &body[10][0], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD110", &body[11][0], (char*)"Body motion misalignment", &status);

    fits_write_key(faptr, TDOUBLE, (char*)"BOD01", &body[0][1], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD11", &body[1][1], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD21", &body[2][1], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD31", &body[3][1], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD41", &body[4][1], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD51", &body[5][1], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD61", &body[6][1], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD71", &body[7][1], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD81", &body[8][1], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD91", &body[9][1], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD101", &body[10][1], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD111", &body[11][1], (char*)"Body motion misalignment", &status);

    fits_write_key(faptr, TDOUBLE, (char*)"BOD02", &body[0][2], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD12", &body[1][2], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD22", &body[2][2], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD32", &body[3][2], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD42", &body[4][2], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD52", &body[5][2], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD62", &body[6][2], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD72", &body[7][2], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD82", &body[8][2], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD92", &body[9][2], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD102", &body[10][2], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD112", &body[11][2], (char*)"Body motion misalignment", &status);

    fits_write_key(faptr, TDOUBLE, (char*)"BOD03", &body[0][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD13", &body[1][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD23", &body[2][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD33", &body[3][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD43", &body[4][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD53", &body[5][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD63", &body[6][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD73", &body[7][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD83", &body[8][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD93", &body[9][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD103", &body[10][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD113", &body[11][3], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD04", &body[0][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD14", &body[1][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD24", &body[2][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD34", &body[3][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD44", &body[4][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD54", &body[5][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD64", &body[6][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD74", &body[7][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD84", &body[8][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD94", &body[9][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD104", &body[10][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD114", &body[11][4], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD05", &body[0][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD15", &body[1][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD25", &body[2][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD35", &body[3][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD45", &body[4][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD55", &body[5][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD65", &body[6][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD75", &body[7][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD85", &body[8][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD95", &body[9][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD105", &body[10][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"BOD115", &body[11][5], (char*)"Body motion misalignment", &status);
    fits_write_key(faptr, TLONG, (char*)"ATMFILE", &natmospherefile, (char*)"Number of atmospherefile", &status);
    for (i=0;i<natmospherefile;i++) {
        sprintf(tempstring, "AFILE%ld", i);
        fits_write_key(faptr, TSTRING, tempstring, (char*)atmospherefile[i].c_str(), (char*)"Atmosphere File", &status);
        sprintf(tempstring, "CFILE%ld", i);
        fits_write_key(faptr, TSTRING, tempstring, (char*)cloudfile[i].c_str(), (char*)"Cloud File", &status);
        sprintf(tempstring, "SEE%ld", i);
        tempf1=seefactor[i]/(pow(1/cos(zenith), 0.6));
        fits_write_key(faptr, TDOUBLE, tempstring, &tempf1, (char*)"Seeing at 5000 angstrom (sigma)", &status);
        sprintf(tempstring, "WIND%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &wind[i], (char*)"Wind speed (m/s)", &status);
        sprintf(tempstring, "WDIR%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &winddir[i], (char*)"Wind direction (degrees)", &status);
        sprintf(tempstring, "OSCL%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &outerscale[i], (char*)"Outer scale (m)", &status);
        sprintf(tempstring, "HGT%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &height[i], (char*)"Height (km)", &status);
        sprintf(tempstring, "DFLU%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &densityfluctuation[i], (char*)"Density fluctuation", &status);
        sprintf(tempstring, "DMEA%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &densitymean[i], (char*)"Density mean", &status);
        sprintf(tempstring, "CMEAN%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &cloudmean[i], (char*)"Mean cloud extinction (mag)", &status);
        sprintf(tempstring, "CVARY%ld", i);
        fits_write_key(faptr, TDOUBLE, tempstring, &cloudvary[i], (char*)"Variation of cloud ext. (mag)", &status);
    }
    fits_write_key(faptr, TLONG, (char*)"ATMDISP", &atmospheric_dispersion, (char*)"Atmos. Dispersion (0=off/1=on)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PRESS", &pressure, (char*)"Air pressure (mmHg)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"H2OPRESS", &water_pressure, (char*)"Water vapor pressure (mmHg)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"TEMPERA", &temperature, (char*)"Ground Temperature (degrees C)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"RELDENS", &raynorm, (char*)"Relative density", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"RELO2", &o2norm, (char*)"Relative O2 fraction", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"RELH2O", &h2onorm, (char*)"Relative H2O fraction", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"AERTAU", &aerosoltau, (char*)"Aerosol optical depth", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"AERIND", &aerosolindex, (char*)"Aerosol index", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"RELO3", &o3norm, (char*)"Relative O3 fraction", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"LASCPR", &miescatter_scat, (char*)"Large angle scattering probability", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"DOMESEE", &domeseeing, (char*)"Dome Seeing (arcseconds FWHM)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"PIXSIZE", &pixsize, (char*)"Pixel Size (microns)", &status);
    fits_write_key(faptr, TSTRING, (char*)"CHIPID", (char*)chipid.c_str(), (char*)"Chip/Amplifier ID", &status);
    fits_write_key(faptr, TLONG, (char*)"PAIRID", &pairid, (char*)"Pair ID", &status);
    fits_write_key(faptr, TSTRING, (char*)"FPFILE", (char*)focalplanefile.c_str(), (char*)"Focal plane file name", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"CENTX", &centerx, (char*)"Chip center x (microns)", &status);
    fits_write_key(faptr, TDOUBLE, (char*)"CENTY", &centery, (char*)"Chip center y (microns)", &status);
    fits_write_key(faptr, TLONG, (char*)"PIXX", &pixelsx, (char*)"Chip x pixels", &status);
    fits_write_key(faptr, TLONG, (char*)"PIXY", &pixelsy, (char*)"Chip y pixels", &status);
    fits_write_key(faptr, TLONG, (char*)"MINX", &minx, (char*)"Minimum x pixel of amplifier", &status);
    fits_write_key(faptr, TLONG, (char*)"MAXX", &maxx, (char*)"Maximum x pixel of amplifier", &status);
    fits_write_key(faptr, TLONG, (char*)"MINY", &miny, (char*)"Minimum y pixel of amplifier", &status);
    fits_write_key(faptr, TLONG, (char*)"MAXY", &maxy, (char*)"Maximum y pixel of amplifier", &status);
    fits_write_key(faptr, TLONG, (char*)"WELDPT", &well_depth, (char*)"Full well depth (electrons)", &status);
    fits_write_key(faptr, TFLOAT, (char*)"NBULK", &nbulk, (char*)"Bulk doping density", &status);
    fits_write_key(faptr, TFLOAT, (char*)"NF", &nf, (char*)"Front side doping density", &status);
    fits_write_key(faptr, TFLOAT, (char*)"NB", &nb, (char*)"Back side doping density", &status);
    fits_write_key(faptr, TFLOAT, (char*)"SF", &sf, (char*)"Front side doping scale", &status);
    fits_write_key(faptr, TFLOAT, (char*)"SB", &sb, (char*)"Back side doping scale", &status);
    fits_write_key(faptr, TFLOAT, (char*)"SITHICK", &siliconthickness, (char*)"Silicon Thickness (microns)", &status);
    fits_write_key(faptr, TFLOAT, (char*)"OVRDEP", &overdepbias, (char*)"Over depletion bias (volts)", &status);
    fits_write_key(faptr, TFLOAT, (char*)"CCDTEMP", &ccdtemp, (char*)"CCD temperature (K):", &status);
    fits_write_key(faptr, TFLOAT, (char*)"QEVAR", &qevariation, (char*)"QE variation:", &status);

    return(0);
}


int Observation::filterTruncateSources () {

    long i, j;
    double filterlow,  filterhigh;
    long lowptr, highptr;
    double lowvalue, highvalue;

    filterlow = 0; filterhigh = 12000;

    for (i = 0; i < nsedptr; i++) {
        highvalue = 0; lowvalue = 0; lowptr = 0; highptr = 0;
        for (j = 0; j < sed_n[i]; j++) {

            if (*(sed_w + sed_ptr[i] + j) < filterlow) {
                lowptr = j;
                lowvalue = *(sed_c + sed_ptr[i] + j);
            }
            if (*(sed_w+sed_ptr[i]+j) < filterhigh) {
                highptr = j;
                highvalue = *(sed_c + sed_ptr[i] + j);
            }

        }

        sed_corr[i]=(*(sed_c + sed_ptr[i] + highptr)-*(sed_c + sed_ptr[i] + lowptr));
        for (j = lowptr; j < highptr + 1; j++) {
            *(sed_c + sed_ptr[i] + j)=(*(sed_c + sed_ptr[i] + j) - lowvalue)/
                (highvalue - lowvalue);
        }
        sed_ptr[i] = sed_ptr[i] + lowptr;
        sed_n[i] = highptr - lowptr;
    }

    for (i = 0; i < nsource; i++) {
        sources.norm[i] = sources.norm[i]*sed_corr[sources.sedptr[i]];
    }

    return(0);

}
