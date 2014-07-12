///
/// @package phosim
/// @file image.cpp
/// @brief image class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <stdlib.h>
#include <stddef.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex>

#include "basic_types.h"
#include "constants.h"
#include "rng_mwc.h"
using namespace RandomNumbers;
#include "helpers.h"
#include "image.h"

#include "atmospheresetup.cpp"
#include "telescopesetup.cpp"
#include "photonmanipulate.cpp"
#include "photonloop.cpp"
#include "cosmicrays.cpp"

void Image::writeImageFile () {

    int status;
    char tempstring[4096], tempstring2[4096], line[4096];
    long naxesa[2];
    float tempfl;
    fitsfile *faptr;
    FILE *indafile;

    status = 0;
    std::string filename = "!"+outputdir+"/"+outputfilename+".fits.gz";
    fits_create_file(&faptr, filename.c_str(), &status);
    naxesa[0] = 1; naxesa[1] = 1;
    fits_create_img(faptr, FLOAT_IMG, 2, naxesa, &status);
    sprintf(tempstring, "RA---TAN");
    fits_write_key(faptr, TSTRING, (char*)"CTYPE1", tempstring, (char*)"", &status);
    tempfl = (float)(-centerx/pixsize+pixelsx/2);
    fits_write_key(faptr, TFLOAT, (char*)"CRPIX1", &tempfl, (char*)"", &status);
    tempfl = (float)(pra*180.0/M_PI);
    fits_write_key(faptr, TFLOAT, (char*)"CRVAL1", &tempfl, (char*)"", &status);
    sprintf(tempstring, "DEC--TAN");
    fits_write_key(faptr, TSTRING, (char*)"CTYPE2", tempstring, (char*)"", &status);
    tempfl = (float)(-centery/pixsize+pixelsy/2);
    fits_write_key(faptr, TFLOAT, (char*)"CRPIX2", &tempfl, (char*)"", &status);
    tempfl = (float)(pdec*180.0/M_PI);
    fits_write_key(faptr, TFLOAT, (char*)"CRVAL2", &tempfl, (char*)"", &status);
    tempfl = (float)(180./M_PI*(9.69627e-7)*cos(rotationangle));
    fits_write_key(faptr, TFLOAT, (char*)"CD1_1", &tempfl, (char*)"", &status);
    tempfl = (float)(180./M_PI*(9.69627e-7)*sin(rotationangle));
    fits_write_key(faptr, TFLOAT, (char*)"CD1_2", &tempfl, (char*)"", &status);
    tempfl = (float)(180./M_PI*(-9.69627e-7)*sin(rotationangle));
    fits_write_key(faptr, TFLOAT, (char*)"CD2_1", &tempfl, (char*)"", &status);
    tempfl = (float)(180./M_PI*9.69627e-7*cos(rotationangle));
    fits_write_key(faptr, TFLOAT, (char*)"CD2_2", &tempfl, (char*)"", &status);
    sprintf(tempstring, "ICRS");
    fits_write_key(faptr, TSTRING, (char*)"RADESYS", tempstring, (char*)"", &status);
    tempfl = 2000.0;
    fits_write_key(faptr, TFLOAT, (char*)"EQUINOX", &tempfl, (char*)"", &status);
    tempfl = 2000.0;
    fits_write_key(faptr, TFLOAT, (char*)"EPOCH", &tempfl, (char*)"", &status);
    header(faptr);
    sprintf(tempstring, "PHOSIM");
    fits_write_key(faptr, TSTRING, (char*)"CREATOR", tempstring, (char*)"", &status);
    sprintf(tempstring, "%s/version", bindir.c_str());
    indafile = fopen(tempstring, "r");
    fgets(line, 4096, indafile);
    sscanf(line, "%s %s", tempstring, tempstring2);
    fits_write_key(faptr, TSTRING, (char*)"VERSION", tempstring2, (char*)"", &status);
    fgets(line, 4096, indafile);
    sscanf(line, "%s", tempstring2);
    fclose(indafile);
    fits_write_key(faptr, TSTRING, (char*)"BRANCH", tempstring2, (char*)"", &status);
    if (date) fits_write_date(faptr, &status);
    naxesa[0] = chip.nampx; naxesa[1] = chip.nampy;
    fits_update_key(faptr, TLONG, (char*)"NAXIS1", &naxesa[0], NULL, &status);
    fits_update_key(faptr, TLONG, (char*)"NAXIS2", &naxesa[1], NULL, &status);
    fits_write_img(faptr, TFLOAT, 1, chip.nampx*chip.nampy, chip.focal_plane, &status);
    fits_close_file(faptr, &status);
    uint32 z = 0, w = 0;
    RngGetSeed(&z, &w);



}

void Image::writeOPD () {

    int status;
    char tempstring[4096];
    long naxesa[2];
    fitsfile *faptr;

    status = 0;
    double total = 0.0;
    double count = 0.0;


    for (long i = 0;i<SCREEN_SIZE/4;i++) {
        for (long j = 0;j<SCREEN_SIZE/4;j++) {
            if (*(opdcount+i*SCREEN_SIZE/4+j)>0) {

                *(opd+i*SCREEN_SIZE/4+j) = (*(opd+i*SCREEN_SIZE/4+j))/(*(opdcount+i*SCREEN_SIZE/4+j));
                total += *(opd+i*SCREEN_SIZE/4+j);
                count+=1.0;
            } else {
                *(opd+i*SCREEN_SIZE/4+j) = 0.0;
            }
        }
    }
    total = total/count;
    for (long i = 0;i<SCREEN_SIZE/4;i++) {
        for (long j = 0;j<SCREEN_SIZE/4;j++) {
            if (*(opdcount+i*SCREEN_SIZE/4+j)>0) {
                *(opd+i*SCREEN_SIZE/4+j) = *(opd+i*SCREEN_SIZE/4+j)-total;
            }
        }
    }

    sprintf(tempstring, "!%s/opd.fits", outputdir.c_str());
    fits_create_file(&faptr, tempstring, &status);
    naxesa[0] = 1; naxesa[1] = 1;
    fits_create_img(faptr, DOUBLE_IMG, 2, naxesa, &status);
    naxesa[0] = SCREEN_SIZE/4; naxesa[1] = SCREEN_SIZE/4;
    fits_update_key(faptr, TLONG, (char*)"NAXIS1", &naxesa[0], NULL, &status);
    fits_update_key(faptr, TLONG, (char*)"NAXIS2", &naxesa[1], NULL, &status);
    fits_write_img(faptr, TDOUBLE, 1, SCREEN_SIZE/4*SCREEN_SIZE/4, opd, &status);
    fits_close_file(faptr, &status);

}

void Image::writeCheckpoint() {

    fitsfile *faptr;
    long naxes[2];
    int status;
    double *tempDynamicTransmission;
    uint32 z = 0, w = 0;

    RngGetSeed(&z, &w);

    tempDynamicTransmission = (double*)calloc((natmospherefile*2+nsurf*2+2)*(901), sizeof(double));
    for (int i = 0;i<(natmospherefile*2+nsurf*2+2);i++) {
        for (int j = 0;j<901;j++) {
            *(tempDynamicTransmission+i*901+j) = dynamicTransmission[i][j];
        }
    }

    status = 0;
    std::string filename = "!"+outputdir+"/"+outputfilename+"_ckptdt.fits.gz";
    fits_create_file(&faptr, filename.c_str(), &status);
    naxes[0] = 1; naxes[1] = 1;
    fits_create_img(faptr, DOUBLE_IMG, 2, naxes, &status);
    naxes[0] = 901; naxes[1] = (natmospherefile*2+nsurf*2+2);
    fits_update_key(faptr, TLONG, (char*)"NAXIS1", &naxes[0], NULL, &status);
    fits_update_key(faptr, TLONG, (char*)"NAXIS2", &naxes[1], NULL, &status);
    fits_update_key (faptr, TUINT, (char*)"M_Z",  &z, NULL, &status);
    fits_update_key (faptr, TUINT, (char*)"M_W",  &w, NULL, &status);
    fits_write_img(faptr, TDOUBLE, 1, (natmospherefile*2+nsurf*2+2)*(901), tempDynamicTransmission, &status);
    fits_close_file(faptr, &status);

    free(tempDynamicTransmission);

    filename = "!"+outputdir+"/"+outputfilename+"_ckptfp.fits.gz";
    fits_create_file(&faptr, filename.c_str(), &status);
    naxes[0] = 1; naxes[1] = 1;
    fits_create_img(faptr, FLOAT_IMG, 2, naxes, &status);
    naxes[0] = chip.nampx; naxes[1] = chip.nampy;
    fits_update_key(faptr, TLONG, (char*)"NAXIS1", &naxes[0], NULL, &status);
    fits_update_key(faptr, TLONG, (char*)"NAXIS2", &naxes[1], NULL, &status);
    fits_write_img(faptr, TFLOAT, 1, chip.nampx*chip.nampy, chip.focal_plane, &status);
    fits_close_file(faptr, &status);


}

void Image::readCheckpoint() {

    fitsfile *faptr;
    long naxes[2];
    int nfound;
    int anynull;
    float nullval;
    int status;
    double *tempDynamicTransmission;
    uint32 z = 0, w = 0;

    tempDynamicTransmission = (double*)calloc((natmospherefile*2+nsurf*2+2)*(901), sizeof(double));

    std::string filename = outputdir+"/"+outputfilename+"_ckptdt.fits.gz";
    status = 0;
    if (fits_open_file(&faptr, filename.c_str(), READONLY, &status)) {printf("Error opening %s\n", filename.c_str()); exit(1);}
    fits_read_keys_lng(faptr, (char*)"NAXIS", 1, 2, naxes, &nfound, &status);
    fits_read_key(faptr, TUINT, "M_Z", &z, NULL, &status);
    fits_read_key(faptr, TUINT, "M_W", &w, NULL, &status);
    fits_read_img(faptr, TDOUBLE, 1, naxes[0]*naxes[1], &nullval, tempDynamicTransmission, &anynull, &status);
    fits_close_file(faptr, &status);
    set_seed_rng32Mwc(z,  w);

    for (int i = 0;i<(natmospherefile*2+nsurf*2+2);i++) {
        for (int j = 0;j<901;j++) {
            dynamicTransmission[i][j] = *(tempDynamicTransmission+i*901+j);
        }
    }

    free(tempDynamicTransmission);


    filename = outputdir+"/"+outputfilename+"_ckptfp.fits.gz";
    status = 0;
    if (fits_open_file(&faptr, filename.c_str(), READONLY, &status)) {printf("Error opening %s\n", filename.c_str()); exit(1);}
    fits_read_keys_lng(faptr, (char*)"NAXIS", 1, 2, naxes, &nfound, &status);
    fits_read_img(faptr, TFLOAT, 1, naxes[0]*naxes[1], &nullval, chip.focal_plane, &anynull, &status);
    fits_close_file(faptr, &status);

}


void Image::cleanup () {


}
