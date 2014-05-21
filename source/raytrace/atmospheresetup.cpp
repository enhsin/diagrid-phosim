///
/// @package phosim
/// @file atmospheresetup.cpp
/// @brief setup for atmosphere (part of image class)
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

#include "ancillary/readtext.h"
#include "constants.h"

using readtext::readText;

int Image::atmSetup () {

    long i,j;
    char tempstring[4095];
    int status;
    long layer;
    std::string dir;

    if (flatdir==0) dir=datadir+"/atmosphere";
    else if (flatdir==1) dir=".";

    sources.sedfilename.clear();
    sources.spatialname.clear();
    sources.dustname.clear();
    sources.dustnamez.clear();

    ranseed=obsseed+pairid;

    if (obsseed==-1)
    RngSetSeedFromTime();
    else
        RngSetSeed32(ranseed);

    RngUnwind(10000);

    if (devtype=="CMOS") {
        exptime=vistime/nsnap;
        timeoffset=0.5*exptime+pairid*(exptime)-0.5*vistime;
    } else {
        exptime=(vistime-(nsnap-1)*devvalue)/nsnap;
        timeoffset=0.5*exptime+pairid*(devvalue+exptime)-0.5*vistime;
    }

    exptime=exptime+shuttererror*random_gaussian();

    // SKY SETUP
    char sersicdata[4096];
    if (flatdir==0) sprintf(sersicdata,"%s/sky/sersic_const.txt",datadir.c_str());
    else if (flatdir==1) sprintf(sersicdata,"sersic_const.txt");
    galaxy.sample_sersic(sersicdata);
    galaxy.sample_sersic_2d(sersicdata);
    dust.setup();
    filterTruncateSources();
    totalnorm=0.0;
    for (i=0;i<nsource;i++) totalnorm+=sources.norm[i];


    // ATMOSPHERIC SETUP
    for (i=0;i<natmospherefile;i++) {
        seefactor[i] = pow(1/cos(zenith),0.6)*seefactor[i];
    }

    screen.large_sizeperpixel=5120.0;
    screen.coarse_sizeperpixel=640.0;
    screen.medium_sizeperpixel=80.0;
    screen.fine_sizeperpixel=10.0;

    screen.wavelengthfactor_nom=pow(0.5,-0.2);

    // Air Setup
    fprintf(stdout,"Creating Air.\n");
    air.opacitySetup(zenith,height,groundlevel,raynorm,o2norm,h2onorm,o3norm,aerosoltau,aerosolindex,natmospherefile,dir,&airmass);

    // ATMOSPHERIC TURBULENCE
    fprintf(stdout,"Generating Turbulence.\n");
    screen.seex_large=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
    screen.seey_large=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
    screen.seex_coarse=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
    screen.seey_coarse=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
    screen.seex_medium=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
    screen.seey_medium=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();

    screen.phase_large=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
    screen.phase_coarse=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
    screen.phase_medium=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
    screen.phase_fine=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();

    screen.phaseh_medium=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();
    screen.phaseh_fine=new float[natmospherefile*SCREEN_SIZE*SCREEN_SIZE]();

    screen.phasescreen=new double[SCREEN_SIZE*SCREEN_SIZE]();
    screen.focalscreen=new double[SCREEN_SIZE*SCREEN_SIZE]();
    screen.tfocalscreen=new double[SCREEN_SIZE*SCREEN_SIZE]();
    screen.outscreen=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*SCREEN_SIZE*SCREEN_SIZE);
    screen.inscreen=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*SCREEN_SIZE*SCREEN_SIZE);
    screen.see_norm=new float[natmospherefile]();
    screen.phase_norm=new float[natmospherefile]();


    for (layer = 0; layer < natmospherefile; layer++) {


        if (cloudvary[layer]!=0 || cloudmean[layer]!=0) {
            screen.cloud[layer]=(float*)calloc(SCREEN_SIZE*SCREEN_SIZE,sizeof(float));

            sprintf(tempstring,"%s.fits",cloudfile[layer].c_str());
            {
                char *ffptr;
                fitsfile *faptr;
                long naxes[2];
                int nfound;
                int anynull;
                float nullval;

                ffptr=tempstring;
                status=0;
                if (fits_open_file(&faptr,ffptr,READONLY,&status)) { printf("Error opening %s\n",ffptr); exit(1);}
                fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
                fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.cloud[layer],&anynull,&status);
                fits_close_file(faptr,&status);
            }

        }

        sprintf(tempstring,"%s_largep.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;
            int keynum;
            char keyname[4096];
            char comment[4096];
            char value[4096];

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) { printf("Error opening %s\n",ffptr); exit(1); }
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            keynum=9;
            fits_read_keyn(faptr,keynum,keyname,value,comment,&status);
            screen.phase_norm[layer]=atof(value);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.phase_large+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }
        sprintf(tempstring,"%s_coarsep.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;
            int keynum;
            char keyname[4096];
            char comment[4096];
            char value[4096];

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) { printf("Error opening %s\n",ffptr); exit(1); }
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            keynum=9;
            fits_read_keyn(faptr,keynum,keyname,value,comment,&status);
            screen.phase_norm[layer]=atof(value);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.phase_coarse+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }
        sprintf(tempstring,"%s_mediump.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;
            int keynum;
            char keyname[4096];
            char comment[4096];
            char value[4096];

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) { printf("Error opening %s\n",ffptr); exit(1); }
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            keynum=9;
            fits_read_keyn(faptr,keynum,keyname,value,comment,&status);
            screen.phase_norm[layer]=atof(value);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.phase_medium+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }
        sprintf(tempstring,"%s_finep.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;
            int keynum;
            char keyname[4096];
            char comment[4096];
            char value[4096];

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) { printf("Error opening %s\n",ffptr); exit(1); }
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            keynum=9;
            fits_read_keyn(faptr,keynum,keyname,value,comment,&status);
            screen.phase_norm[layer]=atof(value);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.phase_fine+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }
        sprintf(tempstring,"%s_mediumh.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;
            int keynum;
            char keyname[4096];
            char comment[4096];
            char value[4096];

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) { printf("Error opening %s\n",ffptr); exit(1); }
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            keynum=9;
            fits_read_keyn(faptr,keynum,keyname,value,comment,&status);
            screen.phase_norm[layer]=atof(value);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.phaseh_medium+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }
        sprintf(tempstring,"%s_finep.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;
            int keynum;
            char keyname[4096];
            char comment[4096];
            char value[4096];

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) { printf("Error opening %s\n",ffptr); exit(1); }
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            keynum=9;
            fits_read_keyn(faptr,keynum,keyname,value,comment,&status);
            screen.phase_norm[layer]=atof(value);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.phaseh_fine+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }


        sprintf(tempstring,"%s_largex.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;
            int keynum;
            char keyname[4096];
            char comment[4096];
            char value[4096];

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) {printf("Error opening %s\n",ffptr); exit(1);}
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            keynum=9;
            fits_read_keyn(faptr,keynum,keyname,value,comment,&status);
            screen.see_norm[layer]=atof(value);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.seex_large+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }
        sprintf(tempstring,"%s_largey.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) { printf("Error opening %s\n",ffptr); exit(1);}
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.seey_large+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }

        sprintf(tempstring,"%s_coarsex.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;
            int keynum;
            char keyname[4096];
            char comment[4096];
            char value[4096];

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) {printf("Error opening %s\n",ffptr); exit(1);}
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            keynum=9;
            fits_read_keyn(faptr,keynum,keyname,value,comment,&status);
            screen.see_norm[layer]=atof(value);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.seex_coarse+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }
        sprintf(tempstring,"%s_coarsey.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) { printf("Error opening %s\n",ffptr); exit(1);}
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.seey_coarse+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }

        sprintf(tempstring,"%s_mediumx.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) {printf("Error opening %s\n",ffptr); exit(1);}
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.seex_medium+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }

        sprintf(tempstring,"%s_mediumy.fits",atmospherefile[layer].c_str());
        {
            char *ffptr;
            fitsfile *faptr;
            long naxes[2];
            int nfound;
            int anynull;
            float nullval;

            ffptr=tempstring;
            status=0;
            if (fits_open_file(&faptr,ffptr,READONLY,&status)) {printf("Error opening %s\n",ffptr); exit(1);}
            fits_read_keys_lng(faptr,(char*)"NAXIS",1,2,naxes,&nfound,&status);
            fits_read_img(faptr,TFLOAT,1,naxes[0]*naxes[1],&nullval,screen.seey_medium+layer*SCREEN_SIZE*SCREEN_SIZE,&anynull,&status);
            fits_close_file(faptr,&status);
        }



        for (i = 0; i < SCREEN_SIZE; i++) {
            for (j = 0; j < SCREEN_SIZE; j++) {
                *(screen.seex_large+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)=(float)
                    (*(screen.seex_large+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)*seefactor[layer]);
                *(screen.seey_large+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)=(float)
                    (*(screen.seey_large+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)*seefactor[layer]);
                *(screen.seex_coarse+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)=(float)
                    (*(screen.seex_coarse+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)*seefactor[layer]);
                *(screen.seey_coarse+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)=(float)
                    (*(screen.seey_coarse+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)*seefactor[layer]);
                *(screen.seex_medium+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)=(float)
                    (*(screen.seex_medium+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)*seefactor[layer]);
                *(screen.seey_medium+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)=(float)
                    (*(screen.seey_medium+layer*SCREEN_SIZE*SCREEN_SIZE+SCREEN_SIZE*i+j)*seefactor[layer]);
            }
        }




    }

    air.air_refraction_adc=64.328+29498.1/(146-1/central_wavelength/central_wavelength)+255.4/(41-1/central_wavelength/central_wavelength);
    air.air_refraction_adc=air.air_refraction_adc*pressure*(1+(1.049-0.0157*temperature)*1e-6*pressure)/720.883/(1+0.003661*temperature);
    air.air_refraction_adc=air.air_refraction_adc-((0.0624-0.000680/central_wavelength/central_wavelength)/(1+0.003661*temperature)*water_pressure);

    return(0);

}
