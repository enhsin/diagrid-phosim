///
/// @package phosim
/// @file observation.h
/// @brief observation header file
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

#include "parameters.h"
#include "vector_type.h"
#include <string>
#include <vector>

enum SourceTypes {POINT=0,IMAGE=1,GAUSSIAN=2,MOVINGPOINT=4,SERSIC=5,SERSIC2D=6,PINHOLE=7};

struct Source {
    std::vector<double> ra;
    std::vector<double> redshift;
    std::vector<double> gamma1;
    std::vector<double> gamma2;
    std::vector<double> kappa;
    std::vector<double> deltara;
    std::vector<double> deltadec;
    std::vector<double> dec;
    double *vx;
    double *vy;
    double *vz;
    double *norm;
    double *mag;
    double **spatialpar;
    double **dustpar;
    double **dustparz;
    std::vector<double> id;
    int *spatialtype;
    int *dusttype;
    int *dusttypez;
    int *type;
    std::vector<std::string> sedfilename;
    std::vector<std::string> spatialname;
    std::vector<std::string> dustname;
    std::vector<std::string> dustnamez;
    long *skysameas;
    long *sedptr;
};

class Observation {

 public:
    double rotationjitter;
    double zenith;
    double azimuth;
    double rotationrate;
    double elevationjitter;
    double azimuthjitter;
    double windjitter;
    double groundlevel;
    double xtelloc;
    double ytelloc;
    double latitude;
    double longitude;
    std::vector<double> seefactor;
    std::vector<double> wind;
    std::vector<double> winddir;
    std::vector<double> outerscale;
    std::vector<double> height;
    std::vector<double> densityfluctuation;
    std::vector<double> densitymean;
    std::vector<double> cloudmean;
    std::vector<double> cloudvary;
    std::vector<std::string> extraCommandString;
    double *dtau;
    double pressure;
    double water_pressure;
    double temperature;
    double raynorm;
    double o2norm;
    double o3norm;
    double h2onorm;
    double aerosoltau;
    double aerosolindex;
    int NZERN; // move from parameters.h, for chebyshev polynomials
    std::string pertType;
    std::vector<std::vector<double> > izernike;
    std::vector<std::vector<double> > body;
    double minr;
    double shuttererror;
    double transtol,backAlpha,backBeta,backRadius,backGamma,backDelta,
        backBuffer,np,finiteDistance;
    double screentol;
    double maxr;
    double domeseeing;
    double pixsize;
    double pra;
    double pdec;
    double spiderangle;
    double platescale;
    double centerx;
    double centery;
    double tai;
    double miescatter_scat;
    double exptime;
    double vistime;
    double timeoffset;
    double rotationangle;
    double *sed_corr;
    double *sed_dwdp;
    double *sed_w;
    double *sed_c;
    double large_scale;
    double coarse_scale;
    double medium_scale;
    double fine_scale;
    double totalnorm;
    double totalseeing;
    double moonalt;
    double moondist;
    double moonra;
    double moondec;
    double phaseang;
    double solarzen;
    double zenith_v;
    double watervar;
    double central_wavelength;
    double domelight;
    double domewave;
    double chipangle;
    double decenterx;
    double decentery;
    int telconfig;
    int checkpointcount;
    int checkpointtotal;
    int impurityvariation;
    int fieldanisotropy;
    int fringeflag;
    int deadlayer;
    int chargesharing;
    int pixelerror;
    int chargediffusion;
    int airrefraction;
    double raydensity;
    double scalenumber;
    double devvalue;
    double airglowvariation;
    float nbulk;
    float nf;
    float nb;
    float sf;
    float sb;
    float siliconthickness;
    float overdepbias;
    float siliconreflectivity;
    float ccdtemp;
    float qevariation;
    float *airglow;
    long long nphot;
    long long *source_xpos;
    long long *source_ypos;
    long long *source_photon;
    long airglowScreenSize;
    long telescope_on;
    long coatingmode;
    long contaminationmode;
    long tracking_on;
    long ranseed;
    long obsseed;
    long zernikemode;
    long atmospheric_dispersion;
    long atmosphericdispcenter;
    long natmospherefile;
    long straylight;
    double straylightcut;
    long detector_on;
    long diffraction_on;
    long aperturemode;
    long filter;
    long saturation;
    long eventfile;
    long opdfile;
    long centroidfile;
    long throughputfile;
    long pixelsx,pixelsy,minx,maxx,miny,maxy;
    long obshistid;
    long pairid;
    long blooming;
    long well_depth;
    long *sed_n;
    long *sed_ptr;
    long nsedptr;
    long sedptr;
    long nsource;
    long nimage;
    long nreallocs;
    long sed_max;
    long ghostonly;
    long atmdebug;
    long large_grid;
    long coarse_grid;
    long medium_grid;
    long fine_grid;
    long nsnap;
    std::vector<int> ghost;
    int flatdir;
    int date;
    std::string devtype;
    std::vector<int> feaflag;
    std::vector<std::string> feafile;


    /* should be part of image but have conflict with settings.c */
    long nsurf;
    long nmirror;
    long npertsurf;
    double airmass;

    /* remainder should be ok */
    std::vector<std::string> atmospherefile;
    std::vector<std::string> cloudfile;
    std::string trackingfile;
    long trackinglines;
    std::string outputfilename;
    std::string chipid;
    std::string focalplanefile;
    std::string outputdir;
    std::string seddir, imagedir;
    std::string datadir, instrdir, bindir;
    std::string eventFitsFileName;
    Vector tpx, tpy, tpz;

    long naxesb[MAX_IMAGE][2];
    float *tempptr[MAX_IMAGE];
    float *cumulativex[MAX_IMAGE];
    float cumulative[MAX_IMAGE];
    float cumulativef[MAX_IMAGE];

    Source sources;

    // functions
     int parser();
    int background();
    int addSource(const std::string & object, int sourcetype);
    int header(fitsfile *faptr);
    int settings();
    int filterTruncateSources();

};
