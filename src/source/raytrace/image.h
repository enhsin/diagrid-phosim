///
/// @package phosim
/// @file image.h
/// @brief header for image class
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <fitsio.h>
#include <fitsio2.h>
#include <fftw3.h>
#include <vector>

#include "galaxy.h"
#include "dust.h"
#include "observation.h"
#include "raytrace.h"
#include "surface.h"
#include "coating.h"
#include "silicon.h"
#include "perturbation.h"
#include "screen.h"
#include "air.h"
#include "medium.h"
#include "obstruction.h"
#include "chip.h"
#include "contamination.h"

struct Photon {
    int direction;
};

class Image : public Observation {

 public:

    // objects
    Galaxy galaxy;
    Dust dust;
    Surface surface;
    Coating coating;
    Silicon silicon;
    Air air;
    Perturbation perturbation;
    Screen screen;
    Medium medium;
    Obstruction obstruction;
    Contamination contamination;
    Chip chip;

    // global optimization data
    double *opd;
    double *opdcount;
    double **dynamicTransmission;
    int *satupmap;
    int *satdownmap;

    Photon photon;
    // photon specific-data (will turn into structure later)
    Vector interact,collect;
    double shiftedAngle;
    double airRefraction;
    double nprev, ncurr;
    double wavelength;
    double wavelengthFactor;
    long indexlx0, indexly0, indexlx1, indexly1;
    long indexcx0, indexcy0, indexcx1, indexcy1;
    long indexmx0, indexmy0, indexmx1, indexmy1;
    long indexfx0, indexfy0, indexfx1, indexfy1;
    double dlx, dly, dcx, dcy, dmx, dmy, dfx, dfy;
    long uuint, vvint, wwint;
    double windx, windy;
    long lindex;
    double xp,yp;
    double xporig, yporig, zporig;
    double opdx,opdy;
    long xPos,yPos;
    double xPosR, yPosR;
    double xpos,ypos;
    double time,prtime;
    double dvr;
    double op;
    long oindex;
    long counter;
    long maxcounter;
    int ghostFlag;
    long sourceOver_m;
    double sourceSaturationRadius;
    double saveRand[MAX_BOUNCE];

    // setup and loop
    int atmSetup();
    int telSetup();
    int photonLoop();

    // physics
    int getWavelengthTime(double *wavelength, double *time, long source);
    int domeSeeing(Vector *angle);
    int tracking(Vector *angle, double time);
    int atmosphericDispersion(Vector *angle);
    int largeAngleScattering(Vector *largeAngle);
    int secondKick(Vector *largeAngle);
    int diffraction(Vector *position, Vector angle, Vector *largeAngle);
    int samplePupil(Vector *position, long long ray);
    int transmissionCheck(double transmission, long surfaceIndex, long waveIndex);
    int transmissionPreCheck(long surfaceIndex, long waveIndex);
    int chooseSurface(long *newSurface, long *oldSurface);
    int findSurface(Vector angle, Vector position, double *distance, long surfaceIndex);
    int goldenBisectSurface(double a,double b, double c,double *z, Vector angle,
                            Vector position, double *distance, long surfaceIndex);
    int getIntercept(double x, double y, double *z, long surfaceIndex);
    int getDeltaIntercept(double x, double y, double *zv, long surfaceIndex);
    int bloom(int saturatedFlag);
    int siliconPropagate(Vector *angle,Vector *position,double lambda,Vector normal,
                         Vector *interact,Vector *collect,double dh, long waveIndex);
    int contaminationSurfaceCheck(Vector position, Vector *angle, long surfaceIndex);
    double airIndexRefraction();
    double surfaceCoating(double wavelength, Vector angle,
                          Vector normal, long surfaceIndex, double *reflection);
    double cloudOpacity(long layer);
    double atmosphereOpacity(Vector angle, long layer);
    double fringing (Vector angle, Vector normal, double wavelength, double nSi, double thickness);
    void getAngle(Vector *angle, double time, long source);
    void getDeltaAngle(Vector *angle,  Vector *position, long source);
    void newRefractionIndex(long surfaceIndex);
    void atmosphereIntercept(Vector *position, Vector angle, long layer, int mode);
    void atmosphereRefraction(Vector *angle, long layer, int mode);
    void atmosphereDiffraction(Vector *angle);
    void transform(Vector *angle, Vector *position, long surfaceIndex);
    void transformInverse(Vector *angle, Vector *position, long surfaceIndex);
    void interceptDerivatives(Vector *normal, Vector position, long surfaceIndex);
    void cosmicRays(long long *raynumber);
    void saturate(long source, Vector *largeAngle);

    // output
    void writeImageFile();
    void writeOPD();
    void readCheckpoint();
    void writeCheckpoint();
    void cleanup();

};
