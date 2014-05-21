#include "vector_type.h"

void vectorCopy(Vector vectorIn, Vector *vectorOut);
void setup_tangent ( double ra, double dec, Vector *tpx, Vector *tpy, Vector *tpz);
void tangent (double ra, double dec, double *x, double *y,Vector *tpx, Vector *tpy, Vector *tpz);
double vectorDot( Vector *vector1, Vector *vector2);
void vectorAdd( Vector vectorA, Vector vectorB, Vector vectorOut);
void vectorInit (Vector *vector);
void vectorSubtract( Vector vectorA, Vector vectorB, Vector vectorOut);
void normalize (Vector *vector);
double modulus (Vector *vector);
void propagate(Vector *position, Vector angle, double distance);
void reflect(Vector *vectorIn, Vector normal);
void refract (Vector *vectorIn, Vector normal, double n_1, double n_2);
double rotateInverseX(Vector *vector, double angle);
double rotateInverseY(Vector *vector, double angle);
void shift_mu (Vector *vector, double mu, double phi);
void asphere (double *surface, double *surface_r, double *surface_normal, long numelements, double radiusofcurv, double spacing, double semidiameter, double holesemi, double conic, double third, double fourth,double fifth, double sixth, double seventh, double eighth, double ninth, double tenth);
void zernikes  (double *zernike_r, double *zernike_phi, double *zernike_r_grid, double *zernike_phi_grid, double *zernike_normal_r, double *zernike_normal_phi, long numelements, long nzernikes);
