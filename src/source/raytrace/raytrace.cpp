///
/// @package phosim
/// @file raytrace.cpp
/// @brief raytrace functions
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @brief Modified by:
/// @author J. Garrett Jernigan (UCB)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <math.h>
#include "raytrace.h"
#include "constants.h"
#include "helpers.h"

void vectorCopy (Vector vectorIn, Vector *vectorOut) {

    vectorOut->x = vectorIn.x;
    vectorOut->y = vectorIn.y;
    vectorOut->z = vectorIn.z;

}

void vectorAdd (Vector vectorA, Vector vectorB, Vector *vectorOut) {

    vectorOut->x = vectorA.x + vectorB.x;
    vectorOut->y = vectorA.y + vectorB.y;
    vectorOut->z = vectorA.z + vectorB.z;

}

void vectorSubtract (Vector vectorA, Vector vectorB, Vector *vectorOut) {

    vectorOut->x = vectorA.x-vectorB.x;
    vectorOut->y = vectorA.y - vectorB.y;
    vectorOut->z = vectorA.z - vectorB.z;

}

void vectorInit (Vector *vector) {

    vector->x = 0.0;
    vector->y = 0.0;
    vector->z = 0.0;

}

void dcrsph (Vector *vector, double longitude, double latitude) {

    vector->x = cos(longitude)*cos(latitude);
    vector->y = sin(longitude)*cos(latitude);
    vector->z = sin(latitude);

}

void dcarsph (Vector *vector, double *longitude, double *latitude) {

    *latitude = asin(vector->z);
    *longitude = atan2(vector->y, vector->x);
    if (*longitude<0.0) *longitude += 2.0*M_PI;

}

double vectorDot(Vector *vectorA, Vector *vectorB) {

    return(vectorA->x*vectorB->x + vectorA->y*vectorB->y + vectorA->z*vectorB->z);
}

void vectorCross(Vector *vectorA, Vector *vectorB, Vector *vectorC) {

    vectorC->x = vectorA->y*vectorB->z - vectorA->z*vectorB->y;
    vectorC->y = vectorA->z*vectorB->x - vectorA->x*vectorB->z;
    vectorC->z = vectorA->x*vectorB->y - vectorA->y*vectorB->x;
}


void setup_tangent ( double ra, double dec, Vector *tpx, Vector *tpy, Vector *tpz) {

    dcrsph(tpz, ra, dec);
    dcrsph(tpx, ra + (M_PI/2.0), 0.0);
    vectorCross(tpz, tpx, tpy);
}

void tangent (double ra, double dec, double *x, double *y, Vector *tpx, Vector *tpy, Vector *tpz) {

    double dotz;
    Vector s;

    dcrsph(&s, ra, dec);
    dotz = vectorDot(&s, tpz);
    *x = vectorDot(&s, tpx)/dotz;
    *y = vectorDot(&s, tpy)/dotz;

    // s.x = tpz->x  +   (*x)*tpx->x   +  (*y)*tpy->x;
    // s.y = tpz->y  +   (*x)*tpx->y   +  (*y)*tpy->y;
    // s.z = tpz->z  +   (*x)*tpx->z   +  (*y)*tpy->z;
    // normalize(&S);

}

void shift_mu (Vector *vector, double mu, double phi) {

    double vx = vector->x;
    double vy = vector->y;
    double vz = vector->z;
    if (vx==0.0 && vy==0.0) vy = vy + 1e-9;
    double vr = sqrt(vx*vx + vy*vy);

    double x = sqrt(1 - mu*mu)*cos(phi);
    double y = sqrt(1 - mu*mu)*sin(phi);

    vector->x = (-vy)/vr*x + vz*vx/vr*y + vx*mu;
    vector->y = vx/vr*x + vz*vy/vr*y + vy*mu;
    vector->z = (-vr)*y + vz*mu;

}

void normalize (Vector *vector) {

    double r;

    r = modulus(vector);
    vector->x /= r;
    vector->y /= r;
    vector->z /= r;

}

double rotateInverseX(Vector *vector, double angle) {

    return(vector->x*cos(angle) + vector->y*sin(angle));

}

double rotateInverseY(Vector *vector, double angle) {

    return(vector->y*cos(angle) - vector->x*sin(angle));

}


double modulus (Vector *vector) {

    return(sqrt((vector->x)*(vector->x) + (vector->y)*(vector->y) + (vector->z)*(vector->z)));

}



void propagate(Vector *position, Vector angle, double distance) {

    position->x = position->x + angle.x*distance;
    position->y = position->y + angle.y*distance;
    position->z = position->z + angle.z*distance;

}

void reflect( Vector *vectorIn, Vector normal) {

    double twicevidotvn;

    twicevidotvn = 2*((vectorIn->x)*normal.x + (vectorIn->y)*normal.y + (vectorIn->z)*normal.z);

    vectorIn->x = vectorIn->x - twicevidotvn*normal.x;
    vectorIn->y = vectorIn->y - twicevidotvn*normal.y;
    vectorIn->z = vectorIn->z - twicevidotvn*normal.z;

}

void refract (Vector *vectorIn, Vector normal, double n1, double n2) {

    double vidotvn;
    double vidotvn2;
    double twicevidotvn;
    double a, b, b2, d2;

    vidotvn = (vectorIn->x)*normal.x + (vectorIn->y)*normal.y + (vectorIn->z)*normal.z;

    if (vidotvn>=0.0) { normal.x = -normal.x; normal.y = -normal.y; normal.z = -normal.z; vidotvn = -vidotvn;}

    vidotvn2 = vidotvn*vidotvn;
    b = n1/n2;
    b2 = b*b;
    d2 = 1.0 - b2*(1.0 - vidotvn2);

    if (d2 >= 0.0) {
        if (vidotvn >= 0.0) {
            a = -b*vidotvn + sqrt(d2);
            vectorIn->x = a*normal.x + b*(vectorIn->x);
            vectorIn->y = a*normal.y + b*(vectorIn->y);
            vectorIn->z = a*normal.z + b*(vectorIn->z);
        } else {
            a = -b*vidotvn-sqrt(d2);
            vectorIn->x = a*normal.x + b*(vectorIn->x);
            vectorIn->y = a*normal.y + b*(vectorIn->y);
            vectorIn->z = a*normal.z + b*(vectorIn->z);
        }
    } else {
        // Total Internal Reflection
        twicevidotvn = vidotvn + vidotvn;
        vectorIn->x = -twicevidotvn*normal.x + (vectorIn->x);
        vectorIn->y = -twicevidotvn*normal.y + (vectorIn->y);
        vectorIn->z = -twicevidotvn*normal.z + (vectorIn->z);
    }

}

void asphere (double *surface, double *surface_r, double *surface_normal, long numelements, double radiusofcurv, double spacing, double semidiameter, double holesemi, double conic, double third, double fourth, double fifth, double sixth, double seventh, double eighth, double ninth, double tenth) {

    long i;

    radiusofcurv = -radiusofcurv;
    third = third*1e3;
    fourth = fourth*1e3;
    fifth = fifth*1e3;
    sixth = sixth*1e3;
    seventh = seventh*1e3;
    eighth = eighth*1e3;
    ninth = ninth*1e3;
    tenth = tenth*1e3;

    for (i = 0;i<numelements;i++) {

        surface_r[i] = holesemi + (double)i*(semidiameter - holesemi)/((double)numelements - 1);

        if (radiusofcurv !=0 ) {

            *(surface + i) = spacing - (pow(surface_r[i], 2.0)/radiusofcurv/(1.0 + sqrt(1.0 - (conic + 1.0)*pow(surface_r[i]/radiusofcurv, 2.0))) + third*pow(surface_r[i], 3.0) + fourth*pow(surface_r[i], 4.0) + fifth*pow(surface_r[i], 5.0) + sixth*pow(surface_r[i], 6.0) + seventh*pow(surface_r[i], 7.0) + eighth*pow(surface_r[i], 8.0) + ninth*pow(surface_r[i], 9.0) + tenth*pow(surface_r[i], 10.0));
            *(surface_normal + i) =  - (surface_r[i]/radiusofcurv/sqrt(1.0 - (conic + 1.0)*pow(surface_r[i]/radiusofcurv, 2.0)) +
                                  third*pow(surface_r[i], 2.0)*3.0 +
                                  fourth*pow(surface_r[i], 3.0)*4.0 +
                                  fifth*pow(surface_r[i], 4.0)*5.0 +
                                  sixth*pow(surface_r[i], 5.0)*6.0 +
                                  seventh*pow(surface_r[i], 6.0)*7.0 +
                                  eighth*pow(surface_r[i], 7.0)*8.0 +
                                  ninth*pow(surface_r[i], 8.0)*9.0 +
                                  tenth*pow(surface_r[i], 9.0)*10.0);
        } else {
            *(surface + i) = spacing;
            *(surface_normal + i) = 0.0;
        }
    }

}


void zernikes  (double *zernike_r, double *zernike_phi, double *zernike_r_grid, double *zernike_phi_grid, double *zernike_normal_r, double *zernike_normal_phi, long numelements, long nzernikes) {

    double r, phi;
    long i;

    for (i = 0;i<numelements;i++) {

        r=((double)i+1)/((double)numelements);
        phi=((double)i)/((double)numelements)*2*M_PI;

        *(zernike_r_grid + i) = r;
        *(zernike_phi_grid + i) = phi;

        // defocus 2/0
        *(zernike_r + 0*numelements + i) = 1.0;
        *(zernike_phi + 0*numelements + i) = 1.0;
        *(zernike_normal_r + 0*numelements + i) = 0.0;
        *(zernike_normal_phi + 0*numelements + i) = 0.0;

        // defocus 2/0
        *(zernike_r + 1*numelements + i) = 2*r;
        *(zernike_phi + 1*numelements + i) = sin(phi);
        *(zernike_normal_r + 1*numelements + i) = 2.0;
        *(zernike_normal_phi + 1*numelements + i) = cos(phi);

        // defocus 2/0
        *(zernike_r + 2*numelements + i) = 2*r;
        *(zernike_phi + 2*numelements + i) = cos(phi);
        *(zernike_normal_r + 2*numelements + i) = 2.0;
        *(zernike_normal_phi + 2*numelements + i) = -sin(phi);

        // defocus 2/0
        *(zernike_r + 3*numelements + i) = sqrt(3.0)*(2*r*r - 1);
        *(zernike_phi + 3*numelements + i) = 1.0;
        *(zernike_normal_r + 3*numelements + i) = sqrt(3.0)*(4*r);
        *(zernike_normal_phi + 3*numelements + i) = 0.0;

        // astigmatism 2/2
        *(zernike_r + 4*numelements + i) = sqrt(6.0)*r*r;
        *(zernike_phi + 4*numelements + i) = cos(2*phi);
        *(zernike_normal_r + 4*numelements + i) = sqrt(6.0)*2*r;
        *(zernike_normal_phi + 4*numelements + i) = -sin(2*phi)*2;

        *(zernike_r + 5*numelements + i) = sqrt(6.0)*r*r;
        *(zernike_phi + 5*numelements + i) = sin(2*phi);
        *(zernike_normal_r + 5*numelements + i) = sqrt(6.0)*2*r;
        *(zernike_normal_phi + 5*numelements + i) = cos(2*phi)*2;

        // coma 3/1
        *(zernike_r + 6*numelements + i) = sqrt(8.0)*(3*r*r*r - 2*r);
        *(zernike_phi + 6*numelements + i) = cos(phi);
        *(zernike_normal_r + 6*numelements + i) = sqrt(8.0)*(9*r*r - 2);
        *(zernike_normal_phi + 6*numelements + i) = -sin(phi);

        *(zernike_r + 7*numelements + i) = sqrt(8.0)*(3*r*r*r - 2*r);
        *(zernike_phi + 7*numelements + i) = sin(phi);
        *(zernike_normal_r + 7*numelements + i) = sqrt(8.0)*(9*r*r - 2);
        *(zernike_normal_phi + 7*numelements + i) = cos(phi);

        // trefoil 3/3
        *(zernike_r + 8*numelements + i) = sqrt(8.0)*(r*r*r);
        *(zernike_phi + 8*numelements + i) = cos(3*phi);
        *(zernike_normal_r + 8*numelements + i) = sqrt(8.0)*(3*r*r);
        *(zernike_normal_phi + 8*numelements + i) = -sin(3*phi)*3;

        *(zernike_r + 9*numelements + i) = sqrt(8.0)*(r*r*r);
        *(zernike_phi + 9*numelements + i) = sin(3*phi);
        *(zernike_normal_r + 9*numelements + i) = sqrt(8.0)*(3*r*r);
        *(zernike_normal_phi + 9*numelements + i) = cos(3*phi)*3;

        // spherical aberration 4/0
        *(zernike_r + 10*numelements + i) = sqrt(5.0)*(6*r*r*r*r - 6*r*r + 1);
        *(zernike_phi + 10*numelements + i) = 1.0;
        *(zernike_normal_r + 10*numelements + i) = sqrt(5.0)*(24*r*r*r - 12*r);
        *(zernike_normal_phi + 10*numelements + i) = 0.0;

        // 2nd astigmatism 4/2
        *(zernike_r + 11*numelements + i) = sqrt(10.0)*(4*r*r*r*r - 3*r*r);
        *(zernike_phi + 11*numelements + i) = cos(2*phi);
        *(zernike_normal_r + 11*numelements + i) = sqrt(10.0)*(16*r*r*r - 6*r);
        *(zernike_normal_phi + 11*numelements + i) = -sin(2*phi)*2;

        *(zernike_r + 12*numelements + i) = sqrt(10.0)*(4*r*r*r*r - 3*r*r);
        *(zernike_phi + 12*numelements + i) = sin(2*phi);
        *(zernike_normal_r + 12*numelements + i) = sqrt(10.0)*(16*r*r*r - 6*r);
        *(zernike_normal_phi + 12*numelements + i) = cos(2*phi)*2;

        // quadfoil 4/4
        *(zernike_r + 13*numelements + i) = sqrt(10.0)*(r*r*r*r);
        *(zernike_phi + 13*numelements + i) = cos(4*phi);
        *(zernike_normal_r + 13*numelements + i) = sqrt(10.0)*(4*r*r*r);
        *(zernike_normal_phi + 13*numelements + i) = -sin(4*phi)*4;

        *(zernike_r + 14*numelements + i) = sqrt(10.0)*(r*r*r*r);
        *(zernike_phi + 14*numelements + i) = sin(4*phi);
        *(zernike_normal_r + 14*numelements + i) = sqrt(10.0)*(4*r*r*r);
        *(zernike_normal_phi + 14*numelements + i) = cos(4*phi)*4;

        // 2nd coma 5/1
        *(zernike_r + 15*numelements + i) = sqrt(12.0)*(10*r*r*r*r*r - 12*r*r*r + 3*r);
        *(zernike_phi + 15*numelements + i) = cos(phi);
        *(zernike_normal_r + 15*numelements + i) = sqrt(12.0)*(50*r*r*r*r - 36*r*r + 3);
        *(zernike_normal_phi + 15*numelements + i) = -sin(phi);

        *(zernike_r + 16*numelements + i) = sqrt(12.0)*(10*r*r*r*r*r - 12*r*r*r + 3*r);
        *(zernike_phi + 16*numelements + i) = sin(phi);
        *(zernike_normal_r + 16*numelements + i) = sqrt(12.0)*(50*r*r*r*r - 36*r*r + 3);
        *(zernike_normal_phi + 16*numelements + i) = cos(phi);

        // 2nd trefoil 5/3
        *(zernike_r + 17*numelements + i) = sqrt(12.0)*(5*r*r*r*r*r - 4*r*r*r);
        *(zernike_phi + 17*numelements + i) = cos(3*phi);
        *(zernike_normal_r + 17*numelements + i) = sqrt(12.0)*(25*r*r*r*r - 12*r*r);
        *(zernike_normal_phi + 17*numelements + i) = -sin(3*phi)*3;

        *(zernike_r + 18*numelements + i) = sqrt(12.0)*(5*r*r*r*r*r - 4*r*r*r);
        *(zernike_phi + 18*numelements + i) = sin(3*phi);
        *(zernike_normal_r + 18*numelements + i) = sqrt(12.0)*(25*r*r*r*r - 12*r*r);
        *(zernike_normal_phi + 18*numelements + i) = cos(3*phi)*3;

        // pentafoil 5/5
        *(zernike_r + 19*numelements + i) = sqrt(12.0)*(r*r*r*r*r);
        *(zernike_phi + 19*numelements + i) = cos(5*phi);
        *(zernike_normal_r + 19*numelements + i) = sqrt(12.0)*(5*r*r*r*r);
        *(zernike_normal_phi + 19*numelements + i) = -sin(5*phi)*5;

        *(zernike_r + 20*numelements + i) = sqrt(12.0)*(r*r*r*r*r);
        *(zernike_phi + 20*numelements + i) = sin(5*phi);
        *(zernike_normal_r + 20*numelements + i) = sqrt(12.0)*(5*r*r*r*r);
        *(zernike_normal_phi + 20*numelements + i) = cos(5*phi)*5;

    }


}

//Chebyshev polynomial of the first kind
double chebyshevT (int n, double x) {
    if (n == 0) return 1.0;
    if (n == 1) return x;
    return 2*x*chebyshevT(n-1,x)-chebyshevT(n-2,x);
}

//differentiation of Chebyshev polynomial of the first kind
double chebyshevT_x (int n, double x) {
    if (n == 0) return 0.0;
    if (n == 1) return 1.0;
    return 2*chebyshevT(n-1,x)+2*chebyshevT_x(n-1,x)-chebyshevT_x(n-2,x);
}

void chebyshevT2D (int n, double x, double y, double *t) {
    int degree, d(n+1);
    for (degree=0; degree<=n; degree++) {
        d -= (degree+1);
        if (d <= 0) {
            d += degree;
            break;
        }
    }
    double tx=chebyshevT(degree-d,x);
    double ty=chebyshevT(d,y);
    double tx_x=chebyshevT_x(degree-d,x);
    double ty_y=chebyshevT_x(d,y);
    t[0]=tx*ty;
    t[1]=tx_x*ty;
    t[2]=tx*ty_y;
}

void chebyshevs (double *r_grid, double *phi_grid, double *chebyshev, double *chebyshev_r, double *chebyshev_phi, long nPoint, long nTerm) {
    for (long j = 0; j < nPoint; j++) {
        double r = r_grid[j];
        for (long l = 0; l < nPoint; l++) {
            double x = r*cos(phi_grid[l]);
            double y = r*sin(phi_grid[l]);
            for (long i = 0; i < nTerm; i++) {
                double t[3], t_r, t_phi;
                chebyshevT2D(i,x,y,&t[0]);
                t_r=t[1]*x/r+t[2]*y/r;         // dz/dr = dz/dx cos(theta) + dz/dy sin(theta)
                t_phi=r*(-t[1]*y/r+t[2]*x/r);  // dz/dtheta = r [ -dz/dx sin(theta) + dz/dy cos(theta) ]
                chebyshev[i*nPoint*nPoint + l*nPoint + j] = t[0];
                chebyshev_r[i*nPoint*nPoint + l*nPoint + j] = t_r;
                chebyshev_phi[i*nPoint*nPoint + l*nPoint + j] = t_phi;
            }
        }
    }
}

