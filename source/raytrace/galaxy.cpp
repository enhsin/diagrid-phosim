///
/// @package phosim
/// @file sersic.cpp
/// @brief sersic photon sampler
///
/// @brief Created by:
/// @author Suzanne Lorenz (Purdue)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "constants.h"
#include "basic_types.h"
#include "rng_mwc.h"
using namespace RandomNumbers;
#include "helpers.h"
#include "parameters.h"

#include "galaxy.h"

double Galaxy::sample_sersic(char *sersicdata) {

    int i;
    int n;
    char temp[20];
    float b_n[99];
    double t= 0;
    double sum = 0;
    FILE *fp;

    F_all=(double**)calloc(99,sizeof(double*));
    for (n=0;n<99;n++) F_all[n]=(double*)calloc(8000,sizeof(double));
    C_all=(double**)calloc(99,sizeof(double*));
    for (n=0;n<99;n++) C_all[n]=(double*)calloc(8000,sizeof(double));
    D_all=(double*)calloc(8000,sizeof(double));

    fp = fopen(sersicdata,"r");
    for (i = 0; i < 99; i++) {
        fgets(temp,20,fp);
        sscanf(temp,"%f", &b_n[i] );
    }
    for (n = 0; n<99; n++) {
        t = 0.0;
        sum = 0.0;
        for (i = 0; i<bin_number; i++) {
            F_all[n][i] = t*t*exp(-b_n[n]*(pow(t,1/((n+1)/scale))-1));
            C_all[n][i] = F_all[n][i] + sum;
            sum = C_all[n][i];
            D_all[i] = t;
            t = t+.01;
        }
        for (i = 0; i< bin_number; i++) {
            C_all[n][i]=C_all[n][i]/sum;
        }
    }

    return 0;
}

int Galaxy::sersic(double a,double b,double c,double alpha,double beta,double n,double *x_out,double *y_out){

    double rand_az  =RngDouble();
    double theta = M_PI*2*(rand_az);
    double rand_z  =RngDouble();
    double z_axis = 2*(rand_z-.5);
    double random  =RngDouble();
    double dummy = n*scale;
    int nn = (int) (dummy);
    long q=0 ;
    find(*(C_all + nn), bin_number, random, &q);//use n here
    int down = q;
    int up = q+1;
    double frak = (random-C_all[nn][down])/(C_all[nn][up] - C_all[nn][down]);
    double r_a = a*(D_all[down]*(1-frak)+D_all[up]*frak);
    double r_b = b*(D_all[down]*(1-frak)+D_all[up]*frak);
    double r_c = c*(D_all[down]*(1-frak)+D_all[up]*frak);
    double x = sqrt(1-pow(z_axis,2))*cos(theta)*r_a;
    double y = sqrt(1-pow(z_axis,2))*sin(theta)*r_b;
    z_axis = z_axis*r_c;
    //double x_prime = x*cos(alpha)+z_axis*sin(alpha);
    double z_prime =x*sin(alpha)+z_axis*cos(alpha);
    double dist = sqrt(pow(z_prime,2.0)+pow(y,2.0));
    if (y > 0) {dist = -dist;}
    double gamma = atan(z_prime/y);
    double angle = beta+gamma;
    double f = dist*(cos(angle));
    double d = dist*(sin(angle));
    *x_out = f;
    *y_out = d;
    return 0;
}

double Galaxy::sample_sersic_2d(char *sersicdata) {
    int i;
    int n;
    char temp[20];
    float b_n[99];
    double t= 0;
    double sum = 0;
    FILE *fp;

    F_all_2d=(double**)calloc(99,sizeof(double*));
    for (n=0;n<99;n++) F_all_2d[n]=(double*)calloc(8000,sizeof(double));
    C_all_2d=(double**)calloc(99,sizeof(double*));
    for (n=0;n<99;n++) C_all_2d[n]=(double*)calloc(8000,sizeof(double));
    D_all_2d=(double*)calloc(8000,sizeof(double));

    fp = fopen(sersicdata,"r");
    for (i = 0; i < 99; i++) {
        fgets(temp,20,fp);
        sscanf(temp,"%f", &b_n[i] );
    }
    for (n = 0; n<99; n++) {
        t = 0.0;
        sum = 0.0;
        for (i = 0; i<bin_number; i++) {
            F_all_2d[n][i] = t*exp(-b_n[n]*(pow(t,1/((n+1)/scale))-1));
            C_all_2d[n][i] = F_all_2d[n][i] + sum;
            sum = C_all_2d[n][i];
            D_all_2d[i] = t;
            t = t+.01;
        }
        for (i = 0; i< bin_number; i++) {
            C_all_2d[n][i]=C_all_2d[n][i]/sum;
        }
    }

    return 0;
}

int Galaxy::sersic2d(double a,double b,double beta,double n,double *x_out,double *y_out){

    double theta=M_PI*2*(RngDouble());
    double random=RngDouble();
    int nn = (int)(n*scale);
    long q=0;
    find(*(C_all_2d + nn), bin_number, random, &q);
    double r=interpolate(D_all_2d,*(C_all_2d+nn),random,q);
    double x = cos(theta)*a*r;
    double y = sin(theta)*b*r;
    double dist = sqrt(x*x+y*y);
    if (y > 0) {dist = -dist;}
    double angle = beta+atan(x/y);
    *x_out = dist*(cos(angle));
    *y_out = dist*(sin(angle));
    return 0;
}
