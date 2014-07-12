///
/// @package phosim
/// @file dlsm.cpp
/// @brief dlsm
///
/// @brief Created by
/// @author John Peterson (Purdue)
///
/// @brief Modified by
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///
const int N=16;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "dlsm.h"
#include "raytrace/constants.h"
#include "raytrace/basic_types.h"
#include "raytrace/rng_mwc.h"


using namespace RandomNumbers;
using namespace readtext;

double asphere(double r, double radiusofcurv, double height, double conic, double third, double fourth, double fifth, double sixth, double seventh, double eighth, double ninth, double tenth) {
    double h;
    if (radiusofcurv != 0) {
        h = (height-(pow(r,2)/radiusofcurv/(1.0 + sqrt(1.0 - (conic + 1.0)*(pow(r/radiusofcurv,2))))+
                      third*pow(r,3)+fourth*pow(r,4)+fifth*pow(r,5)+
                      sixth*pow(r,6)+seventh*pow(r,7)+eighth*pow(r,8)+ninth*pow(r,9)+tenth*pow(r,10)))/1e3;
    } else {
        h = height/1e3;
    }
    return h;
}

void dlsm (int surfaceIndex, const std::string & links, long seed, const std::string & obsID, int filter, const std::string & instrdir,
           double zenith, double temperature) {

    double q2 = 0.8;
    double tol = 1e-1;

    double *x, *y, *z, *xo, *yo, *zo, *vx, *vy, *vz;
    int *e, *cx, *cy, *cz, *done;
    double *exx, *exy, *exz, *eyx, *eyy, *eyz, *ezx, *ezy, *ezz;

    x=(double*)malloc(N*N*N*sizeof(double));
    y=(double*)malloc(N*N*N*sizeof(double));
    z=(double*)malloc(N*N*N*sizeof(double));
    xo=(double*)malloc(N*N*N*sizeof(double));
    yo=(double*)malloc(N*N*N*sizeof(double));
    zo=(double*)malloc(N*N*N*sizeof(double));
    vx=(double*)malloc(N*N*N*sizeof(double));
    vy=(double*)malloc(N*N*N*sizeof(double));
    vz=(double*)malloc(N*N*N*sizeof(double));
    e=(int*)malloc(N*N*N*sizeof(int));
    cx=(int*)malloc(N*N*N*sizeof(int));
    cy=(int*)malloc(N*N*N*sizeof(int));
    cz=(int*)malloc(N*N*N*sizeof(int));
    done=(int*)malloc(N*N*N*sizeof(int));
    exx=(double*)malloc(N*N*N*sizeof(double));
    exy=(double*)malloc(N*N*N*sizeof(double));
    exz=(double*)malloc(N*N*N*sizeof(double));
    eyy=(double*)malloc(N*N*N*sizeof(double));
    eyz=(double*)malloc(N*N*N*sizeof(double));
    ezz=(double*)malloc(N*N*N*sizeof(double));

    double radiusofcurv[2];
    double height[2];
    double outerradius[2];
    double innerradius[2];
    double conic[2];
    double third[2];
    double fourth[2];
    double fifth[2];
    double sixth[2];
    double seventh[2];
    double eighth[2];
    double ninth[2];
    double tenth[2];

    if (seed == -1) RngSetSeedFromTime();
    else RngSetSeed32(seed);
    RngUnwind(10000);

    std::vector<int> linkedSurface;
    std::istringstream iss(links);
    int nn;
    while(iss>>nn) linkedSurface.push_back(nn);

    std::ostringstream opticsFile;
    opticsFile << instrdir << "/optics_" << filter << ".txt";
    readText opticsPars(opticsFile.str());
    int nsurf(0), surfIdx(0);
    double runningz(0);
    for (size_t t(0); t < opticsPars.getSize(); t++){
        std::istringstream iss(opticsPars[t]);
        std::string surfaceName, surfacetype, coatingFile, mediumFile;
        double dz, radiusCurvature;
        iss >> surfaceName >> surfacetype >> radiusCurvature >> dz;
        runningz += dz;
        //only calculate a pair of surfaces currently
        if (surfacetype != "none") {
            if (nsurf == surfaceIndex || nsurf == linkedSurface[0]) {
                radiusofcurv[surfIdx] = - radiusCurvature;
                height[surfIdx] = runningz;
                iss >> outerradius[surfIdx]
                    >> innerradius[surfIdx]
                    >> conic[surfIdx]
                    >> third[surfIdx]
                    >> fourth[surfIdx]
                    >> fifth[surfIdx]
                    >> sixth[surfIdx]
                    >> seventh[surfIdx]
                    >> eighth[surfIdx]
                    >> ninth[surfIdx]
                    >> tenth[surfIdx];
                surfIdx++;
            }
            nsurf++;
        }
    }
    printf("curv %g %g out %g %g height %g %g\n",radiusofcurv[0],radiusofcurv[1],outerradius[0],outerradius[1],height[0],height[1]);
    double height0 = (height[0] + height[1]) * 0.5;
    height[0] -= height0;
    height[1] -= height0;
    printf("height %g %g d %g mm \n",height[0],height[1],2.0*(outerradius[0])/float(N-1));

    //properties of material
    double density=2.23*(1e-3)*1e6;
    double youngs = 7.2e10;
    double nu = 0.16;
    double alpha = 0.5e-6;

    //gravity & thermal inputs
    double grav0 = 1.0;
    double grav = 9.8;
    double theta = zenith;
    double ophi = 0.0;
    double temperature0 = 20.0;
    double dtemp = temperature - temperature0;

    //derived
    double d = 2.0*(outerradius[0]/1e3)/float(N-1); //in m
    double dmm = d*1e3; // in mm
    double a3d = 15.0;
    double kn = 3.0*youngs/(1.0-2.0*nu)*d/a3d;
    double ks = 3.0*(1.0-4.0*nu)*youngs/(1.0+nu)/(1.0-2.0*nu)*d/a3d;
    double m=density*d*d*d;
    double dt=sqrt(m/(youngs*d/3.0/(1-2.0*nu)+4./3.*youngs*d/2.0/(1+nu)))*0.5;

    long kk1=-N;
    long kk2=N;
    long jj1=-N;
    long jj2=N;
    long ii1=-N;
    long ii2=N;
    long ia,ja,ka;

    for (long i=0; i<N; i++) {
        for (long j=0; j<N; j++) {
            for (long k=0; k<N; k++) {
                e[i*N*N+j*N+k]=0;
                cx[i*N*N+j*N+k]=0;
                cy[i*N*N+j*N+k]=0;
                cz[i*N*N+j*N+k]=0;
                vx[i*N*N+j*N+k]=0.0;
                vy[i*N*N+j*N+k]=0.0;
                vz[i*N*N+j*N+k]=0.0;
                double xx=-(outerradius[0]/1e3)+2*(outerradius[0]/1e3)*((double)i)/((double)(N-1));
                double yy=-(outerradius[0]/1e3)+2*(outerradius[0]/1e3)*((double)j)/((double)(N-1));
                double zz=-(outerradius[0]/1e3)+2*(outerradius[0]/1e3)*((double)k)/((double)(N-1));
                double rr=sqrt(pow(xx,2)+pow(yy,2))*1e3;
                double hh[2];
                for (int s=0; s<2; s++) {
                    hh[s] = asphere(rr,radiusofcurv[s],height[s],conic[s],third[s],fourth[s],fifth[s],
                                    sixth[s],seventh[s],eighth[s],ninth[s],tenth[s]);
                }
                //double hmax = (hh[1] > hh[0]) ? hh[1]+d/2: hh[0]+d/2;
                //double hmin = (hh[1] > hh[0]) ? hh[0]-d/2: hh[1]-d/2;
                double hmax = (hh[1] > hh[0]) ? hh[1]+d: hh[0]+d;
                double hmin = (hh[1] > hh[0]) ? hh[0]-d: hh[1]-d;
                //               if (rr <= outerradius[0]+dmm/2 && rr >= innerradius[0]-dmm/2 && zz >= hmin && zz <= hmax) {
                if (rr <= outerradius[0]+dmm && rr >= innerradius[0]-dmm && zz >= hmin && zz <= hmax) {
                    e[i*N*N+j*N+k]=1;
                    x[i*N*N+j*N+k]=xx;
                    y[i*N*N+j*N+k]=yy;
                    z[i*N*N+j*N+k]=zz;
                    xo[i*N*N+j*N+k]=xx;
                    yo[i*N*N+j*N+k]=yy;
                    zo[i*N*N+j*N+k]=zz;
                    if (k > kk1) kk1=k;
                    if (k < kk2) kk2=k;
                    if (j > jj1) jj1=j;
                    if (j < jj2) jj2=j;
                    if (i > ii1) ii1=i;
                    if (i < ii2) ii2=i;
                }
            }
        }
    }
    long ir=ii1-ii2+1;
    long jr=jj1-jj2+1;
    long kr=kk1-kk2+1;

    //constraints
    double dphi=M_PI/16.;
    double fraction=0.8;
    for (long i=0; i<N; i++) {
        for (long j=0; j<N; j++) {
            double rr=sqrt(pow(x[i*N*N+j*N+kk1],2)+pow(y[i*N*N+j*N+kk1],2));
            double phi=atan2(j-(N-1)/2,i-(N-1)/2);
            if (phi<0) phi=phi+2*M_PI;
            if (rr > outerradius[0]/1e3*fraction) {
                if (fabs(phi-0.0) < dphi || fabs(phi-8.*M_PI/12.) < dphi || fabs(phi-16.*M_PI/12.) < dphi || fabs(phi-24.*M_PI/12.) < dphi) {
                    cx[i*N*N+j*N+kk1]=1;
                    cy[i*N*N+j*N+kk1]=1;
                    cz[i*N*N+j*N+kk1]=1;
                }
            }

        }
    }

    long counter=0;
    double vpeak=0.0;
    double vmag=1.0;
    double vc;

    while (vmag >= tol*vpeak) {

        printf("%ld %lf %lf\n",counter,vmag,tol*vpeak);
        vc=0.0;
        vmag=0.0;

        for (long i=ii2; i<=ii1; i++) {
            for (long j=jj2; j<=jj1; j++) {
                for (long k=kk2; k<=kk1; k++) {

                    done[i*N*N+j*N+k]=0;

                    if (e[i*N*N+j*N+k]==1) {

                        long ip=i+1;
                        if (ip > N-1) ip=N-1;
                        if (e[ip*N*N+j*N+k] != 1) ip=i;
                        long im=i-1;
                        if (im < 0) im=0;
                        if (e[im*N*N+j*N+k] != 1) im=i;
                        long jp=j+1;
                        if (jp > N-1) jp=N-1;
                        if (e[i*N*N+jp*N+k] != 1) jp=j;
                        long jm=j-1;
                        if (jm < 0) jm=0;
                        if (e[i*N*N+jm*N+k] != 1) jm=j;
                        long kp=k+1;
                        if (kp > N-1) kp=N-1;
                        if (e[i*N*N+j*N+kp] != 1) kp=k;
                        long km=k-1;
                        if (km < 0) km=0;
                        if (e[i*N*N+j*N+km] != 1) km=k;

                        exx[i*N*N+j*N+k]=0.0;
                        exy[i*N*N+j*N+k]=0.0;
                        exz[i*N*N+j*N+k]=0.0;
                        eyy[i*N*N+j*N+k]=0.0;
                        eyz[i*N*N+j*N+k]=0.0;
                        ezz[i*N*N+j*N+k]=0.0;

                        if (im != ip) {
                            exx[i*N*N+j*N+k]+=((x[ip*N*N+j*N+k]-xo[ip*N*N+j*N+k])-(x[im*N*N+j*N+k]-xo[im*N*N+j*N+k]))/(xo[ip*N*N+j*N+k]-xo[im*N*N+j*N+k]);
                            exy[i*N*N+j*N+k]+=0.5*(((y[ip*N*N+j*N+k]-yo[ip*N*N+j*N+k])-(y[im*N*N+j*N+k]-yo[im*N*N+j*N+k]))/(xo[ip*N*N+j*N+k]-xo[im*N*N+j*N+k]));
                            exz[i*N*N+j*N+k]+=0.5*(((z[ip*N*N+j*N+k]-zo[ip*N*N+j*N+k])-(z[im*N*N+j*N+k]-zo[im*N*N+j*N+k]))/(xo[ip*N*N+j*N+k]-xo[im*N*N+j*N+k]));
                        }
                        if (jm != jp) {
                            eyy[i*N*N+j*N+k]+=((y[i*N*N+jp*N+k]-yo[i*N*N+jp*N+k])-(y[i*N*N+jm*N+k]-yo[i*N*N+jm*N+k]))/(yo[i*N*N+jp*N+k]-yo[i*N*N+jm*N+k]);
                            exy[i*N*N+j*N+k]+=0.5*(((x[i*N*N+jp*N+k]-xo[i*N*N+jp*N+k])-(x[i*N*N+jm*N+k]-xo[i*N*N+jm*N+k]))/(yo[i*N*N+jp*N+k]-yo[i*N*N+jm*N+k]));
                            eyz[i*N*N+j*N+k]+=0.5*(((z[i*N*N+jp*N+k]-zo[i*N*N+jp*N+k])-(z[i*N*N+jm*N+k]-zo[i*N*N+jm*N+k]))/(yo[i*N*N+jp*N+k]-yo[i*N*N+jm*N+k]));
                        }
                        if (km != kp) {
                            ezz[i*N*N+j*N+k]+=((z[i*N*N+j*N+kp]-zo[i*N*N+j*N+kp])-(z[i*N*N+j*N+km]-zo[i*N*N+j*N+km]))/(zo[i*N*N+j*N+kp]-zo[i*N*N+j*N+km]);
                            exz[i*N*N+j*N+k]+=0.5*(((x[i*N*N+j*N+kp]-xo[i*N*N+j*N+kp])-(x[i*N*N+j*N+km]-xo[i*N*N+j*N+km]))/(zo[i*N*N+j*N+kp]-zo[i*N*N+j*N+km]));
                            eyz[i*N*N+j*N+k]+=0.5*(((y[i*N*N+j*N+kp]-yo[i*N*N+j*N+kp])-(y[i*N*N+j*N+km]-yo[i*N*N+j*N+km]))/(zo[i*N*N+j*N+kp]-zo[i*N*N+j*N+km]));
                        }

                    }
                }
            }
        }

        for (long iii = 0; iii < (ir*jr*kr); iii++) {
        redo:;
            long i=floor(drand48()*ir)+ii2;
            long j=floor(drand48()*jr)+jj2;
            long k=floor(drand48()*kr)+kk2;

            if (done[i*N*N+j*N+k] == 1) goto redo;
            done[i*N*N+j*N+k]=1;

            if (e[i*N*N+j*N+k] == 1) {

                double forcex=0.0;
                double forcey=0.0;
                double forcez=0.0;

                double x1d=x[i*N*N+j*N+k]-xo[i*N*N+j*N+k];
                double y1d=y[i*N*N+j*N+k]-yo[i*N*N+j*N+k];
                double z1d=z[i*N*N+j*N+k]-zo[i*N*N+j*N+k];

                for (int l = 0; l < 18; l++) {

                    if (l == 0) {ia=i-1; ja=j; ka=k;}
                    if (l == 1) {ia=i+1; ja=j; ka=k;}
                    if (l == 2) {ia=i; ja=j-1; ka=k;}
                    if (l == 3) {ia=i; ja=j+1; ka=k;}
                    if (l == 4) {ia=i; ja=j; ka=k-1;}
                    if (l == 5) {ia=i; ja=j; ka=k+1;}
                    if (l == 6) {ia=i-1; ja=j+1; ka=k;}
                    if (l == 7) {ia=i-1; ja=j-1; ka=k;}
                    if (l == 8) {ia=i-1; ja=j; ka=k+1;}
                    if (l == 9) {ia=i-1; ja=j; ka=k-1;}
                    if (l == 10) {ia=i; ja=j-1; ka=k-1;}
                    if (l == 11) {ia=i; ja=j-1; ka=k+1;}
                    if (l == 12) {ia=i; ja=j+1; ka=k-1;}
                    if (l == 13) {ia=i; ja=j+1; ka=k+1;}
                    if (l == 14) {ia=i+1; ja=j+1; ka=k;}
                    if (l == 15) {ia=i+1; ja=j-1; ka=k;}
                    if (l == 16) {ia=i+1; ja=j; ka=k+1;}
                    if (l == 17) {ia=i+1; ja=j; ka=k-1;}

                    if (ia >= 0 && ja >= 0 && ka >= 0 && ia <= N-1 && ja <= N-1 && ka <= N-1) {

                        if (e[ia*N*N+ja*N+ka] == 1) {

                            double r=sqrt(pow(x[i*N*N+j*N+k]-x[ia*N*N+ja*N+ka],2)+
                                          pow(y[i*N*N+j*N+k]-y[ia*N*N+ja*N+ka],2)+pow(z[i*N*N+j*N+k]-z[ia*N*N+ja*N+ka],2));
                            double ro=sqrt(pow(xo[i*N*N+j*N+k]-xo[ia*N*N+ja*N+ka],2)+
                                           pow(yo[i*N*N+j*N+k]-yo[ia*N*N+ja*N+ka],2)+pow(zo[i*N*N+j*N+k]-zo[ia*N*N+ja*N+ka],2));

                            double nx=(x[ia*N*N+ja*N+ka]-x[i*N*N+j*N+k])/r;
                            double ny=(y[ia*N*N+ja*N+ka]-y[i*N*N+j*N+k])/r;
                            double nz=(z[ia*N*N+ja*N+ka]-z[i*N*N+j*N+k])/r;

                            double x2d=x[ia*N*N+ja*N+ka]-xo[ia*N*N+ja*N+ka];
                            double y2d=y[ia*N*N+ja*N+ka]-yo[ia*N*N+ja*N+ka];
                            double z2d=z[ia*N*N+ja*N+ka]-zo[ia*N*N+ja*N+ka];

                            double ebxx=0.5*(exx[i*N*N+j*N+k]+exx[ia*N*N+ja*N+ka]);
                            double ebxy=0.5*(exy[i*N*N+j*N+k]+exy[ia*N*N+ja*N+ka]);
                            double ebxz=0.5*(exz[i*N*N+j*N+k]+exz[ia*N*N+ja*N+ka]);
                            double ebyy=0.5*(eyy[i*N*N+j*N+k]+eyy[ia*N*N+ja*N+ka]);
                            double ebyz=0.5*(eyz[i*N*N+j*N+k]+eyz[ia*N*N+ja*N+ka]);
                            double ebzz=0.5*(ezz[i*N*N+j*N+k]+ezz[ia*N*N+ja*N+ka]);

                            double ednx=(ebxx*nx+ebxy*ny+ebxz*nz)*ro;
                            double edny=(ebxy*nx+ebyy*ny+ebyz*nz)*ro;
                            double ednz=(ebxz*nx+ebyz*ny+ebzz*nz)*ro;

                            double edndn=ednx*nx+edny*ny+ednz*nz;

                            double ux=x2d-x1d;
                            double uy=y2d-y1d;
                            double uz=z2d-z1d;
                            double udn=ux*nx+uy*ny+uz*nz;

                            forcex+=kn*udn*nx-kn*alpha*dtemp*d*nx+ks*(ednx-edndn*nx);
                            forcey+=kn*udn*ny-kn*alpha*dtemp*d*ny+ks*(edny-edndn*ny);
                            forcez+=kn*udn*nz-kn*alpha*dtemp*d*nz+ks*(ednz-edndn*nz);

                        }
                    }
                }

                forcex+=-m*grav*sin(theta)*cos(ophi)+m*grav0*grav;
                forcey+=-m*grav*sin(theta)*sin(ophi);
                forcez+=-m*grav*cos(theta);

                vx[i*N*N+j*N+k]+=forcex/m*dt-q2*fabs(forcex)/m*dt*vx[i*N*N+j*N+k]/(fabs(vx[i*N*N+j*N+k])+1e-10);
                vy[i*N*N+j*N+k]+=forcey/m*dt-q2*fabs(forcey)/m*dt*vy[i*N*N+j*N+k]/(fabs(vy[i*N*N+j*N+k])+1e-10);
                vz[i*N*N+j*N+k]+=forcez/m*dt-q2*fabs(forcez)/m*dt*vz[i*N*N+j*N+k]/(fabs(vz[i*N*N+j*N+k])+1e-10);

                if (cx[i*N*N+j*N+k] == 0) x[i*N*N+j*N+k]+=vx[i*N*N+j*N+k]*dt;
                if (cy[i*N*N+j*N+k] == 0) y[i*N*N+j*N+k]+=vy[i*N*N+j*N+k]*dt;
                if (cz[i*N*N+j*N+k] == 0) z[i*N*N+j*N+k]+=vz[i*N*N+j*N+k]*dt;

                if (cx[i*N*N+j*N+k] == 0 && cy[i*N*N+j*N+k] == 0 && cz[i*N*N+j*N+k] == 0) {
                    vmag+=sqrt(pow(vx[i*N*N+j*N+k],2)+pow(vy[i*N*N+j*N+k],2)+pow(vz[i*N*N+j*N+k],2));
                    vc+=1.0;
                }

            }
        }

        vmag=vmag/vc;
        if (vmag>vpeak) vpeak=vmag;

        counter++;

    }

    std::ostringstream surface1, surface2;
    if (height[0] > height[1]) {
        surface1 << "fea_"  << obsID << "_" << surfaceIndex << ".txt";
        surface2 << "fea_"  << obsID << "_" << linkedSurface[0] << ".txt";
    } else {
        surface1 << "fea_"  << obsID << "_" << linkedSurface[0] << ".txt";
        surface2 << "fea_"  << obsID << "_" << surfaceIndex<< ".txt";
    }

    std::ofstream output1(surface1.str().c_str());
    std::ofstream output2(surface2.str().c_str());
    for (long i=ii2; i<=ii1; i++) {
        for (long j=jj2; j<=jj1; j++) {
            for (long k=kk2; k<=kk1; k++) {
                int top(0), bottom(0);
                long idx = i*N*N+j*N+k;
                if (e[idx] == 1) {
                    if (k + 1 < N) {
                        if (e[idx+1]==0) top = 1;
                    } else top = 1;
                    if (k > 0) {
                        if (e[idx-1]==0) bottom = 1;
                    } else bottom = 1;
                    if (top == 1 && bottom == 1) {
                        std::cout<<"Warning: top = bottom\n";
                    }
                    if (top == 1 || bottom == 1) {
                        double x0 = xo[idx] * 1e3; // in mm
                        double y0 = yo[idx] * 1e3;
                        double rr = sqrt(pow(x0,2)+pow(y0,2));
                        double dx = (x[idx] - xo[idx]) * 1e3;
                        double dy = (y[idx] - yo[idx]) * 1e3;
                        double dz = (z[idx] - zo[idx]) * 1e3;
                        if (top == 1) {
                            int s = (height[0] > height[1]) ? 0:1;
                            double z0 = asphere(rr,radiusofcurv[s],height[s],conic[s],third[s],fourth[s],fifth[s],
                                                sixth[s],seventh[s],eighth[s],ninth[s],tenth[s]) * 1e3 + height0;
                            output1 << std::scientific
                                    << x0 << " " << y0 << " " << z0 << " "
                                    << dx << " " << dy << " " << dz << " 0 0 0\n";
                        }
                        if (bottom == 1) {
                            int s = (height[0] > height[1]) ? 1:0;
                            double z0 = asphere(rr,radiusofcurv[s],height[s],conic[s],third[s],fourth[s],fifth[s],
                                                sixth[s],seventh[s],eighth[s],ninth[s],tenth[s]) * 1e3 + height0;
                            output2 << std::scientific
                                    << x0 << " " << y0 << " " << z0 << " "
                                    << dx << " " << dy << " " << dz << " 0 0 0\n";
                        }
                    }
                }
            }
        }
    }
    output1.close();
    output2.close();

}

