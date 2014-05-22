///
/// @package phosim
/// @file atmosphere_creator.cpp
/// @brief atmosphere parameter distribution calculator
///
/// @brief Created by:
/// @author Mallory Young (Purdue)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
/// @author En-Hsin Peng (Purdue)
/// @author J. Chiang (SLAC)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <cstdio>
#include <cmath>
#include <cstdlib>
#include <cctype>
#include <cstring>

#include "raytrace/basic_types.h"
#include "raytrace/rng_mwc.h"
#include "ancillary/readtext.h"

using readtext::readText;

#include "atmosphere/atmosphere_creator.h"

using namespace RandomNumbers;
// namespace atmosphere

#define pi 3.141592653589793238462643

namespace {
    float mean(float nums[], int howmany) {
        // calculates the mean of a given array of values*/
        float sum=0.0;
        int i=0;
        for(i=0;i<howmany;i++) {
            sum=sum+nums[i];
        }
        float mean=0.0;
        mean=sum/howmany;
        return mean;
    }

    float stdev( float nums[],int howmany) {
        // calculates the standard deviation of a given array of values*/
        int n=0;
        float stdv=0;
        float ave=0.0;
        float sum=0.0;
        float sumsqr=0.0;
        for(n=0;n<howmany;n++) {
            sum=sum+nums[n];
        }
        ave=sum/howmany;
        for(n=0;n<howmany;n++) {
            sumsqr=sumsqr+pow(nums[n]-ave,2);
        }
        stdv=sumsqr/howmany;
        stdv=sqrt(stdv);
        return stdv;
    }

    float interpol( float xx, float xs[], float ys[], int hm) {
        // linearly interpolates between two known data arrays*/
        float x1=0, x2=0,h1=0, h2=0, h=0;
        int i=0,j=0;
        for (i=0; i<hm; i++) {
            if ( xs[i] >= xx ) {
                x1 = xs[i];
                j = i;
            }
        }
        x2=xs[j+1];
        h1 = ys[j];
        h2 = ys[j+1];
        h = h1+((xx-x1)*(h2-h1))/(x2-x1);
        return h;
    }
} // anonymous namespace

namespace atmosphere {

    AtmosphereCreator::
    AtmosphereCreator(int numlevel, float groundlevel, const std::string & datadir, const std::string & instrdir)
        : m_numlevel(numlevel),
          m_groundlevel(groundlevel),
          m_datadir(datadir),
          m_instrdir(instrdir),
          m_osests(numlevel, 0),
          m_altitudes(numlevel, 0),
          m_jests(numlevel, 0) {}

    void AtmosphereCreator::run(float monthnum, float constrainseeing,
                        const std::string & outputfilename, const std::vector<int> & cloudscreen, long seed) {

        /// Generates estimates for the following atmospheric parameters: wind speed, wind
        /// direction, Cn^2, outer scale as a function of height
        /// Method: Please see the corresponding paper entitled "Atmospheric Parameters for the Image
        /// Simulator Based on the Site of the Large Synoptic Survey Telescope" by Young for further
        /// discussion.

        RngSetSeed32(seed);
        RngUnwind(10000);

        float totalj;
        float totalseeing;
        float bestseeing=0.0;
        std::vector<float> relseeing(m_numlevel, 0);

        m_altitudes[m_numlevel-1] = 20.0 + m_groundlevel;
        for (int i=0; i < m_numlevel-1; i++) {
            m_altitudes[i] = pow(2.0,(m_numlevel-1)-i/(((float)m_numlevel-1)/6.0))
                *(16000.0/pow(2.0,m_numlevel-1)) + m_groundlevel;
        }

        std::vector<float> magests(m_numlevel, 0);
        std::vector<float> dirests(m_numlevel, 0);
        for (int i=0; i < m_numlevel; i++) {
            magcalc(monthnum, m_altitudes[i], magests[i], dirests[i]);
            m_osests[i] = outerscale(m_altitudes[i]);
            if (m_osests[i]>(m_altitudes[i]-m_groundlevel)/2) m_osests[i]=(m_altitudes[i]-m_groundlevel)/2;
        }
        if (constrainseeing == -1) { /*program will run as normal*/
            ccalc();
            totalj = 0.0;
            for (int b=0; b < m_numlevel; b++) {
                totalj += m_jests[b];
            }
            totalseeing = 5.25/(pow(5000e-10,0.2))/(2.0*sqrt(2.0*log(2.0)))
                /(M_PI/180/3600)*pow(totalj,0.6);
            for (int b=0; b < m_numlevel; b++) {
                relseeing[b] = pow(m_jests[b]/totalj, 0.6);
            }
            float renorm = 0.0;
            for (int b=0; b < m_numlevel; b++) {
                renorm += relseeing[b]*relseeing[b];
            }
            renorm = sqrt(renorm);
            for (int b=0; b < m_numlevel; b++) {
                relseeing[b] = relseeing[b]*totalseeing/renorm;
            }
        } else if (constrainseeing>0.0) {
            int bestindex=0;
            float seeingdiff=0.0, minseeingdiff=0.0;
            std::vector<std::vector<float> > jestarray(100);
            float besttotalj,totalj;
            for (int c=0; c < 100; c++) {
                // 100 different versions of cests are created
                jestarray[c].resize(m_numlevel, 0);
                ccalc();
                for (int b=0; b < m_numlevel; b++){
                    jestarray[c][b] = m_jests[b];
                }
                totalj = 0.0;
                for (int b=0; b < m_numlevel; b++) {
                    totalj += jestarray[c][b];
                }
                totalseeing = 5.25/(pow(5000e-10,0.2))/(2.0*sqrt(2.0*log(2.0)))
                    /(M_PI/180/3600)*pow(totalj,0.6);
                // here the totalseeing value closest to constrainseeing is
                // found by comparing differences between the two
                seeingdiff = fabs(constrainseeing - totalseeing);
                if (c == 0 || seeingdiff < minseeingdiff) {
                    minseeingdiff=seeingdiff;
                    bestseeing=totalseeing;
                    bestindex=c;
                    besttotalj=totalj;
                }
            }
            totalseeing = constrainseeing;
            // the totalseeing value that is closest to constrainseeing is
            // used to generate the relseeing at each level
            for (int b = 0; b < m_numlevel; b++){
                relseeing[b] = pow(jestarray[bestindex][b]/besttotalj, 0.6);
            }
            float renorm=0.0;
            for (int b=0; b < m_numlevel; b++) {
                renorm += relseeing[b]*relseeing[b];
            }
            renorm=sqrt(renorm);
            for (int b=0; b < m_numlevel; b++) {
                relseeing[b] = relseeing[b]*constrainseeing/renorm;
            }
        } else {
            for(int c = 0; c < m_numlevel; c++) {
                relseeing[c] = 0.0;
                m_jests[c] = 0.0;
            }
            totalseeing = 0.0;
        }

        double cloudmean1, cloudsigma1, cloudmean2, cloudsigma2, cloudvariationscale;
        readText pars(m_instrdir+"/site.txt");
        for (size_t t(0); t<pars.getSize(); t++) {
            readText::get(pars[t],"cloudmean1",cloudmean1);
            readText::get(pars[t],"cloudsigma1",cloudsigma1);
            readText::get(pars[t],"cloudmean2",cloudmean2);
            readText::get(pars[t],"cloudsigma2",cloudsigma2);
            readText::get(pars[t],"cloudvariationscale",cloudvariationscale);
        }
        FILE *afile;
        if ((afile = fopen(outputfilename.c_str(),"wt"))!=NULL) {
            float value=0.0,rvalue[4];
            rvalue[0]=1.0+.01*random_gaussian();
            rvalue[1]=1.0+.01*random_gaussian();
            rvalue[2]=exp(0.18*random_gaussian());
            rvalue[3]=exp(0.20*random_gaussian());
            value=.01*sqrt(m_numlevel);
            fprintf(afile,"natmospherefile %i\n",m_numlevel);
            fprintf(afile,"totalseeing %f\n",totalseeing*(2.0*sqrt(2.0*log(2.0))));
            fprintf(afile,"reldensity %f\n",rvalue[0]);
            fprintf(afile,"relo2 %f\n",rvalue[1]);
            fprintf(afile,"relh2o %f\n",rvalue[2]);
            fprintf(afile,"relo3 %f\n\n",rvalue[3]);

            for (int i = 0; i < m_numlevel; i++) {
                fprintf(afile, "height %d %f\n", i, (m_altitudes[i]-m_groundlevel)/1000.0);
                fprintf(afile, "wind %d %f\n", i, magests[i]);
                fprintf(afile, "winddir %d %f\n", i, dirests[i]*(180.0/pi));
                fprintf(afile, "outerscale %d %f\n", i, m_osests[i]);
                fprintf(afile, "seeing %d %f\n", i, relseeing[i]);
                fprintf(afile, "densityfluc %d %f\n\n", i, value);
                if (cloudscreen[i]) {
                    float cloudmean, cloudvary;
                    if (i == 1) cloudmean = exp(log(cloudmean1)+cloudsigma1*random_gaussian());
                    if (i == 2) cloudmean = exp(log(cloudmean2)+cloudsigma2*random_gaussian());
                    cloudvary = RngDouble()*cloudmean/cloudvariationscale;
                    fprintf(afile, "cloudmean %d %f\n\n", i, cloudmean);
                    fprintf(afile, "cloudvary %d %f\n\n", i, cloudvary);
                }
            }
            fclose(afile);
        } else {
            printf("error: could not write to file.\n");
        }
    }

    void AtmosphereCreator::ccalc() {
        /// @fn void AtmosphereCreator::ccalc()
        /// @brief Interpolated form of Tokovinin and Travouillon (2006) for Cerro Pachon site
        /// Assume independent ground layer but correlated free atmosphere
        /// JRP replaces earlier model in Mallory document

        char filename[4096];
        float data[7][11];
        FILE *file;
        sprintf(filename,"%s/site.txt", m_instrdir.c_str());
        file = fopen(filename, "r");
        int i, j;
        for(i = 0; i < 7; i++) {
            fscanf(file,"%f %f %f %f %f %f %f %f %f %f %f\n",&data[i][0],
                   &data[i][1], &data[i][2], &data[i][3], &data[i][4], &data[i][5],
                   &data[i][6], &data[i][7], &data[i][8], &data[i][9], &data[i][10]);
        }
        fclose(file);
        float j1grid[7];
        float j2grid[7];
        float jtgrid[7];
        std::vector<float> hlow(m_numlevel, 0);
        std::vector<float> hhigh(m_numlevel, 0);
        float r1,r2,xl,xh,overlap;

        r1 = random_gaussian();
        r2 = random_gaussian();

        for (i=0; i < 7; i++) {

            if (data[i][4] != 0.0) {
                j1grid[i]=exp(log(data[i][4])+data[i][6]*random_gaussian());
            } else {
                j1grid[i]=0.0;
            }
            if (j1grid[i] < 0) {
                j1grid[i]=0.0;
            }
            if (data[i][8] != 0.0) {
                j2grid[i] = exp(log(data[i][8])+data[i][10]*random_gaussian());
            } else {
                j2grid[i] = 0.0;
            }
            if (j2grid[i] < 0) {
                j2grid[i] = 0.0;
            }
            jtgrid[i] = (j1grid[i] + j2grid[i])*1e-13;
        }
        for (i = 0; i < m_numlevel; i++) {
            if (i > 0) {
                hhigh[i] = 0.5*(m_altitudes[i-1] + m_altitudes[i]);
            }
            if (i < m_numlevel - 1) {
                hlow[i]=0.5*(m_altitudes[i]+m_altitudes[i+1]);
            }
        }
        hlow[m_numlevel-1]=0.0;
        hhigh[0] = 20000 + m_groundlevel;

        for (i=0; i < m_numlevel; i++) {
            m_jests[i]=0.0;
            for (j=0;j<7;j++) {
                if (data[j][0] > hlow[i]) {
                    xl=data[j][0];
                } else {
                    xl=hlow[i];
                }
                if (data[j][2] < hhigh[i]) {
                    xh=data[j][2];
                } else {
                    xh=hhigh[i];
                }
                overlap=(xh-xl)/(data[j][2]-data[j][0]);
                if (data[j][2] < hlow[i]) {
                    overlap=0.0;
                }
                if (data[j][0] > hhigh[i]) {
                    overlap=0.0;
                }
                m_jests[i] = m_jests[i] + overlap*jtgrid[j];
            }
        }
    }

    void AtmosphereCreator::magcalc(float monthnum, float altitude, float & magest, float & direst) {

        /// @fn void AtmosphereCreator::magcalc(float monthnum, float altitude, float & magest, float & direst)
        /// @brief creates wind speed and direction estimates.  Wind speed is based
        /// on a Rayleigh distribution and interpolation of known data values.
        /// The speed in given in m/s.  The wind direction is based on a Gaussian distribution and
        /// is ultimately converted to degrees, where 0 degrees points directly east.
        /// "pressure_altitude_relationship.txt" - contains corresponding pressures for each altitude
        /// "vmonthly10.txt" through "vmonthly1000.txt" - these 17 data files contain velocity data
        /// values for the v-wind componet for each pressure level.  This data is taken from the website
        /// www.noaa.gov using their NCAR Reanalysis data.
        /// "umonthly10.txt" through "umonthly1000.txt" - these files contain data for the u-component
        /// of the wind.


        char filename[4096];
        static int nsig(17);
        float sigmas[nsig],udata[13][61],mags[61],sqmags[61],meanmag[nsig];
        float uvel[61], vvel[61], stddevs[nsig], avedirs[nsig],dirs[61],
            vdata[13][61];
        float sum=0.0,sigma=0.0,u=0.0, v=0.0,dir=0.0,uave=0.0,vave=0.0;
        char* pres[]={(char*)"10",(char*)"20",(char*)"30",(char*)"50",
                      (char*)"70",(char*)"100",(char*)"150",(char*)"200",
                      (char*)"250",(char*)"300",(char*)"400",(char*)"500",
                      (char*)"600", (char*)"700",(char*)"850",(char*)"925",
                      (char*)"1000"};
        int i=0;
        for(i=0;i<nsig;i++) {
            sprintf(filename,"%s/atmosphere/vmonthly%s.txt",m_datadir.c_str(),pres[i]);
            FILE *f=fopen(filename,"r");
            int k=0;
            for(k=0; k<61; k++) {/*reads in v velocities for a given pressure level*/
                fscanf(f,"%f %f %f %f %f %f %f %f %f %f %f %f %f",&vdata[0][k],
                       &vdata[1][k],&vdata[2][k],&vdata[3][k],&vdata[4][k],
                       &vdata[5][k],&vdata[6][k],&vdata[7][k],&vdata[8][k],
                       &vdata[9][k],&vdata[10][k],&vdata[11][k],&vdata[12][k]);
            }
            fclose(f);
            sprintf(filename,"%s/atmosphere/umonthly%s.txt",m_datadir.c_str(),pres[i]);
            FILE *f2=fopen(filename,"r");
            for(k=0; k<61; k++) {/*reads in u velocities*/
                fscanf(f2,"%f %f %f %f %f %f %f %f %f %f %f %f %f",&udata[0][k],
                       &udata[1][k],&udata[2][k],&udata[3][k],&udata[4][k],
                       &udata[5][k],&udata[6][k],&udata[7][k],&udata[8][k],
                       &udata[9][k],&udata[10][k],&udata[11][k],&udata[12][k]);
            }
            fclose(f2);
            for(k=0;k<61;k++){
                // finds the appropriate velocities for the chosen month
                uvel[k]=udata[(int) monthnum][k];
                vvel[k]=vdata[(int) monthnum][k];
            }
            for(k=0;k<61;k++) {
                // finds the magnitudes of those velocity components
                mags[k]=sqrt(pow(uvel[k],2)+pow(vvel[k],2));
            }
            uave=mean(uvel,61);
            vave=mean(vvel,61);
            int j=0;
            for(j=0;j<61;j++) {
                // finds the direction in radians of each velocity component pair*/
                u=uvel[j];
                v=vvel[j];
                if(u>0.0) {
                    if(v>0.0) {
                        dirs[j]=atan(v/u);
                    } else {
                        dirs[j]=atan(v/u)+2*pi;
                    }
                } else {
                    if(v>0.0) {
                        dirs[j]=atan(v/u)+pi;
                    } else {
                        dirs[j]=atan(v/u)+pi;
                    }
                }
            }
            // finds the standard deviation of the directions
            stddevs[i]=stdev( dirs,61);
            float dirs2[61];
            float stddevs2[nsig];
            int n=0;
            for(n=0;n<61;n++) {
                dirs2[n]=dirs[n];
            }
            for(n=0;n<61;n++) {
                if(dirs[n]<= pi) {
                    dirs2[n]=dirs2[n]+2*pi;
                }
            }
            stddevs2[i]=stdev( dirs2, 61);
            if (stddevs2[i]<stddevs[i]) {
                // finds the true standard deviation by taking the lower of
                // the two ie. when taken from 0 to 2pi or pi to 3pi
                stddevs[i]=stddevs2[i];
            }
            if (uave > 0.0) {/*finds the direction of the average magnitude*/
                if (vave > 0.0) {
                    dir = atan(vave/uave);
                } else {
                    dir = atan(vave/uave) + 2*pi;
                }
            } else {
                if (vave > 0) {
                    dir = atan(vave/uave) + pi;
                } else {
                    dir = atan(vave/uave) + pi;
                }
            }
            avedirs[i]=dir;
            meanmag[i]=mean(mags,61);/*creates an array of mean magnitudes for each height*/
            int b=0;
            sum=0.0;
            for(b=0;b<61;b++) {
                sqmags[b]=pow(mags[b],2);
                sum=sum+sqmags[b];
            }
            sigma=sqrt(sum/122.0);
            sigmas[i]=sigma;
        }
        float hdata[2][nsig];
        float alts[nsig];
        FILE *hfile;
        sprintf(filename,"%s/atmosphere/pressure_altitude_relationship.txt",m_datadir.c_str());
        hfile=fopen(filename,"r");
        for(i=0; i<nsig; i++) { /*gets the altitudes for each set pressure level*/
            fscanf(hfile,"%f %f ",&hdata[0][i],&hdata[1][i]);
        }
        fclose(hfile);
        for(i=0;i<nsig;i++) {
            alts[i]=hdata[1][i];
        }

        // interpolates a sigma value for the desired height
        float mean = interpol(altitude,alts,sigmas,nsig);
        float r = RngDouble();
        float squareroot = sqrt(-2*log(r));
        magest = mean*squareroot;

        // interpolates a direction for the desired altitude
        float intdir = interpol(altitude, alts, avedirs, nsig);
        // interpolates a standard deviation for the desired altitude
        float intstd = interpol(altitude, alts, stddevs, nsig);
        r = random_gaussian();
        direst = intdir + r*intstd;
        if (direst < 0.0) {
            direst = direst+2*pi;
        }
    }

    float AtmosphereCreator::outerscale(float altitude) const {
        /// @fn float AtmosphereCreator::outerscale(float altitude) const
        //  @brief creates an outerscale estimate based on a lognormal distribution.
        /// The estimates are created from known data values taken from the 1999 Gemini Site
        /// Testing Campaign and also employs the Tartarski formulas for outerscale dependence
        /// on altitude.

        float meanos, mean, median;

        // float num = (altitude-8500.0)/2500.0;
        // The following may be over-modeling
        // so it has been removed
        //
        // if (altitude <= 2000.0) {
        //     meanos = 3.21*(1/(pow(altitude, .11)));
        // } else {
        //     if (altitude <= 17000.0) {
        //         meanos = 4.0/(1.0+(pow(num,2)));
        //     } else {
        //         meanos = 0.307-0.0324*((altitude/1000)-17)
        //             + 0.00167*pow(((altitude/1000)-17),2)
        //             + 0.000476*pow(((altitude/1000)-17),3);
        //     }
        // }
        // meanos = meanos/1.27;
        // meanos = meanos*26.7;

        readText pars(m_instrdir+"/site.txt");
        for (size_t t(0); t<pars.getSize(); t++) {
            readText::get(pars[t],"meanos",meanos);
            readText::get(pars[t],"mean",mean);
            readText::get(pars[t],"median",median);
        }

        float mu = log(median);
        float sigma = pow(2.0*(log(mean)-mu), .500);
        mu = log(meanos) - (1/2)*(pow(sigma, 2.0));
        float r = random_gaussian();
        float osest = exp(mu+sigma*r);
        return osest;
    }

} // namespace atmosphere
