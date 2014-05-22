///
/// @package phosim
/// @file atmosphere.cpp
/// @brief atmosphere parameter and screen program
///
/// @brief Created by
/// @author En-Hsin Peng (Purdue)
///
/// @brief Modified by
/// @author John R. Peterson (Purdue)
/// @author J. Chiang (SLAC)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <cstdio>
#include <cstring>
#include <iostream>

#include "ancillary/readtext.h"
#include "atmosphere/atmosphere_creator.h"
#include "atmosphere/airglow.h"
#include "atmosphere/cloud.h"
#include "atmosphere/turb2d.h"

#include "raytrace/basic_types.h"
#include "raytrace/rng_mwc.h"
using namespace RandomNumbers;

using readtext::readText;
using atmosphere::AtmosphereCreator;
using atmosphere::turb2d;
using atmosphere::cloud;
using atmosphere::airglow;

int main(void) {
    // Default values.
    int obshistid = 9999;
    double outerx = 500000.0;
    double pix = 256.0; // pixel size in cm
    int telconfig = 0;
    float constrainseeing = -1;
    int numlevel = 7;

    double seeing;
    char tempstring[4096];
    FILE *fptr;
    char outputfilename[4096];
    double zenith, altitude;
    float groundlevel;
    float monthnum;
    long seed = 0;


    static int maxlevel = 100;
    std::vector<int> cloudscreen(maxlevel, 0);
    cloudscreen[1] = 1;
    cloudscreen[2] = 1;
    std::vector<float> outerscale(maxlevel, -1);


    // Set some default values.
    std::string datadir("../data");
    std::string instrdir("../data/lsst");
    std::string outputdir(".");

    std::cout << "------------------------------------------------------------------------------------------\n"
              << "Atmosphere Creator\n"
              << "------------------------------------------------------------------------------------------" << std::endl;

    // Read parameters from stdin.
    readText pars(std::cin);
    for (size_t t(0); t < pars.getSize(); t++) {
        std::string line(pars[t]);
        readText::get(line, "outputdir", outputdir);
        readText::get(line, "instrdir", instrdir);
        readText::get(line, "datadir", datadir);
        readText::get(line, "numlevel", numlevel);
        readText::get(line, "telconfig", telconfig);
        readText::get(line, "obshistid", obshistid);
        readText::get(line, "outerx", outerx);
        readText::get(line, "cloudscreen", cloudscreen);
        readText::get(line, "outerscale", outerscale);
        readText::get(line, "obsseed", seed);
        readText::get(line, "monthnum", monthnum);

        if (readText::getKey(line, "constrainseeing", seeing)) constrainseeing = seeing/2.35482;
        else if (readText::getKey(line, "altitude", altitude)) zenith = 90-altitude;
        readText::get(line, "zenith", zenith);
    }

    readText locationPars(instrdir + "/location.txt");
    for (size_t t(0); t<locationPars.getSize(); t++) {
        std::string line(locationPars[t]);
        readText::get(line, "groundlevel", groundlevel);
    }

    sprintf(outputfilename,  "%s/atmosphere_%d.pars",  outputdir.c_str(),  obshistid);

    AtmosphereCreator creator(numlevel, groundlevel, datadir, instrdir);
    creator.run(monthnum, constrainseeing, outputfilename, cloudscreen, seed);
    const std::vector<float> & altitudes(creator.altitudes());
    std::vector<float> osests(creator.osests());

    /// overwrite outer scales if they have been provided.
    for (size_t m(0); m<osests.size(); m++) {
        if(outerscale[m]>=0) osests[m] = outerscale[m];
    }


    fptr = fopen(outputfilename, "a+");
    for (int i = 0; i < numlevel; i++) {
        printf("Creating layer %d.\n", i);
        double outers = osests[i]*100.0;
        sprintf(tempstring, "%s/atmospherescreen_%d_%d", outputdir.c_str(), obshistid, i);
        turb2d(seed*10+i, seeing, outerx, outers, zenith, 0.5, tempstring);
        fprintf(fptr, "atmospherefile %d atmospherescreen_%d_%d\n", i, obshistid, i);
    }
    for (int i = 0; i < numlevel; i++) {
        if (cloudscreen[i]) {
            double height = (altitudes[i] - groundlevel)/1000.0;
            sprintf(tempstring, "%s/cloudscreen_%d_%d", outputdir.c_str(), obshistid, i);
            cloud(seed*10+i, height, pix, tempstring);
            fprintf(fptr, "cloudfile %d cloudscreen_%d_%d\n", i, obshistid, i);
        }
    }

    sprintf(tempstring, "%s/airglowscreen_%d", outputdir.c_str(), obshistid);
    airglow(seed*10, tempstring);

    RngSetSeed32(seed);
    RngUnwind(10000);
    if (telconfig == 0) {
        fprintf(fptr, "zenith_v %.3f\n", 22.08+0.9*random_gaussian());
    } else {
        fprintf(fptr, "zenith_v 10000.0\n");
    }
    fclose(fptr);

    return 0;
}
