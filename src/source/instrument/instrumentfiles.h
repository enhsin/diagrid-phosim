///
/// @package phosim
/// @file instrumentfiles.h
/// @brief Class to make instrument files for program insturment
///
/// @brief Created by:
/// @author Glenn Sembroski
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#ifndef INSTRUMENTFILES_H
#define INSTRUMENTFILES_H

#include <cstdlib>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
#include <sstream>
#include "../ancillary/phosim_parser.h"
#include "../ancillary/readtext.h"
#include "../raytrace/basic_types.h"
#include "../raytrace/rng_mwc.h"
#include "../raytrace/parameters.h"
using namespace RandomNumbers;
using ancillary::PhosimParser;
using readtext::readText;

class InstrumentFiles {
 public:
    InstrumentFiles();
    ~InstrumentFiles();

    void makeTrackingFile(std::string trackingFileName, double vistime, double jittertime);
    int  readControlFile(PhosimParser& pars, PhosimParser& controlPars);
    void readActuatorFile(PhosimParser& pars, PhosimParser& controlPars,
                          std::vector < std::vector < double> >& actuatorMatrix,
                          std::vector < double >& actuatorDistance,
                          std::vector < double >& actuatorError);
    void writeBodyFile(std::map< int, std::vector< double>* >& body, std::string opticsFileName, std::string obsID);
    void readoutPars(std::string focalPlaneLayoutFileName,
                     std::string segmentationFileName,
                     std::string readoutString, int camConfig);
    void focalPlanePars(std::string focalPlaneLayoutFileName,
                        std::string outChipString, int camConfig);
    void makeSurfaceMap(std::string opticsFile);
    bool getDeviceIndex(std::string deviceStr, int& deviceIndex);
 private:
    std::map< std::string, int > fSurfaceMap;
    std::map< std::string, int >::iterator fSurfaceMapPos;
    int fLastDevice;
    int fLastSurface;
    bool haveBody;
};
#endif
