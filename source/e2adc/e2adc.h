///
/// @package phosim
/// @file e2adc.h
/// @brief e2adc header file
///
/// @brief Created by:
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <vector>
#include <string>

class E2adc {

public:
    void setup();
    void setHotpixels();
    void convertADC();
    void writeFitsImage();

    std::string instrdir;
    std::string focalplanefile;
    std::string instr;
    std::ostringstream infile;

    //currently global
    double nonlinear;
    long well_depth;
    double parallelcte;
    double serialcte;
    int readorder;
    std::vector<double> adcerr;
    long minx, miny;
    long namp;
    long obshistid;
    long exposureid;
    std::string chipid;
    long obsseed;
    std::string devtype;
    double devvalue;
    double vistime;
    double exptime;
    long nsnap;

    //individual amp
    std::vector<std::string> outchipid;
    std::vector<long> outminx;
    std::vector<long> outmaxx;
    std::vector<long> outminy;
    std::vector<long> outmaxy;
    std::vector<long> nx;
    std::vector<long> ny;
    std::vector<int> parallelread;
    std::vector<int> serialread;
    std::vector<int> parallelPrescan;
    std::vector<int> serialOverscan;
    std::vector<int> serialPrescan;
    std::vector<int> parallelOverscan;
    std::vector<double> bias;
    std::vector<double> gain;
    std::vector<double> readnoise;
    std::vector<double> darkcurrent;
    std::vector<double> hotpixelrate;
    std::vector<double> hotcolumnrate;
    std::vector<std::vector<unsigned short> > fullReadoutMap;
    std::vector<std::vector<float> >  crosstalk;

    long onaxes[2];
    std::vector<float> emap;
    std::vector<float> adcmap;
    std::vector<unsigned int> hotpixelmap;
    int flatdir;

};
