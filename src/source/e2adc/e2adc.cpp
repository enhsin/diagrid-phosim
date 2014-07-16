///
/// @package phosim
/// @file e2adc.cpp
/// @brief electron digitization code
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

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <iostream>

#include "e2adc.h"
#include "ancillary/readtext.h"
#include "ancillary/fits_utils.h"

#include "raytrace/basic_types.h"
#include "raytrace/rng_mwc.h"

using namespace RandomNumbers;
using namespace fits_utils;
using readtext::readText;

void E2adc::setup() {

    std::cout << "------------------------------------------------------------------------------------------" << std::endl;
    std::cout << "Electron to ADC Image Converter" << std::endl;
    std::cout << "------------------------------------------------------------------------------------------" << std::endl;

    chipid = "R22_S11";
    instrdir = "../data/lsst";
    filter = 0;
    flatdir = 0;
    tarfile = 0;
    readorder = 1;
    serialcte = 0.999995;
    parallelcte = 1.0;
    nonlinear = 0.0;
    vistime = 33.0;
    nsnap = 2;
    well_depth = 100000;
    minx = 0;
    miny = 0;

    // Read parameters from stdin.
    readText pars(std::cin);

    for (size_t t(0); t < pars.getSize(); t++) {
        std::string line(pars[t]);
        readText::get(line, "instrdir", instrdir);
        readText::get(line, "flatdir", flatdir);
        readText::get(line, "tarfile", tarfile);
        readText::get(line, "filter", filter);
        readText::get(line, "nonlinear", nonlinear);
        readText::get(line, "welldepth", well_depth);
        readText::get(line, "parallelcte", parallelcte);
        readText::get(line, "serialcte", serialcte);
        readText::get(line, "readorder", readorder);
        readText::get(line, "adcerr", adcerr);
        readText::get(line, "obshistid", obshistid);
        readText::get(line, "exposureid", exposureid);
        readText::get(line, "chipid", chipid);
        readText::get(line, "obsseed", obsseed);
        readText::get(line, "vistime", vistime);
        readText::get(line, "nsnap", nsnap);
        readText::get(line, "parallelread", parallelread);
        readText::get(line, "serialread", serialread);
        readText::get(line, "parallelprescan", parallelPrescan);
        readText::get(line, "serialoverscan", serialOverscan);
        readText::get(line, "serialprescan", serialPrescan);
        readText::get(line, "paralleloverscan", parallelOverscan);
        readText::get(line, "bias", bias);
        readText::get(line, "gain", gain);
        readText::get(line, "readnoise", readnoise);
        readText::get(line, "darkcurrent", darkcurrent);
        readText::get(line, "hotpixelrate", hotpixelrate);
        readText::get(line, "hotcolumnrate", hotcolumnrate);
    }

    instr = "";
    unsigned pos = instrdir.rfind("/") + 1;
    for (unsigned i = pos; i<instrdir.length(); i++) instr += instrdir[i];

    if (flatdir == 1) instrdir = ".";

    if (adcerr.size() == 0) adcerr.resize(16, 0);

    std::istringstream focalplanePars(readText::get(instrdir + "/focalplanelayout.txt", chipid));
    double centerx_t, centery_t, pixsize_t;
    long pixelsx_t, pixelsy_t;
    focalplanePars >> centerx_t>> centery_t >> pixsize_t >> pixelsx_t >> pixelsy_t >> devtype >> devvalue;

    focalplanefile = instrdir + "/segmentation.txt";
    std::vector<std::string> amplifiers;
    readText::readSegmentation(focalplanefile, chipid, amplifiers);
    namp = (long) amplifiers.size();
    outchipid.resize(namp);
    outminx.resize(namp);
    outmaxx.resize(namp);
    outminy.resize(namp);
    outmaxy.resize(namp);
    crosstalk.resize(namp);
    for (long j = 0; j < namp; j++) {
        std::istringstream segmentationPars(amplifiers[j]);
        segmentationPars >> outchipid[j] >> outminx[j] >> outmaxx[j] >> outminy[j] >> outmaxy[j];
        double v;
        for (int i = 0; i < 16; i++) segmentationPars >> v;
        crosstalk[j].resize(namp);
        for (long k = 0; k < namp; k++)  segmentationPars >> crosstalk[j][k];
    }

    long seed = obsseed+exposureid;
    if (obsseed ==- 1)
        RngSetSeedFromTime();
    else
        RngSetSeed64(seed);

    RngUnwind(10000);


    if (devtype == "CMOS") {
        exptime = vistime/nsnap;
    } else {
        exptime = (vistime - (nsnap -1 )*devvalue)/nsnap;
    }

    infile << instr << "_e_"  << obshistid << "_f"<< filter << "_" << chipid << "_E" << std::setfill('0') << std::setw(3) << exposureid << ".fits.gz";
    read_fits_image(emap, onaxes, infile.str());
    adcmap.reserve(onaxes[0]*onaxes[1]);

    uint32 seedchip = 0;
    for(size_t m(0); m < chipid.size(); m++) seedchip += (uint32)(((int)chipid.c_str()[m]%10)*pow(10,m));
    RngSetSeed32_reseed(seedchip);

}

void E2adc::setHotpixels() {

    hotpixelmap.resize(onaxes[0]*onaxes[1],0);

    for (long l = 0;l < namp; l++) {
        for (long j = (outminy[l] - miny); j <= (outmaxy[l] - miny); j++) {
            for (long i = (outminx[l] - minx); i<= (outmaxx[l] - minx); i++) {
                if (RngDouble_reseed() < hotpixelrate[l]) {
                    hotpixelmap[onaxes[0]*j + i] = 1;
                }
            }
        }
    }
    for (long l = 0; l < namp; l++) {
        for (long j = (outminy[l] - miny); j <= (outmaxy[l] - miny); j++) {
            if (parallelread[l]==-1) {
                for (long i = (outminx[l] - minx); i <= (outmaxx[l] - minx); i++) {
                    if (RngDouble_reseed() < hotcolumnrate[l]/(outmaxx[l] - outminx[l] + 1)*2) {
                        for (long ii = i;ii <= (outmaxx[l] - minx); ii++) {
                            hotpixelmap[onaxes[0]*j + ii] = 1;
                        }
                        break;
                    }
                }
            } else {
                for (long i = (outmaxx[l] - minx);i >= (outminx[l] - minx); i--) {
                    if (RngDouble_reseed() < hotcolumnrate[l]/(outmaxx[l] - outminx[l] + 1)*2) {
                        for (long ii = i;ii >= (outminx[l] - minx); ii--) {
                            hotpixelmap[onaxes[0]*j + ii] = 1;
                        }
                        break;
                    }
                }
            }
        }
    }
}

void E2adc::convertADC() {

    std::vector<std::vector<float> > readoutmap;
    readoutmap.resize(namp);
    nx.resize(namp);
    ny.resize(namp);

    for (long l = 0; l < namp; l++) {
        std::cout << "Reading out chip " << chipid << " with amplifier chain " << outchipid[l] << std::endl;
        double nbackground = darkcurrent[l]*exptime;
        for (long i = (outminx[l] - minx); i <= (outmaxx[l] - minx); i++) {
            for (long j = (outminy[l] - miny); j <= (outmaxy[l] - miny); j++) {
                /* readnoise and dark current */
                long factor = (long)(random_gaussian()*sqrt(nbackground) + nbackground);
                if (factor < 0) factor = 0;
                adcmap[onaxes[0]*j + i] = emap[onaxes[0]*j + i] + factor;
                if (adcmap[onaxes[0]*j + i] > well_depth) adcmap[onaxes[0]*j + i] = well_depth;
            }
        }

        for (long i = (outminx[l] - minx); i <= (outmaxx[l] - minx); i++) {
            for (long j = (outminy[l] - miny); j <= (outmaxy[l] - miny); j++) {
                if (hotpixelmap[onaxes[0]*j + i] == 1) {
                    adcmap[onaxes[0]*j + i] = well_depth;
                }
            }
        }

        // small adc map
        long naxes[2];
        naxes[0] = outmaxx[l] - outminx[l] + 1;
        naxes[1] = outmaxy[l] - outminy[l] + 1;
        std::vector<float> smalladcmap(naxes[0]*naxes[1]);
        for (long i = 0; i < naxes[0]; i++) {
            for (long j = 0; j < naxes[1]; j++) {
                smalladcmap[naxes[0]*j + i] = adcmap[onaxes[0]*(j + (outminy[l] - miny))+(i + (outminx[l] - minx))];
            }
        }

        // Charge transfer inefficiency   a*b+(1-a)*b+(1-b)=a*b+b-a*b+1-b=1
        if (parallelread[l] == 1 && serialread[l] == -1) {
            for (long i = naxes[0] - 1;i >= 0; i--) {
                for (long j = 0; j < naxes[1]; j++) {
                    float origct = smalladcmap[naxes[0]*j + i];
                    for (long k = 0; k < ((long)origct); k++) {
                        if (RngDouble() < parallelcte) {
                            if (RngDouble() >= serialcte) {
                                if (j != naxes[1] - 1) {
                                    smalladcmap[naxes[0]*(j + 1) + i] += 1.0;
                                    smalladcmap[naxes[0]*j + i]-=1.0;
                                }
                            }
                        } else {
                            if (i!=0) {
                                smalladcmap[naxes[0]*j + (i - 1)] += 1.0;
                                smalladcmap[naxes[0]*j + i]-=1.0;
                            }
                        }
                    }
                }
            }
        }

        if (parallelread[l] == 1 && serialread[l] == 1) {
            for (long i = naxes[0] - 1; i >= 0; i--) {
                for (long j = naxes[1] - 1; j >= 0; j--) {
                    float origct = smalladcmap[naxes[0]*j + i];
                    for (long k = 0; k < ((long)origct); k++) {
                        if (RngDouble() < parallelcte) {
                            if (RngDouble() >= serialcte) {
                                if (j != 0) {
                                    smalladcmap[naxes[0]*(j - 1) + i] += 1.0;
                                    smalladcmap[naxes[0]*j + i]-=1.0;
                                }
                            }
                        } else {
                            if (i != 0) {
                                smalladcmap[naxes[0]*j + (i - 1)] += 1.0;
                                smalladcmap[naxes[0]*j + i] -= 1.0;
                            }
                        }
                    }
                }
            }
        }

        if (parallelread[l] == -1 && serialread[l] == -1) {
            for (long i = 0; i < naxes[0]; i++) {
                for (long j = 0; j < naxes[1]; j++) {
                    float origct = smalladcmap[naxes[0]*j + i];
                    for (long k = 0; k < ((long)origct); k++) {
                        if (RngDouble() < parallelcte) {
                            if (RngDouble() >= serialcte) {
                                if (j != naxes[1] - 1) {
                                    smalladcmap[naxes[0]*(j + 1) + i] += 1.0;
                                    smalladcmap[naxes[0]*j + i]-=1.0;
                                }
                            }
                        } else {
                            if (i != naxes[0] - 1) {
                                smalladcmap[naxes[0]*j + (i + 1)] += 1.0;
                                smalladcmap[naxes[0]*j + i]-=1.0;
                            }
                        }
                    }
                }
            }
        }

        if (parallelread[l] == -1 && serialread[l] == 1) {
            for (long i = 0;i < naxes[0]; i++) {
                for (long j = naxes[1] - 1; j >= 0; j--) {
                    float origct = smalladcmap[naxes[0]*j + i];
                    for (long k = 0; k < ((long)origct); k++) {
                        if (RngDouble() < parallelcte) {
                            if (RngDouble() >= serialcte) {
                                if (j != 0) {
                                    smalladcmap[naxes[0]*(j - 1) + i] += 1.0;
                                    smalladcmap[naxes[0]*j + i]-=1.0;
                                }
                            }
                        } else {
                            if (i != naxes[0] - 1) {
                                smalladcmap[naxes[0]*j + (i + 1)] += 1.0;
                                smalladcmap[naxes[0]*j + i]-=1.0;
                            }
                        }
                    }
                }
            }
        }

        // electron to ADC conversion
        for (long i = 0; i < naxes[0]; i++) {
            for (long j = 0; j < naxes[1]; j++) {
                smalladcmap[naxes[0]*j + i] = smalladcmap[naxes[0]*j + i]/(gain[l])/
                    (1 + nonlinear*(smalladcmap[naxes[0]*j + i]/well_depth)) + bias[l] + (long)(random_gaussian()*readnoise[l]);
            }
        }

        nx[l] = naxes[0] + serialOverscan[l] + serialPrescan[l];
        ny[l] = naxes[1] + parallelPrescan[l] + parallelOverscan[l];
        std::vector<float> fulladcmap(nx[l]*ny[l]);

        if (parallelread[l] == 1 && serialread[l] == -1) {
            for (long i = 0; i < naxes[0]; i++)
                for (long j = 0;j < naxes[1]; j++)
                    fulladcmap[nx[l]*(j + parallelPrescan[l]) + i + serialPrescan[l]] = smalladcmap[naxes[0]*j + i];
            for (long i = (nx[l] - serialOverscan[l]); i < nx[l]; i++)
                for (long j = 0; j < ny[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = 0; i < serialPrescan[l]; i++)
                for (long j = 0;j<ny[l];j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = 0; i < nx[l]; i++)
                for (long j = 0; j < parallelPrescan[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = 0; i < nx[l]; i++)
                for (long j = (ny[l] - parallelOverscan[l]); j < ny[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
        } else if (parallelread[l] == -1 && serialread[l] == -1) {
            for (long i = 0; i < naxes[0]; i++)
                for (long j = 0; j < naxes[1]; j++)
                    fulladcmap[nx[l]*(j + parallelPrescan[l]) + i + serialOverscan[l]] = smalladcmap[naxes[0]*j + i];
            for (long i = 0; i < serialOverscan[l]; i++)
                for (long j = 0; j < ny[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = (nx[l] - serialPrescan[l]); i < nx[l]; i++)
                for (long j = 0; j < ny[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = 0; i < nx[l]; i++)
                for (long j = 0; j < parallelPrescan[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = 0; i < nx[l]; i++)
                for (long j = (ny[l] - parallelOverscan[l]); j < ny[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
        } else if (parallelread[l] == 1 && serialread[l] == 1) {
            for (long i = 0;i < naxes[0]; i++)
                for (long j = 0; j < naxes[1]; j++)
                    fulladcmap[nx[l]*(j + parallelOverscan[l]) + i + serialPrescan[l]] = smalladcmap[naxes[0]*j + i];
            for (long i = (nx[l] - serialOverscan[l]);i < nx[l]; i++)
                for (long j = 0;j < ny[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = 0;i < serialPrescan[l]; i++)
                for (long j = 0; j < ny[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = 0; i < nx[l]; i++)
                for (long j = 0; j < parallelOverscan[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = 0; i < nx[l]; i++)
                for (long j = (ny[l] - parallelPrescan[l]); j < ny[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
        } else if (parallelread[l] == -1 && serialread[l] == 1) {
            for (long i = 0; i < naxes[0]; i++)
                for (long j = 0; j < naxes[1]; j++)
                    fulladcmap[nx[l]*(j + parallelOverscan[l]) + i + serialOverscan[l]] = smalladcmap[naxes[0]*j + i];
            for (long i = 0; i < serialOverscan[l]; i++)
                for (long j = 0; j < ny[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = (nx[l] - serialPrescan[l]); i < nx[l]; i++)
                for (long j = 0; j < ny[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = 0; i < nx[l]; i++)
                for (long j = 0; j < parallelOverscan[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
            for (long i = 0; i < nx[l]; i++)
                for (long j = (ny[l] - parallelPrescan[l]); j < ny[l]; j++)
                    fulladcmap[nx[l]*j + i] = bias[l] + (long)(random_gaussian()*readnoise[l]);
        }

        // ADC digitization
        std::vector<float> readoutmap_orig(nx[l]*ny[l]);
        for (long i = 0; i < nx[l]; i++) {
            for (long j = 0; j < ny[l]; j++) {
                readoutmap_orig[nx[l]*j + i] = 0.0;
                for (int k = 0; k < 16; k++) {
                    long bit = (((int)(fulladcmap[nx[l]*j + i]/(pow(2,k)+adcerr[k]))) % 2);
                    readoutmap_orig[nx[l]*j + i] +=  bit*pow(2,k);
                }
            }
        }

        // change to readout order
        readoutmap[l].resize(nx[l]*ny[l]);
        if (readorder == 1) {
            if (parallelread[l] == 1 && serialread[l] == -1) {
                for (long i = 0; i < nx[l]; i++)
                    for (long j = 0;j < ny[l]; j++)
                        readoutmap[l][ny[l]*i + j] = readoutmap_orig[nx[l]*j + (nx[l] - 1 - i)];
            } else if (parallelread[l] == -1 && serialread[l] == -1) {
                for (long i = 0; i < nx[l]; i++)
                    for (long j = 0; j < ny[l]; j++)
                        readoutmap[l][ny[l]*i + j] = readoutmap_orig[nx[l]*j + i];
            } else if (parallelread[l] == 1 && serialread[l] == 1) {
                for (long i = 0; i < nx[l]; i++)
                    for (long j = 0;j < ny[l]; j++)
                        readoutmap[l][ny[l]*i + j] = readoutmap_orig[nx[l]*(ny[l] - 1 - j)+(nx[l] - 1 - i)];
            } else if (parallelread[l] == -1 && serialread[l] == 1) {
                for (long i = 0; i < nx[l]; i++)
                    for (long j = 0; j < ny[l]; j++)
                        readoutmap[l][ny[l]*i + j] = readoutmap_orig[nx[l]*(ny[l] - 1 - j) + i];
            }
        } else {
            for (long i = 0; i < nx[l]*ny[l]; i++)  readoutmap[l][i] = readoutmap_orig[i];
        }
    }
    hotpixelmap.clear();
    adcmap.clear();
    emap.clear();

    // crosstalk
    fullReadoutMap.resize(namp);
    for (long l = 0; l < namp; l++) {
        size_t mapSize = readoutmap[l].size();
        fullReadoutMap[l].resize(mapSize);
        for (size_t i = 0; i < mapSize; i++) {
            float sum = 0.0;
            for (long k = 0; k < namp; k++)
                if (i<readoutmap[k].size()) sum += crosstalk[l][k]*readoutmap[k][i];
            fullReadoutMap[l][i] = (unsigned short) sum;
        }
    }
}

void E2adc::writeFitsImage() {

    std::string tarFiles(infile.str());
    for (long l = 0; l < namp; l++) {
        keyword_map keywords;
        keywords.insert(keyword_map::value_type("BZERO", keyProperties("TLONG", "32768", "offset data range to that of unsigned short")));
        keywords.insert(keyword_map::value_type("BSCALE", keyProperties("TLONG", "1", "default scaling factor")));
        keywords.insert(keyword_map::value_type("E2AIFILE", keyProperties("TSTRING", "", "E2adc input filename")));
        keywords.insert(keyword_map::value_type("E2AOFILE", keyProperties("TSTRING", "", "E2adc output filename")));
        keywords.insert(keyword_map::value_type("BIAS", keyProperties("TDOUBLE", bias[l], "Bias")));
        keywords.insert(keyword_map::value_type("GAIN", keyProperties("TDOUBLE", gain[l], "Gain")));
        keywords.insert(keyword_map::value_type("SCTE", keyProperties("TDOUBLE", serialcte, "Serial CTE")));
        keywords.insert(keyword_map::value_type("PCTE", keyProperties("TDOUBLE", parallelcte, "Parallel CTE")));
        keywords.insert(keyword_map::value_type("NONLIN", keyProperties("TDOUBLE", nonlinear, "Non-linear gain")));
        keywords.insert(keyword_map::value_type("E2AWLDP", keyProperties("TLONG", well_depth, "E2adc well depth")));
        keywords.insert(keyword_map::value_type("SATURATE", keyProperties("TFLOAT", (float)(well_depth/gain[l]+bias[l])*0.95, "Conservative saturation estimate")));
        keywords.insert(keyword_map::value_type("PREAD", keyProperties("TINT", parallelread[l], "Parallel read out direction")));
        keywords.insert(keyword_map::value_type("SREAD", keyProperties("TINT", serialread[l], "Serial read out direction")));
        keywords.insert(keyword_map::value_type("PSCANP", keyProperties("TINT", parallelPrescan[l], "Pre-scan parallel")));
        keywords.insert(keyword_map::value_type("OSCANS", keyProperties("TINT", serialOverscan[l], "Over-scan serial")));
        keywords.insert(keyword_map::value_type("PSCANS", keyProperties("TINT", serialPrescan[l], "Pre-scan serial")));
        keywords.insert(keyword_map::value_type("OSCANP", keyProperties("TINT", parallelOverscan[l], "Over-scan parallel")));
        for (int i = 0; i < 16; i++ ) {
            std::ostringstream ss;
            ss << "ADCER" << i;
            keywords.insert(keyword_map::value_type(ss.str(), keyProperties("TDOUBLE", adcerr[i], "ADC error bit")));
        }
        keywords.insert(keyword_map::value_type("E2AFPFL", keyProperties("TSTRING", focalplanefile, "E2adc focalplanefile")));
        keywords.insert(keyword_map::value_type("E2AICHI", keyProperties("TSTRING", chipid, "E2adc input chip ID")));
        keywords.insert(keyword_map::value_type("E2AOCHI", keyProperties("TSTRING", outchipid[l], "E2adc output chip ID")));
        unsigned pos = outchipid[l].find_last_of("_");
        keywords.insert(keyword_map::value_type("CCDID", keyProperties("TSTRING", outchipid[l].substr(0,pos), "CCD ID")));
        keywords.insert(keyword_map::value_type("AMPID", keyProperties("TSTRING", outchipid[l].substr(pos+1), "Amplifier ID")));
        keywords.insert(keyword_map::value_type("RDORDER", keyProperties("TINT", readorder, "0=CCS; 1=readorder")));
        keywords.insert(keyword_map::value_type("RDNOISE", keyProperties("TDOUBLE", (double)readnoise[l]/gain[l], "Readout noise (ADU/pixel)")));
        keywords.insert(keyword_map::value_type("DRKCURR", keyProperties("TDOUBLE", darkcurrent[l], "Dark Current (e-/pixel/s)")));

        if (readorder == 1) {
            std::ostringstream ss;
            ss.str("");
            ss<<"["<<std::setw(4)<<1<<":"<<std::setw(4)<<ny[l]
              <<","<<std::setw(4)<<1<<":"<<std::setw(4)<<serialOverscan[l]<<"]";
            keywords.insert(keyword_map::value_type("BIASSEC", keyProperties("TSTRING", ss.str(), "Scan section of amplifier")));
            ss.str("");
            ss<<"["<<std::setw(4)<<parallelPrescan[l]+1<<":"<<std::setw(4)<<ny[l]-parallelOverscan[l]
              <<","<<std::setw(4)<<serialOverscan[l]+1<<":"<<std::setw(4)<<nx[l]-serialPrescan[l]<<"]";
            keywords.insert(keyword_map::value_type("TRIMSEC", keyProperties("TSTRING", ss.str(), "Trimmed section of amplifier")));
            keywords.insert(keyword_map::value_type("DATASEC", keyProperties("TSTRING", ss.str(), "Data section of amplifier")));
        }

        std::ostringstream outfile;
        outfile << instr << "_a_" << obshistid << "_f"<< filter << "_" << outchipid[l] << "_E" << std::setfill('0') << std::setw(3) << exposureid << ".fits.gz";
        int fflag[] = {1, serialread[l], parallelread[l]};
        if (readorder == 0) fflag[0] = 0;
        float dcrpix[] = {-outminx[l], -outminy[l]};
        write_fits_image_cphead(fullReadoutMap[l], ny[l], nx[l], "!"+outfile.str(), infile.str(), dcrpix, fflag, keywords);
        tarFiles += " " + outfile.str();
    }
    if ( tarfile == 1 ) {
        std::ostringstream tarName;
        tarName << instr << "_ " << obshistid << "_f"<< filter << "_" << chipid << "_E" << std::setfill('0') << std::setw(3) << exposureid << ".tar";
        std::cout<<"Tarring "<<tarName.str()<<std::endl;
        std::string tarCommand = "tar cf " + tarName.str() + " " + tarFiles + " --remove-files";
        system(tarCommand.c_str());
    }
}


int main(void) {

    E2adc e2adc;

    e2adc.setup();
    e2adc.setHotpixels();
    e2adc.convertADC();
    e2adc.writeFitsImage();

    return 0;

}
