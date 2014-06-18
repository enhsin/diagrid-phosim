///
/// @package phosim
/// @file trim.cpp
/// @brief trim program: removes sources that have no chance of producing photons on a chip
///
/// @brief Created by:
/// @author Alan Meert (Purdue)
///
/// @brief Modified by:
/// @author Justin Bankert (Purdue)
/// @author John R. Peterson (Purdue)
/// @author En-Hsin Peng (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <stdio.h>
#include <zlib.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <string>
#include <iostream>
#include <fstream>

#include "trim/trim.h"
#include "raytrace/constants.h"
#include "ancillary/readtext.h"

using readtext::readText;

void Trim::xyPosition(double alpha, double delta, double *x, double *y ) {

    double a = cos(delta)*cos(alpha - pointingRA);
    double f = focalLength/(sin(pointingDec)*sin(delta) + a*cos(pointingDec));
    double yp = f*(cos(pointingDec)*sin(delta) - a*sin(pointingDec));
    double xp = f*cos(delta)*sin(alpha - pointingRA);
    *x = (xp*cos(rotationAngle) + yp*sin(rotationAngle));
    *y = (-xp*sin(rotationAngle) + yp*cos(rotationAngle));

}

void Trim::readCatalog() {

    gzFile gzfp[maxCatalog];
    FILE *fp[maxCatalog];
    FILE *fp2[maxChip];
    char line[4096];
    long sourceCount[maxChip];
    char ignore[2][4096];
    double ignore2;
    double ra, dec, mag;
    double x, y;
    double dx, dy, xp, yp;
    double currentBuffer;

    for (int d = 0; d < nChip; d++) {
        sprintf(outputFilename[d], "trimcatalog_%ld_%s.pars", obshistid, chipid[d].c_str());
    }

    for (int d = 0; d < nChip; d++) {
        fp2[d] = fopen(outputFilename[d], "wt");
        if (fp2[d]==NULL) {
            std::cout << "Cannot open for writing output file " << outputFilename[d] << "\n.";
            exit(1);
        }
        fprintf(fp2[d], "\n");
        fprintf(fp2[d], "\n");
        sourceCount[d] = 0;
    }

    for (int c = 0; c < nCatalog; c++) {
        if (strstr(catalog[c].c_str(), ".gz")==NULL) {

            fp[c] = fopen(catalog[c].c_str(), "rt");
            if (fp[c]==NULL){
                std::cout <<"Cannot open catalog file " << catalog[c] << ".\n";
                exit(1);
            }
            while (fgets( line, 4096, fp[c] )) {
                sscanf(line, "%s %lf %lf %lf %lf %s", ignore[0], &ignore2, &ra, &dec, &mag, ignore[1]);
                ra = ra*DEGREE;
                dec = dec*DEGREE;
                for (int d = 0; d < nChip; d++) {
                    currentBuffer = (buffer*pixelSize[d] + scale*pow(2.5, 17 - mag));
                    if (strayLight == 0 && currentBuffer > extendedBuffer + buffer) currentBuffer = extendedBuffer + buffer;
                    xyPosition(ra, dec, &x, &y);
                    dx = x - xPosition[d] - deltaX[d];
                    dy = y - yPosition[d] - deltaY[d];
                    xp = cos(angle[d])*dx + sin(angle[d])*dy;
                    yp = -sin(angle[d])*dx + cos(angle[d])*dy;
                    if ((fabs(xp) <= xDimension[d] + currentBuffer) &&
                        (fabs(yp) <= yDimension[d] + currentBuffer) ){
                        fprintf( fp2[d], "%s\n", line);
                        sourceCount[d]++;
                    }
                }
            }
            fclose(fp[c]);

        } else {

            gzfp[c] = gzopen( catalog[c].c_str(), "r" );
            if (gzfp[c] == NULL) {
                std::cout <<"Cannot open catalog file " << catalog[c] << ".\n";
                exit(1);
            }
            while (gzgets( gzfp[c] , line, 4096)) {
                sscanf(line, "%s %lf %lf %lf %lf %s", ignore[0], &ignore2, &ra, &dec, &mag, ignore[1]);
                ra = ra*DEGREE;
                dec = dec*DEGREE;
                for (int d = 0; d < nChip; d++) {
                    currentBuffer = (buffer*pixelSize[d] + scale*pow(2.5, 17 - mag));
                    if (strayLight == 0 && currentBuffer > extendedBuffer + buffer) currentBuffer = extendedBuffer + buffer;
                    xyPosition(ra, dec, &x, &y);
                    dx = x - xPosition[d] - deltaX[d];
                    dy = y - yPosition[d] - deltaY[d];
                    xp = cos(angle[d])*dx + sin(angle[d])*dy;
                    yp = -sin(angle[d])*dx + cos(angle[d])*dy;
                     if ((fabs(xp) <= xDimension[d] + currentBuffer) &&
                        (fabs(yp) <= yDimension[d] + currentBuffer)){
                        fprintf(fp2[d], "%s\n", line);
                        sourceCount[d]++;
                    }
                }
            }
            gzclose(gzfp[c]);

        }
    }

    for (int d = 0; d < nChip; d++) {
        fclose(fp2[d]);
    }

    long totalSourceCount = 0;
    for (int d = 0; d < nChip; d++) {
        totalSourceCount += sourceCount[d];
    }
    if (totalSourceCount>=minSource) {
        std::cout << "------------------------------------------------------------------------------------------\n";
        std::cout << "Trim Catalog\n";
        std::cout << "------------------------------------------------------------------------------------------\n";
        for (int d = 0; d < nChip; d++) {
            if (sourceCount[d] >= minSource) {
                std::cout << "Found " << sourceCount[d] << " source(s) for chip " << chipid[d] << ".\n";
            }
        }
    }

}

void Trim::getDetectorProperties( int d ){

    std::istringstream focalplanePars(readText::get(instrdir + "/focalplanelayout.txt", chipid[d]));
    focalplanePars >> xPosition[d] >> yPosition[d] >> pixelSize[d] >> xDimension[d] >> yDimension[d];
    xDimension[d] *= pixelSize[d]/2.0;
    yDimension[d] *= pixelSize[d]/2.0;
    std::string tempstring;
    double temp;
    focalplanePars >> tempstring;
    focalplanePars >> temp;
    focalplanePars >> tempstring;
    focalplanePars >> angle[d];
    angle[d] *= M_PI/180.;
    focalplanePars >> temp;
    focalplanePars >> temp;
    focalplanePars >> deltaX[d];
    focalplanePars >> deltaY[d];
    deltaX[d] *= 1000.0;
    deltaY[d] *= 1000.0;
}

void Trim::setup() {


    instrdir = "../data/lsst";
    flatDirectory = 0;

    readText pars(std::cin);

    for (size_t t(0); t<pars.getSize(); t++) {
        std::string line(pars[t]);
        readText::get(line, "pointingra", pointingRA);
        readText::get(line, "pointingdec", pointingDec);
        readText::get(line, "rotationangle", rotationAngle);
        readText::get(line, "filter", filter);
        readText::get(line, "obshistid", obshistid);
        readText::get(line, "instrdir", instrdir);
        readText::get(line, "flatdir", flatDirectory);
        readText::get(line, "minsource", minSource);
        readText::get(line, "chipid", chipid);
        readText::get(line, "catalog", catalog);
    }
    pointingRA *= DEGREE;
    pointingDec *= DEGREE;
    rotationAngle *= DEGREE;

    nCatalog = catalog.size();
    nChip = chipid.size();
    if (flatDirectory == 1) instrdir = ".";

    std::istringstream wavelengthPars(readText::get(instrdir + "/central_wavelengths.txt", filter));
    double wavelength;
    wavelengthPars >> wavelength >> plateScale;

    buffer = 100;
    strayLight = 1;
    focalLength = plateScale/DEGREE;
    scale = (0.2*ARCSEC/DEGREE)*(focalLength*DEGREE);
    extendedBuffer = (60.0*ARCSEC/DEGREE)*(focalLength*DEGREE);

    if (nChip > maxChip || nCatalog > maxCatalog ){
        std::cout << "Cannot split into that many files.\n";
        exit(1);
    }

}

int main() {

    Trim trim;

    trim.setup();
    for (int d = 0; d < trim.nChip; d++) {
        trim.getDetectorProperties(d);
    }
    trim.readCatalog();
    return(0);

}
