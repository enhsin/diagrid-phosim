///
/// @package phosim
/// @file counter.cpp
/// @brief various logging functions
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <stdio.h>
#include <stdlib.h>

#include "counter.h"

void counterInit(Clog *counterLog) {


    fprintf(stdout,"------------------------------------------------------------------------------------------\n");
    fprintf(stdout,"Type                Sources         Photons  (Sat,Rem,Rej,Acc)%%  Time (s)       Photons/s\n");
    fprintf(stdout,"------------------------------------------------------------------------------------------\n");

    counterLog->rejected=0;
    counterLog->accepted=0;
    counterLog->removed=0;
    counterLog->removed_dt=0;
    counterLog->totalPhoton=0;
    counterLog->previousTime=TicksToSec(GetNewTick());

}

void counterCheck(Clog *counterLog, long sourcecounter, char *name) {

    double newTime;
    long rate;
    char sourceString[4096];
    char photonString[4096];
    char rateString[4096];

    newTime=TicksToSec(clock());
    rate=(long)(counterLog->totalPhoton/(newTime-counterLog->previousTime+1e-2));

    if (sourcecounter<1000) sprintf(sourceString,"%7ld",sourcecounter);
    if (sourcecounter>=1000) sprintf(sourceString,"%3ld,%03ld",sourcecounter/1000,sourcecounter%1000);

    if (counterLog->totalPhoton<1000) sprintf(photonString,"%15lld",counterLog->totalPhoton);
    if (counterLog->totalPhoton>=1000 && counterLog->totalPhoton<1000000)
        sprintf(photonString,"%11lld,%03lld",counterLog->totalPhoton/1000,counterLog->totalPhoton%1000);
    if (counterLog->totalPhoton>=1000000 && counterLog->totalPhoton<1000000000)
        sprintf(photonString,"%7lld,%03lld,%03lld",counterLog->totalPhoton/1000000,(counterLog->totalPhoton/1000)%1000,counterLog->totalPhoton%1000);
    if (counterLog->totalPhoton>=1000000000)
        sprintf(photonString,"%3lld,%03lld,%03lld,%03lld",(counterLog->totalPhoton/1000000000),(counterLog->totalPhoton/1000000)%1000,(counterLog->totalPhoton/1000)%1000,counterLog->totalPhoton%1000);

    if (rate<1000) sprintf(rateString,"%15ld",rate);
    if (rate>=1000 && rate<1000000) sprintf(rateString,"%11ld,%03ld",rate/1000,rate%1000);
    if (rate>=1000000 && rate<1000000000) sprintf(rateString,"%7ld,%03ld,%03ld",rate/1000000,(rate/1000)%1000,rate%1000);
    if (rate>=1000000000) sprintf(rateString,"%3ld,%03ld,%03ld,%03ld",(rate/1000000000),(rate/1000000)%1000,(rate/1000)%1000,rate%1000);

    fprintf(stdout,"%s %s %s  (%3.0f,%3.0f,%3.0f,%3.0f)  %9.2f %s\n",name,sourceString,photonString,(double)counterLog->rejected/(double)counterLog->totalPhoton*100,(double)counterLog->removed_dt/(double)counterLog->totalPhoton*100,(double)counterLog->removed/(double)counterLog->totalPhoton*100,(double)counterLog->accepted/(double)counterLog->totalPhoton*100,newTime-counterLog->previousTime,rateString);

    counterLog->previousTime=newTime;
    counterLog->rejected=0;
    counterLog->accepted=0;
    counterLog->removed=0;
    counterLog->removed_dt=0;
    counterLog->totalPhoton=0;
}


void countGood(Clog *counterLog, long long photons, long long *ray) {

    counterLog->accepted+=1;
    counterLog->rejected+=(photons-1);
    counterLog->totalPhoton+=photons;
    *ray+=photons;

}

void countBad(Clog *counterLog, long long photons, long long *ray) {

    counterLog->removed+=photons;
    counterLog->totalPhoton+=photons;
    *ray+=photons;

}
void countBad_dt(Clog *counterLog, long long photons, long long *ray) {

    counterLog->removed_dt+=photons;
    counterLog->totalPhoton+=photons;
    *ray+=photons;

}

void addThroughput (Tlog *throughputlog, long surf, long waveindex, long long sourceover) {

    throughputlog->throughput[(surf+1)*901+waveindex]+=(double)sourceover;
}

void initThroughput (Tlog *throughputlog, long nsurf) {

    throughputlog->throughput=(double*)calloc((nsurf+2)*901,sizeof(double));

}

void writeThroughputFile (const std::string & outputdir, const std::string & outputfilename, Tlog *throughputlog, long nsurf) {

    FILE *outdafile;
    long i,k;
    char tempstring[4096];

        sprintf(tempstring,"%s/throughput_%s.txt",outputdir.c_str(),outputfilename.c_str());
        outdafile=fopen(tempstring,"w");
        for (k=0;k<901;k++) {
            fprintf(outdafile,"%ld ",k);
            for (i=0;i<nsurf+2;i++) {
                fprintf(outdafile,"%lf ",throughputlog->throughput[i*901+k]);
            }
            fprintf(outdafile,"\n");
        }
        fclose(outdafile);


}

void writeCentroidFile (const std::string & outputdir, const std::string & outputfilename, long long *source_saturation, long long *source_xpos, long long *source_ypos, std::vector<double> source_id, long nsource) {

    FILE *outdafile;
    long k;
    char tempstring[4096];


        sprintf(tempstring,"%s/centroid_%s.txt",outputdir.c_str(),outputfilename.c_str());
        outdafile=fopen(tempstring,"w");
        fprintf(outdafile,"SourceID Photons AvgX AvgY\n");
        for (k=0;k<nsource;k++) {
            fprintf(outdafile,"%lf %lld %lf %lf\n",source_id[k],source_saturation[k],((double)source_xpos[k])/((double)source_saturation[k]),((double)source_ypos[k])/((double)source_saturation[k]));
        }
        fclose(outdafile);

}
