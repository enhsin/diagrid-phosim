///
/// @package phosim
/// @file counter.h
/// @brief header for counter
///
/// @brief Created by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <time.h>
#include <string>
#include <vector>

struct Clog {
    long long rejected;
    long long removed;
    long long accepted;
    long long removed_dt;
    long long totalPhoton;
    double previousTime;
};

struct Tlog {
    double *throughput;
};

void counterInit (Clog *counterLog);
void counterCheck (Clog *counterLog, long sourcecounter, char *name);
void countBad (Clog *counterLog, long long photons, long long *ray);
void countBad_dt (Clog *counterLog, long long photons, long long *ray);
void countGood (Clog *counterLog, long long photons, long long *ray);

void writeThroughputFile (const std::string & outputdir, const std::string & outputfilename, Tlog *throughpuTlog, long nsurf);
void addThroughput (Tlog *throughpuTlog, long surf, long waveindex, long long sourceover);
void initThroughput (Tlog *throughpuTlog, long nsurf);
void writeCentroidFile (const std::string & outputdir, const std::string & outputfilename,
                        long long *source_saturation, long long *source_xpos, long long *source_ypos,
                        std::vector<double> source_id, long nsource);

// Waits for a new tick. Returns number of ticks since program start.
inline clock_t GetNewTick() {
    clock_t prev = clock();
    clock_t cur = clock();
    while(cur == prev)
        cur = clock();
    return cur;
}

// Converts ticks to seconds
inline double TicksToSec(clock_t ticks) {
    return (double)ticks/CLOCKS_PER_SEC;
}
