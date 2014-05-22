///
/// @package phosim
/// @file rng_mwc.cpp
/// @brief random number generator functions
///
/// @brief Created by:
/// @author Kreso Cosic (Purdue)
///
/// @brief Modified by:
/// @author John R. Peterson (Purdue)
///
/// @warning This code is not fully validated
/// and not ready for full release.  Please
/// treat results with caution.
///

#include <assert.h>
#include <time.h>
#include <math.h>
#include "basic_types.h"

#include "rng_mwc.h"

namespace RandomNumbers {

    // returns a random number from
    // ======== normal distribution ==========
    // with mean=0 and variance=1
    double random_gaussian(void) {

        double u1 = 0.0, u2 = 0.0, v1 = 0.0, v2 = 0.0, s = 2.0;
        while(s >= 1) {
            u1 = RngDouble();
            u2 = RngDouble();
            v1 = 2.0*u1-1.0;
            v2 = 2.0*u2-1.0;
            s = pow(v1,2) + pow(v2,2);
        }
        double x1 = v1*sqrt((-2.0*log(s))/s);
        return x1;
    }

    // returns a random number from
    // ========= exponential distribution ==========
    // with rate=1
    double random_exp() {
        return -log(RngDouble());
    }


    // Poisson distributor (JRP)
    long long poisson(long long lambda) {
        long long count;
        double f,g;
        long long d;

        if (lambda<10) {
            f = RngDouble();
            d = 0;
            while (f >= exp(-(double)lambda)) {
                g = RngDouble();
                f = f*g;
                d++;
            }
            count = d;
        } else {
            count = (long long)((double)lambda + sqrt((double)lambda)*random_gaussian());
            if (count<0) count = 0;
        }
        return(count);
    }



    //double MWC rng state
    static uint32 m_z = 1234;
    static uint32 m_w = 42;

    //Marsaglia's double MWC random number generator
    uint32 RngUint32() {
        m_z = 36969 * (m_z & 65535) + (m_z >> 16);
        m_w = 18000 * (m_w & 65535) + (m_w >> 16);
        return ((m_z << 16) + m_w);  /* 32-bit result */
    }

    //derived
    // ========= uniform distribution ==========
    //returns number from interval <0,1>
    //conversion bias removed
    float RngFloat() {
        const uint32 rng_floatMask = 0x007FFFFF;
        const float intToFloat = float(1.0/(rng_floatMask + 1.0));
        const float plusBiasFloat = intToFloat/2;

        uint32 rn = RngUint32()&rng_floatMask;
        return rn*intToFloat + plusBiasFloat;
    }

    //derived
    // ========= uniform distribution ==========
    //returns number from interval <0,1>
    //conversion bias removed
    double RngDouble() {
        const double intToDouble = 1.0/(uint32_MAX + 1.0);
        const double plusBiasDouble = intToDouble/2;

        uint32 rn = RngUint32();
        return rn*intToDouble + plusBiasDouble;
    }

    // Draw 'count' numbers from rng
    void RngUnwind(int count) {
        for (int i = 0;i<count;i++)
            RngUint32();
    }

    time_t timeSec1970() {    return time(NULL); }

    //Note: values of m_z and m_w must not be zero (neither one)!!!!!
    void set_seed_rng32Mwc(uint32 z, uint32 w) {

        assert(w!=0 && z!=0);
        m_z = z;
        m_w = w;

    }

    // set seed from 64-bit value
    void RngSetSeed64(uint64 seed) {

        uint32 z = uint32(seed>>32);
        uint32 w = uint32(seed&0xFFFFFFFF);
        if (z==0) z = 4232482;
        if (w==0) w = 1234628;

        set_seed_rng32Mwc(z,w);
    }

    // set seed from 32-bit value
    // mixes up bits from the seed to produce a 64-bit seed
    void RngSetSeed32(uint32 seed) {

        uint32 lobits = seed;
        uint32 hibits = seed*91*53;
        uint64 seed64 = (uint64(hibits)<<32) + lobits;

        RngSetSeed64(seed64);
    }

    // set seed from current time and clock()
    void RngSetSeedFromTime() {

        uint32 lobits = uint32(clock()) + uint32(timeSec1970());
        uint32 hibits = uint32(clock())*53 + uint32(timeSec1970())*91;
        uint64 seed64 = ((uint64)hibits<<32) + lobits;

        RngSetSeed64(seed64);
    }

    // SECOND SET OF RANDOM NUMBER FUNCTIONS WHERE WE MESS WITH SEED


    // returns a random number from
    // ======== normal distribution ==========
    // with mean=0 and variance=1
    double random_gaussian_reseed(void) {

        double u1 = 0.0, u2 = 0.0, v1 = 0.0, v2 = 0.0, s = 2.0;
        while (s >= 1) {
            u1 = RngDouble_reseed();
            u2 = RngDouble_reseed();
            v1 = 2.0*u1 - 1.0;
            v2 = 2.0*u2 - 1.0;
            s = pow(v1,2) + pow(v2,2);
        }
        double x1 = v1*sqrt((-2.0*log(s))/s);
        return x1;

    }

    // returns a random number from
    // ========= exponential distribution ==========
    // with rate=1
    double random_exp_reseed() {
        return -log(RngDouble_reseed());
    }

    //double MWC rng state
    static uint32 m_z_reseed = 1234;
    static uint32 m_w_reseed = 42;

    //Marsaglia's double MWC random number generator
    uint32 RngUint32_reseed() {
        m_z_reseed  =  36969 * (m_z_reseed & 65535) + (m_z_reseed >> 16);
        m_w_reseed  =  18000 * (m_w_reseed & 65535) + (m_w_reseed >> 16);
        return ((m_z_reseed << 16) + m_w_reseed);  /* 32-bit result */
    }

    //derived
    // ========= uniform distribution ==========
    //returns number from interval <0,1>
    //conversion bias removed
    float RngFloat_reseed() {
        const uint32 rng_floatMask = 0x007FFFFF;
        const float intToFloat = float(1.0/(rng_floatMask + 1.0));
        const float plusBiasFloat = intToFloat/2;

        uint32 rn = RngUint32_reseed()&rng_floatMask;
        return rn*intToFloat + plusBiasFloat;
    }

    //derived
    // ========= uniform distribution ==========
    //returns number from interval <0,1>
    //conversion bias removed
    double RngDouble_reseed() {
        const double intToDouble = 1.0/(uint32_MAX + 1.0);
        const double plusBiasDouble = intToDouble/2;

        uint32 rn = RngUint32_reseed();
        return rn*intToDouble + plusBiasDouble;
    }

    // Draw 'count' numbers from rng
    void RngUnwind_reseed(int count) {

        for (int i = 0;i<count;i++)
            RngUint32_reseed();
    }

    //Note: values of m_z and m_w must not be zero (neither one)!!!!!
    void set_seed_rng32Mwc_reseed(uint32 z, uint32 w) {
        assert(w != 0 && z != 0);
        m_z_reseed = z;
        m_w_reseed = w;
    }

    // set seed from 64-bit value
    void RngSetSeed64_reseed(uint64 seed) {

        uint32 z = uint32(seed>>32);
        uint32 w = uint32(seed&0xFFFFFFFF);
        if (z==0) z=4232482;
        if (w==0) w = 1234628;
        set_seed_rng32Mwc_reseed(z, w);

    }

    // set seed from 32-bit value
    // mixes up bits from the seed to produce a 64-bit seed
    void RngSetSeed32_reseed(uint32 seed) {

        uint32 lobits = seed;
        uint32 hibits = seed*91*53;
        uint64 seed64 = (uint64(hibits)<<32) + lobits;

        RngSetSeed64_reseed(seed64);
    }

    // set seed from current time and clock()
    void RngSetSeedFromTime_reseed() {

        uint32 lobits = uint32(clock()) + uint32(timeSec1970());
        uint32 hibits = uint32(clock())*53 + uint32(timeSec1970())*91;
        uint64 seed64 = ((uint64)hibits<<32) + lobits;

        RngSetSeed64_reseed(seed64);
    }

    //return the current seed
    void RngGetSeed(uint32* z, uint32* w) {
        *z = m_z;
        *w = m_w;
    }



} // namespace RandomNumbers ends
