namespace RandomNumbers
{

    // ========= uniform distribution ==========
    //returns number from interval <0,1>
    float RngFloat();
    double RngDouble();

    // returns a random number from
    // ======== normal distribution ==========
    // with mean=0 and variance=1
    double random_gaussian(void);

    // returns a random number from
    // ========= exponential distribution ==========
    // with rate=1
    double random_exp();

    long long poisson(long long lambda);

    void set_seed_rng32Mwc(uint32 z, uint32 w);

    // set seed from 64-bit value
    void RngSetSeed64(uint64 seed);

    // set seed from 32-bit value
    // mixes up bits from the seed to produce a 64-bit seed
    void RngSetSeed32(uint32 seed);

    // set seed by using time and clock functions for entrpy
    void RngSetSeedFromTime();

    // Draw 'count' numbers from rng
    void RngUnwind(int count);

    // ========= uniform distribution ==========
    //returns number from interval <0,1>
    float RngFloat_reseed();
    double RngDouble_reseed();

    // returns a random number from
    // ======== normal distribution ==========
    // with mean=0 and variance=1
    double random_gaussian_reseed(void);

    // returns a random number from
    // ========= exponential distribution ==========
    // with rate=1
    double random_exp_reseed();

    // set seed from 64-bit value
    void RngSetSeed64_reseed(uint64 seed);

    // set seed from 32-bit value
    // mixes up bits from the seed to produce a 64-bit seed
    void RngSetSeed32_reseed(uint32 seed);

    // set seed by using time and clock functions for entrpy
    void RngSetSeedFromTime_reseed();

    // Draw 'count' numbers from rng
    void RngUnwind_reseed(int count);

    //return the current seed
    void RngGetSeed(uint32* z, uint32* w);

} // namespace RandomNumbers ends
