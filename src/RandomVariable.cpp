#include <cmath>
#include <iostream>
#include <limits>
#include "Msg.hpp"
#include "RandomVariable.hpp"



RandomVariable::RandomVariable(void) {

    uint32_t seed = (uint32_t)time(NULL);
    initialize(seed);
}

RandomVariable::RandomVariable(uint32_t seed) {

    initialize(seed);
}

uint32_t RandomVariable::extractU32(void) {

    int i = index;
    if (index >= N)
        {
        twist();
        i = index;
        }

    uint32_t y = mt[i];
    index = i + 1;

    y ^= (mt[i] >> U);
    y ^= (y << S) & B;
    y ^= (y << T) & C;
    y ^= (y >> L);

    return y;
}


void RandomVariable::initialize(uint32_t seed) {

    mt[0] = seed;
    for (uint32_t i=1; i<N; i++)
        {
        mt[i] = (F * (mt[i - 1] ^ (mt[i - 1] >> 30)) + i);
        }
    index = N;
    
    for (size_t i=0; i<10000; i++)
        extractU32();
}



void RandomVariable::twist(void) {

    for (uint32_t i=0; i<N; i++)
        {
        uint32_t x = (mt[i] & MASK_UPPER) + (mt[(i + 1) % N] & MASK_LOWER);
        uint32_t xA = x >> 1;

        if ( x & 0x1 )
            xA ^= A;

        mt[i] = mt[(i + M) % N] ^ xA;
        }
    index = 0;
}

double RandomVariable::uniformRv(void) {

    // strictly inside (0,1): dividing by UINT32_MAX gave exactly 1.0 (and 0.0) about once
    // in 2^32 draws, so that (int)(u * n) indexed past the end of an array and -log(u)
    // was infinite
    return ((double)extractU32() + 0.5) / 4294967296.0;
}
