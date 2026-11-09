#ifndef RandomVariable_hpp
#define RandomVariable_hpp

#include <cstdint>
#include <vector>

enum {

    N = 624,
    M = 397,
    R = 31,
    A = 0x9908B0DF,
    F = 1812433253,
    U = 11,
    S = 7,
    B = 0x9D2C5680,
    T = 15,
    C = 0xEFC60000,
    L = 18,
    MASK_LOWER = (1ull << R) - 1,
    MASK_UPPER = (1ull << R)
};

class RandomVariable {

    public:
        static RandomVariable&  randomVariableInstance(void)
                                    {
                                    static RandomVariable singleRandomVariable;
                                    return singleRandomVariable;
                                    }
        void                    setSeed(uint32_t seed) { initialize(seed); }
        double                  uniformRv(void);

    private:
                                RandomVariable(void);
                                RandomVariable(RandomVariable& r);
                                RandomVariable(uint32_t seed);
        RandomVariable&         operator=(const RandomVariable&);
        uint32_t                extractU32(void);
        void                    initialize(uint32_t seed);
        void                    twist(void);
        uint16_t                index;
        uint32_t                mt[N];
};

#endif
