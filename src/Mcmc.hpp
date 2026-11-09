#ifndef Mcmc_hpp
#define Mcmc_hpp

#include <chrono>
#include <fstream>
#include <string>
#include <vector>
#include "UpdateInfo.hpp"
class Model;
class UserSettings;

typedef std::chrono::high_resolution_clock::time_point Timer;



class Mcmc {

    public:
                        Mcmc(Model** m, UserSettings* s);
                       ~Mcmc(void) { UpdateInfo::setActive(nullptr); }
        void            run(void);
        void            run(int numChains, double temperature);
        void            burnin(void);
        void            burnin(int numChains, double temperature);
    
    private:
        void            closeOutputFiles(void);
        int             findColdChain(std::vector<int>& indices);
        std::string     formattedTime(Timer& t1, Timer& t2);
        void            openOutputFiles(void);
        double          power(int idx, double temperature);
        void            printToScreen(int coldChainIdx, int n, double lnL, Timer& timePt, Timer& start);
        void            sample(int coldChainIdx, int n, double lnL);
        void            tuneProposals(int n);
        std::ofstream   parmStrm;
        std::ofstream   treeStrm;
        std::vector<int> chainIndex;
        std::vector<UpdateInfo> tuning;     // proposal tuning for each chain temperature (0 = cold)
        int             numCycles;
        int             burnInCycles;
        int             printFrequency;
        int             sampleFrequency;
        int             tuningFrequency;
        Model**         model;
        UserSettings*   settings;
};

#endif
