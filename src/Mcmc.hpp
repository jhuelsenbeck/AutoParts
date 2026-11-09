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
        void            run(int numChains);
        void            burnin(void);
        void            burnin(int numChains);
    
    private:
        void            closeOutputFiles(void);
        int             findColdChain(std::vector<int>& indices);
        std::string     formattedTime(Timer& t1, Timer& t2);
        void            openOutputFiles(void);
        double          power(int idx) { return powers[idx]; }
        void            printToScreen(int coldChainIdx, int n, double lnL, Timer& timePt, Timer& start);
        void            sample(int coldChainIdx, int n, double lnL);
        void            tuneProposals(int n);
        void            attemptSwap(std::vector<double>& lnL);
        void            printSwapSummary(void);
        void            tuneTemperature(void);
        std::ofstream   parmStrm;
        std::ofstream   treeStrm;
        std::vector<int> chainIndex;
        std::vector<UpdateInfo> tuning;     // proposal tuning for each chain temperature (0 = cold)
        std::vector<double> powers;         // MC^3 heating: the likelihood of the chain at temperature i is raised to powers[i]
        std::vector<double> logGaps;        // log(log(powers[i]) - log(powers[i+1])), what tuning adjusts
        bool            isTemperatureTuned;
        int             numTemperatureTunings;
        std::vector<int> swapTries;         // swaps tried and accepted between temperatures i and i+1
        std::vector<int> swapAccepts;
        std::vector<double> swapProbSum;    // for tuning: sum of the swap acceptance probabilities of each pair of
        int             swapProbCount;      // neighbouring temperatures, computed every cycle, and the number of cycles
        int             numCycles;
        int             burnInCycles;
        int             printFrequency;
        int             sampleFrequency;
        int             tuningFrequency;
        Model**         model;
        UserSettings*   settings;
};

#endif
