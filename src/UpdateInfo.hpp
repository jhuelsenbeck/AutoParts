#ifndef UpdateInfo_hpp
#define UpdateInfo_hpp

#include "Msg.hpp"
#include <string>
#include <array>

struct UpdateStats {
    int numTries;
    int numAccepted;
};

enum UpdateType {
    SHAPE_SCALE = 0,
    BASE_FREQUENCY_BETA = 1,
    BASE_FREQUENCY_DIRICHLET = 2,
    RATE_BETA = 3,
    RATE_DIRICHLET = 4,
    TBR = 5,
    LOCAL = 6,
    BRANCH_PROPORTIONS = 7,
    TREE_LENGTH = 8,
    NO_MOVE = 9
};


inline const std::array<std::string, 9>& UpdateNames() {
    static std::array<std::string, 9> s = {
        "Gamma Shape",
        "Base Frequencies (Beta)",
        "Base Frequencies (Dirichlet)",
        "Exchangeability Rates (Beta)",
        "Exchangeability Rates (Dirichlet)",
        "Tree (TBR)",
        "Tree (LOCAL)",
        "Tree (Branch Proportions)",
        "Tree (Tree Length)"
    };
    return s;
}

// Acceptance statistics and tuning parameters for the proposals. There is one UpdateInfo
// for each chain temperature (heated chains need larger steps than the cold chain); the
// one in use is the "active" one, which Mcmc sets before it updates a chain, and which
// updateInfo() returns.
//
// The tuning parameter of a proposal on a table's value is scaled by the number of sites
// seated at that table (Restaurant sets it, with setNumSites(), before each proposal):
// the posterior of a subset's value narrows as the number of sites grows, so a fixed
// proposal would be too bold for large tables and too timid for small ones. For moves
// whose tuning parameter is a Dirichlet or Beta concentration (larger = smaller steps),
// the value used is c * (n / nRef); for moves whose tuning parameter is a window width,
// it is c * sqrt(nRef / n), where n is the number of sites at the table and nRef the mean
// number of sites in a data subset. The constants c are what tune() adjusts, and what the
// -t options set (they are the values used for a table of nRef sites).

class UpdateInfo {

    public:
                                            UpdateInfo(void);
        static UpdateInfo&                  updateInfo(void) { return *activeInfo(); }
        static void                         setActive(UpdateInfo* u) { activeInfo() = (u == nullptr ? defaultInfo() : u); }
        static void                         setReferenceNumSites(double n) { referenceNumSites() = n; }
        void                                accept(void);
        double                              attempt(UpdateType m);
        void                                reject(void);
        UpdateType                          getAttemptedUpdate(void) { return attemptedUpdate; }
        double                              getTunableParam(UpdateType m) { return tunableParams[m]; }
        void                                setAttemptedUpdate(UpdateType m) { attemptedUpdate = m; }
        void                                setNumSites(int n) { numSites = n; }
        void                                resetStats();
        void                                setTunableParam(UpdateType m, double v)
                                                {
                                                if(m < tunableParams.size()) tunableParams[m] = v;
                                                else Msg::error("");
                                                }
        void                                summary(void);
        void                                tune(void);
        
    private:
        static UpdateInfo*&                 activeInfo(void);
        static UpdateInfo*                  defaultInfo(void);
        static double&                      referenceNumSites(void);
        std::array<double, 9>               tunableParams;
        std::array<UpdateStats, 9>          stats;
        UpdateType                          attemptedUpdate = UpdateType::NO_MOVE;
        int                                 numSites = 0;       // sites at the table being updated (0 = do not scale)
        int                                 numTunings = 0;     // number of times tune() has adjusted the parameters
};

#endif
