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

class UpdateInfo {

    public:
        static UpdateInfo&                  updateInfo(void)
                                                {
                                                static UpdateInfo singleUpdateInfo;
                                                return singleUpdateInfo;
                                                }
        void                                accept(void);
        double                              attempt(UpdateType m) 
                                                { 
                                                attemptedUpdate = m; 
                                                if(m < tunableParams.size()) return tunableParams[m];
                                                else return 1.0;
                                                }
        void                                reject(void);
        UpdateType                          getAttemptedUpdate(void) { return attemptedUpdate; }
        void                                setAttemptedUpdate(UpdateType m) { attemptedUpdate = m; }
        void                                resetStats();
        void                                setTunableParam(UpdateType m, double v)
                                                {
                                                if(m < tunableParams.size()) tunableParams[m] = v;
                                                else Msg::error("");
                                                }
        void                                summary(void);
        void                                tune(void);
        
    private:
            
                                            UpdateInfo(void) {}
                                            UpdateInfo(UpdateInfo& r);
        std::array<double, 9>               tunableParams;
        std::array<UpdateStats, 9>          stats;
        UpdateType                          attemptedUpdate = UpdateType::NO_MOVE;
};

#endif
