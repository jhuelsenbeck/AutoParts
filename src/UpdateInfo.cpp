#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include "UpdateInfo.hpp"



namespace {

    // whether a move's tuning parameter is a concentration (larger = smaller steps), as
    // opposed to a window width (larger = larger steps)
    bool isConcentration(int m) {

        return m == BASE_FREQUENCY_BETA || m == BASE_FREQUENCY_DIRICHLET || m == RATE_BETA ||
               m == RATE_DIRICHLET || m == BRANCH_PROPORTIONS;
    }

    // the acceptance rate tune() aims for: about 0.44 suits a one-dimensional move, and a
    // little less a move of several values at once (Roberts and Rosenthal 2001)
    double targetAcceptance(int m) {

        if (m == BASE_FREQUENCY_DIRICHLET || m == RATE_DIRICHLET)
            return 0.30;
        return 0.44;
    }
}

UpdateInfo::UpdateInfo(void) {

    tunableParams.fill(1.0);
    resetStats();
}

UpdateInfo* UpdateInfo::defaultInfo(void) {

    // holds the tuning parameters given by the user; active when no chain is running
    static UpdateInfo info;
    return &info;
}

UpdateInfo*& UpdateInfo::activeInfo(void) {

    static UpdateInfo* active = defaultInfo();
    return active;
}

double& UpdateInfo::referenceNumSites(void) {

    static double n = 0.0;
    return n;
}

double UpdateInfo::attempt(UpdateType m) {

    attemptedUpdate = m;
    if (m >= tunableParams.size())
        return 1.0;
    double c = tunableParams[m];
    double nRef = referenceNumSites();
    if (numSites <= 0 || nRef <= 0.0 || m == TBR || m == LOCAL)
        return c;
    if (isConcentration(m) == true)
        return c * (numSites / nRef);
    return c * std::sqrt(nRef / numSites);
}

void UpdateInfo::accept(void) {
    if (attemptedUpdate == UpdateType::NO_MOVE)
        Msg::error("Tried to accept a move when there was none!");

    auto& moveStats = stats[attemptedUpdate];
    moveStats.numAccepted++;
    moveStats.numTries++;
        
    attemptedUpdate = UpdateType::NO_MOVE;
}

void UpdateInfo::reject(void) {
    if (attemptedUpdate == UpdateType::NO_MOVE)
        Msg::error("Tried to accept a move when there was none!");

    auto& moveStats = stats[attemptedUpdate];
    moveStats.numTries++;
        
    attemptedUpdate = UpdateType::NO_MOVE;
}

void UpdateInfo::resetStats(void) {
    for(auto& s : stats){
        s.numAccepted = 0;
        s.numTries = 0;
    }
}

void UpdateInfo::summary(void) {

    const auto& updateNames = UpdateNames();

    int paramWidth = 0;
    for (const auto& name : updateNames)
        paramWidth = std::max(paramWidth, (int)name.size());

    int triesWidth = 8;
    int accWidth = 8;
    int freqWidth  = 8;

    std::cout << "   * MCMC Acceptance Rates:\n";
    std::cout << "   * "
              << std::left << std::setw(paramWidth) << "Parameter"
              << " "
              << std::right << std::setw(triesWidth) << "Tries"
              << std::setw(accWidth)                 << "Accep"
              << std::setw(freqWidth)                << "Freq"
              << std::setw(12)                       << "Tuning"
              << "\n";

    std::cout << "   * "
              << std::string(paramWidth, '-')
              << " "
              << std::string(triesWidth + accWidth + freqWidth + 12, '-')
              << "\n";

    for (int i = 0; i < updateNames.size(); i++)
        {
        const auto& name = updateNames[i];
        auto& s = stats[i];

        std::cout << "   * "
                  << std::left  << std::setw(paramWidth) << name
                  << " "
                  << std::right << std::setw(triesWidth) << s.numTries
                  << std::setw(accWidth)                 << s.numAccepted
                  << std::fixed << std::setprecision(1)
                  << std::setw(freqWidth-1)              // Minus 1 for % sign
                  << (s.numTries > 0 ?
                      (double)s.numAccepted / s.numTries * 100.0 : 0.0)
                  << "%"
                  << std::setprecision(4) << std::setw(12) << tunableParams[i]
                  << "\n";
        }

    std::cout << std::endl;
}


void UpdateInfo::tune(void) {

    // A stochastic approximation (Robbins-Monro) step on the log of each tuning
    // parameter, toward the target acceptance rate. The steps shrink as tuning goes on,
    // so the parameters settle rather than jump about on the noise of each batch. Tuning
    // only happens during the burn-in, so the sampling afterwards is a fixed, valid chain.
    // A move with too few tries in this batch keeps its counts for the next one.
    const int minTries = 20;
    const double gain = 3.0;
    const double maxLogChange = std::log(3.0);
    numTunings++;
    double gamma = 1.0 / std::sqrt((double)numTunings);
    for (int i=0; i<9; i++)
        {
        auto& s = stats[i];
        if (i == TBR || i == LOCAL)
            {
            // these are not tuned on their acceptance rates (only their counts are reset)
            s.numAccepted = 0;
            s.numTries = 0;
            continue;
            }
        if (s.numTries < minTries)
            continue;
        double acceptRate = (double)s.numAccepted / s.numTries;
        double logChange = gain * gamma * (acceptRate - targetAcceptance(i));
        logChange = std::max(-maxLogChange, std::min(maxLogChange, logChange));
        // accepting too often means the steps are too small: make them larger
        if (isConcentration(i) == true)
            tunableParams[i] *= std::exp(-logChange);
        else
            tunableParams[i] *= std::exp(logChange);
        s.numAccepted = 0;
        s.numTries = 0;
        }
}
