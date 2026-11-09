#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "Mcmc.hpp"
#include "Chunk.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"



Mcmc::Mcmc(Model** m, UserSettings* s) {

    model = m;
    settings = s;
    
    numCycles = settings->getNumMcmcCycles();
    burnInCycles = settings->getBurnIn();
    tuningFrequency = settings->getTuningFrequency();
    printFrequency = settings->getPrintFrequency();
    sampleFrequency = settings->getSampleFrequency();

    // Each chain temperature tunes its proposals separately, starting from the values
    // given by the user. The step sizes are scaled by the number of sites at a table
    // relative to the mean number of sites in a data subset (see UpdateInfo).
    int numChains = settings->getNumChains();
    UpdateInfo::setActive(nullptr);
    tuning.assign(numChains, UpdateInfo::updateInfo());
    double totalSites = 0.0;
    for (int i=0; i<model[0]->getNumSubsets(); i++)
        totalSites += model[0]->getChunk(i)->getNumObservedSites();
    UpdateInfo::setReferenceNumSites(totalSites / model[0]->getNumSubsets());
    UpdateInfo::setActive(&tuning[0]);

    // the powers start at 1/(1+i*T), T given by -temp; when they are tuned, each gap
    // between neighbouring powers (on the log scale) is adjusted separately
    double temperature = settings->getTemperature();
    powers.resize(numChains);
    for (int i=0; i<numChains; i++)
        powers[i] = 1.0 / (1.0 + i * temperature);
    logGaps.resize(std::max(numChains - 1, 0));
    for (int i=0; i<numChains-1; i++)
        logGaps[i] = std::log(std::log(powers[i]) - std::log(powers[i+1]));
    isTemperatureTuned = settings->getTuneTemperature();
    numTemperatureTunings = 0;
    swapTries.assign(std::max(numChains - 1, 0), 0);
    swapAccepts.assign(std::max(numChains - 1, 0), 0);
    swapProbSum.assign(std::max(numChains - 1, 0), 0.0);
    swapProbCount = 0;
}

void Mcmc::closeOutputFiles(void) {

    parmStrm.close();

    treeStrm << "end;" << std::endl; // Close the NEXUS block started at the start of sampling.
    treeStrm.close();
}

int Mcmc::findColdChain(std::vector<int>& indices) {

    for (int i=0; i<indices.size(); i++)
        {
        if (indices[i] == 0)
            return i;
        }
    Msg::error("Could not find cold chain index");
    return -1;
}

std::string Mcmc::formattedTime(Timer& t1, Timer& t2) {

    std::chrono::duration<double> durationSecs  = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
    int s = (int)durationSecs.count();
    int m = s / 60;
    int h = s / 3600;
        
    std::string tStr = "";
    if (h > 0)
        {
        tStr += std::to_string(h) + "h:";
        m -= h * 60;
        s -= h * 60 * 60;
        }
    if (m > 0 || (m == 0 && h > 0))
        {
        tStr += std::to_string(m) + "m:";
        s -= m * 60;
        }
    tStr += std::to_string(s) + "s";
    
    return tStr;
}

void Mcmc::openOutputFiles(void) {

    // open files for samples
    std::string outPath = settings->getOutputFile();
    std::string parmFileName = outPath + ".tsv";
    std::string treeFileName = outPath + ".tre";

    parmStrm.open( parmFileName.c_str(), std::ios::out );
    if (!parmStrm)
        Msg::error("Cannot open file \"" + parmFileName + "\"");
    treeStrm.open( treeFileName.c_str(), std::ios::out );
    if (!treeStrm)
        Msg::error("Cannot open file \"" + treeFileName + "\"");


    Tree* t = model[0]->getTree(0); // The tree shouldn't really matter for this
    treeStrm << "begin trees;" << std::endl; // Write the start of the NEXUS block that is closed when we close the stream.
    treeStrm << "   translate" << std::endl;
    std::vector<std::string> taxonNames = t->getTaxonNames();
    for (int i=0; i<taxonNames.size(); i++)
        {
        treeStrm << "   " << i+1 << " " << taxonNames[i];
        if (i == taxonNames.size()-1)
            treeStrm << ";" << std::endl;
        else
            treeStrm << "," << std::endl;
        }

}

void Mcmc::printToScreen(int coldChainIdx, int n, double lnL, Timer& t2, Timer& t1) {

    std::cout << "   * " << std::setw(6) << n << " -- " << std::fixed << std::setprecision(2) << lnL << " -- ";
    std::cout << "[" << coldChainIdx << "] ";
    std::cout << "tt(" << model[coldChainIdx]->getDegreeTree() << ") ";
    std::cout << "tl(" << model[coldChainIdx]->getDegreeTreeLength() << ") ";
    std::cout << "rv(" << model[coldChainIdx]->getDegreeGammaShape() << ") ";
    std::cout << "bf(" << model[coldChainIdx]->getDegreeBaseFrequencies() << ") ";
    std::cout << "er(" << model[coldChainIdx]->getDegreeRates() << ") -- ";

    std::chrono::duration<double> durationSecs  = std::chrono::duration_cast<std::chrono::seconds>(t2 - t1);
    double timePerCycle = (double)durationSecs.count() / n;
    if (timePerCycle == 0)
        timePerCycle = 1.0 / printFrequency;
    int s = (int)((numCycles - n) * timePerCycle);
    int m = s / 60;
    int h = s / 3600;

    if (h > 0)
        {
        std::cout << h << "h";
        std::cout << ":";
        m -= h * 60;
        s -= h * 60 * 60;
        }
    if (m > 0 || (m == 0 && h > 0))
        {
        if (m < 10)
            std::cout << "0";
        std::cout << m << "m";
        std::cout << ":";
        s -= m * 60;
        }
    if (s < 10)
        std::cout << "0";
    std::cout << s << "s";
    std::cout << " remaining  ";
    
    std::cout << std::endl;
}

// Single chain MCMC
void Mcmc::run(void) {

    std::cout << "   Markov Chain Monte Carlo Sampling:" << std::endl;
    // open files
    openOutputFiles();

    // run the chain
    auto start = std::chrono::high_resolution_clock::now();
    for (int n=1; n<=numCycles; n++)
        {
        // propose a new state for the chain
        UpdateInfo::setActive(&tuning[0]);
        double lnL = model[0]->update(1.0);
        
        if (n % printFrequency == 0)
            {
            auto timePt = std::chrono::high_resolution_clock::now();
            printToScreen(0, n, lnL, timePt, start);
            }
                    
        if (n == 1 || n % sampleFrequency == 0)
            sample(0, n, lnL);
        }
    std::cout << std::endl;
    tuning[0].summary();
    
    // close file
    closeOutputFiles();
}

// MCMCMC
void Mcmc::run(int numChains) {

    if (numChains == 1)
        {
        run();
        return;
        }
    
    // the swap statistics shown at the end are for the sampling only
    std::fill(swapTries.begin(), swapTries.end(), 0);
    std::fill(swapAccepts.begin(), swapAccepts.end(), 0);

    // This allows us to continue from burn-in, but also jump straight into sampling if we need to
    if(chainIndex.empty()){
        for (int i=0; i<numChains; i++)
            chainIndex.push_back(i);
    }
        
    std::cout << "   Markov Chain Monte Carlo Sampling:" << std::endl;
    // open files
    openOutputFiles();

    // run the chain
    auto start = std::chrono::high_resolution_clock::now();
    std::vector<double> lnL(numChains);
    for (int n=1; n<=numCycles; n++)
        {
        for (int k=0; k<numChains; k++)
            {
            // propose a new state for the chain
            UpdateInfo::setActive(&tuning[chainIndex[k]]);
            lnL[k] = model[k]->update( power(chainIndex[k]) );
            }
        
        attemptSwap(lnL);
            
        if (n % printFrequency == 0)
            {
            int coldChainOffset = findColdChain(chainIndex);
            auto timePt = std::chrono::high_resolution_clock::now();
            printToScreen(coldChainOffset, n, lnL[coldChainOffset], timePt, start);
            }
                    
        if (n == 1 || n % sampleFrequency == 0)
            {
            int coldChainOffset = findColdChain(chainIndex);
            sample(coldChainOffset, n, lnL[coldChainOffset]);
            }
        }
    std::cout << std::endl;
    tuning[0].summary();
    printSwapSummary();
    
    // close file
    closeOutputFiles();
}

void Mcmc::burnin(void) {

    std::cout << "   Markov Chain Monte Carlo Burn-In:" << std::endl;

    for (int n=1; n<=burnInCycles; n++)
        {
        // propose a new state for the chain
        UpdateInfo::setActive(&tuning[0]);
        double lnL = model[0]->update(1.0);
        
        if(n % printFrequency == 0 || n == 1)
            {
            std::cout << "   *    " << n << " -- " << lnL << std::endl;
            }

        if (n % tuningFrequency == 0)
            tuneProposals(n);
        }
    std::cout << std::endl;
    for (UpdateInfo& u : tuning)
        u.resetStats();
}

void Mcmc::burnin(int numChains) {

    if (numChains == 1)
        {
        burnin();
        return;
        }

    if(chainIndex.empty()){
        for (int i=0; i<numChains; i++)
            chainIndex.push_back(i);
    }
    
    std::cout << "   Markov Chain Monte Carlo Burn-In:" << std::endl;

    std::vector<double> lnL(numChains);
    for (int n=1; n<=burnInCycles; n++)
        {
        // propose a new state for the chain
        for (int k=0; k<numChains; k++)
            {
            UpdateInfo::setActive(&tuning[chainIndex[k]]);
            lnL[k] = model[k]->update( power(chainIndex[k]) );
            }
            
        attemptSwap(lnL);
        
        if(n % printFrequency == 0 || n == 1)
            {
            int coldChainOffset = findColdChain(chainIndex);
            std::cout << "   *    " << n << " -- " << lnL[coldChainOffset] << std::endl;
            }

        if (n % tuningFrequency == 0)
            tuneProposals(n);
        }
    std::cout << std::endl;
    for (UpdateInfo& u : tuning)
        u.resetStats();
    if (isTemperatureTuned == true)
        {
        std::cout << "   MC^3 powers after tuning =" << std::fixed << std::setprecision(4);
        for (double p : powers)
            std::cout << " " << p;
        std::cout << std::endl << std::endl;
        }
}

void Mcmc::sample(int coldChainIdx, int n, double lnL) {

    // sample the tree
    Tree* t = model[coldChainIdx]->getTree(0);
    treeStrm << "   tree t_" << n << " = " << t->getNewick() << std::endl;

    // sample the other parameters
    if (n == 1)
        {
        std::string parmHeader;
        model[coldChainIdx]->getHeader(parmHeader);
        parmStrm << "Gen" << '\t' << "lnL" << '\t' << parmHeader << std::endl;
        }
    std::string parmStr;
    model[coldChainIdx]->getParameterValues(parmStr);
    parmStrm << n << '\t' << lnL << '\t' << parmStr << std::endl;
}

void Mcmc::tuneProposals(int n) {

    // show the cold chain's acceptance rates for this batch now and again, then tune
    // the proposals of every chain temperature
    if (n % (10 * tuningFrequency) == 0 || n + tuningFrequency > burnInCycles)
        {
        std::cout << "\n   Tuning Parameters (cold chain, last " << tuningFrequency << " cycles)..." << std::endl;
        tuning[0].summary();
        printSwapSummary();
        }
    for (UpdateInfo& u : tuning)
        u.tune();
    tuneTemperature();
}

void Mcmc::attemptSwap(std::vector<double>& lnL) {

    // Propose to swap the states of the chains at two neighbouring temperatures, i and
    // i+1 (swaps between distant temperatures are seldom accepted). Only the likelihood
    // is heated (the target of a chain is L^power * prior), so the prior cancels in the
    // ratio.
    int numChains = (int)lnL.size();
    if (numChains < 2)
        return;
    RandomVariable& rng = RandomVariable::randomVariableInstance();

    // for tuning the powers, the probability that a swap between each pair of neighbouring
    // temperatures would be accepted (much less noisy than whether the one swap tried is)
    if (isTemperatureTuned == true)
        {
        std::vector<int> chainAt(numChains);
        for (int k=0; k<numChains; k++)
            chainAt[chainIndex[k]] = k;
        for (int t=0; t<numChains-1; t++)
            {
            double r = (lnL[chainAt[t]] - lnL[chainAt[t+1]]) * (power(t+1) - power(t));
            swapProbSum[t] += (r >= 0.0) ? 1.0 : std::exp(r);
            }
        swapProbCount++;
        }
    int t = (int)(rng.uniformRv() * (numChains - 1));
    int idx0 = -1, idx1 = -1;
    for (int k=0; k<numChains; k++)
        {
        if (chainIndex[k] == t)
            idx0 = k;
        else if (chainIndex[k] == t + 1)
            idx1 = k;
        }
    double pow0 = power(t);
    double pow1 = power(t + 1);
    double lnR = (lnL[idx0] - lnL[idx1]) * (pow1 - pow0);
    swapTries[t]++;
    if (std::log(rng.uniformRv()) < lnR)
        {
        std::swap(chainIndex[idx0], chainIndex[idx1]);
        swapAccepts[t]++;
        }
}

void Mcmc::printSwapSummary(void) {

    if (swapTries.empty() == true)
        return;
    std::cout << "   * MC^3 swaps between neighbouring temperatures:" << std::endl;
    for (int t=0; t<swapTries.size(); t++)
        {
        std::cout << "   *    power " << std::setprecision(3) << std::fixed << power(t) << " <-> " << power(t + 1) << ": ";
        std::cout << swapAccepts[t] << " / " << swapTries[t] << " accepted";
        if (swapTries[t] > 0)
            std::cout << " (" << std::setprecision(1) << 100.0 * swapAccepts[t] / swapTries[t] << "%)";
        std::cout << std::endl;
        }
    std::cout << std::endl;
}

void Mcmc::tuneTemperature(void) {

    // During the burn-in, adjust the heating so that the probability of accepting a swap
    // between each pair of neighbouring temperatures is about 30%: lower, and states rarely
    // pass between the chains; higher, and the heated chains are hardly heated. Each gap
    // between neighbouring powers is adjusted on its own (with a single temperature the
    // gap next to the cold chain is the hardest to cross). As for the proposals, this is
    // a Robbins-Monro step on the log scale, with steps that shrink; the powers are fixed
    // once the burn-in is over.
    const double target = 0.30;
    const double gain = 3.0;
    const double maxLogChange = std::log(2.0);
    if (isTemperatureTuned == false || swapProbSum.empty() == true || swapProbCount == 0)
        return;
    numTemperatureTunings++;
    double gamma = 1.0 / std::sqrt((double)numTemperatureTunings);
    for (int t=0; t<swapProbSum.size(); t++)
        {
        double acceptProb = swapProbSum[t] / swapProbCount;
        // accepting too often means the powers are too close: widen the gap
        double logChange = gain * gamma * (acceptProb - target);
        logChange = std::max(-maxLogChange, std::min(maxLogChange, logChange));
        logGaps[t] += logChange;
        logGaps[t] = std::max(std::log(1e-4), std::min(std::log(10.0), logGaps[t]));
        swapProbSum[t] = 0.0;
        }
    swapProbCount = 0;
    for (int t=0; t<logGaps.size(); t++)
        powers[t+1] = powers[t] * std::exp(-std::exp(logGaps[t]));
    std::fill(swapTries.begin(), swapTries.end(), 0);
    std::fill(swapAccepts.begin(), swapAccepts.end(), 0);
}
