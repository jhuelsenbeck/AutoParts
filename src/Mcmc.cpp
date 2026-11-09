#include <cmath>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "Mcmc.hpp"
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

double Mcmc::power(int idx, double temperature) {

    return 1.0 / (1 + idx*temperature);
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
    UpdateInfo::updateInfo().summary();
    
    // close file
    closeOutputFiles();
}

// MCMCMC
void Mcmc::run(int numChains, double temperature) {

    if (numChains == 1)
        {
        run();
        return;
        }
    
    RandomVariable& rng = RandomVariable::randomVariableInstance();

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
            lnL[k] = model[k]->update( power(chainIndex[k], temperature) );
            }
        
        // attempt swap
        int idx0 = (int)(rng.uniformRv() * numChains);
        int idx1 = idx0;
        while (idx1 == idx0)
            idx1 = (int)(rng.uniformRv() * numChains);
        double pow0 = power(chainIndex[idx0], temperature);
        double pow1 = power(chainIndex[idx1], temperature);
        // only the likelihood is heated (target of a chain is L^power * prior), so the
        // prior cancels in the swap ratio
        double lnR = (lnL[idx0] - lnL[idx1]) * (pow1 - pow0);
        if (std::log(rng.uniformRv()) < lnR)
            {
            int tempIdx = chainIndex[idx0];
            chainIndex[idx0] = chainIndex[idx1];
            chainIndex[idx1] = tempIdx;
            }
            
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
    UpdateInfo::updateInfo().summary();
    
    // close file
    closeOutputFiles();
}

void Mcmc::burnin(void) {

    std::cout << "   Markov Chain Monte Carlo Burn-In:" << std::endl;

    for (int n=1; n<=burnInCycles; n++)
        {
        // propose a new state for the chain
        double lnL = model[0]->update(1.0);
        
        if(n % printFrequency == 0 || n == 1)
            {
            std::cout << "   *    " << n << " -- " << lnL << std::endl;
            }

        if (n % tuningFrequency == 0)
            {
            std::cout << "\n   Tuning Parameters..." << std::endl;
            UpdateInfo::updateInfo().summary();
            UpdateInfo::updateInfo().tune();
            }
        }
    std::cout << std::endl;
    UpdateInfo::updateInfo().resetStats();
}

void Mcmc::burnin(int numChains, double temperature) {

    if (numChains == 1)
        {
        burnin();
        return;
        }

    RandomVariable& rng = RandomVariable::randomVariableInstance();

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
            lnL[k] = model[k]->update( power(chainIndex[k], temperature) );
            }
            
        // swap
        int idx0 = (int)(rng.uniformRv() * numChains);
        int idx1 = idx0;
        while (idx1 == idx0)
            idx1 = (int)(rng.uniformRv() * numChains);
        double pow0 = power(chainIndex[idx0], temperature);
        double pow1 = power(chainIndex[idx1], temperature);
        // only the likelihood is heated (target of a chain is L^power * prior), so the
        // prior cancels in the swap ratio
        double lnR = (lnL[idx0] - lnL[idx1]) * (pow1 - pow0);
        if (std::log(rng.uniformRv()) < lnR)
            {
            int tempIdx = chainIndex[idx0];
            chainIndex[idx0] = chainIndex[idx1];
            chainIndex[idx1] = tempIdx;
            }
        
        if(n % printFrequency == 0 || n == 1)
            {
            int coldChainOffset = findColdChain(chainIndex);
            std::cout << "   *    " << n << " -- " << lnL[coldChainOffset] << std::endl;
            }

        if (n % tuningFrequency == 0)
            {
            std::cout << "\n   Tuning Parameters..." << std::endl;
            UpdateInfo::updateInfo().summary();
            UpdateInfo::updateInfo().tune();
            }
        }
    std::cout << std::endl;
    UpdateInfo::updateInfo().resetStats();
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
