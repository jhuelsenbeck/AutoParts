#include <cmath>
#include <ctime>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>
#include "Msg.hpp"
#include "ThreadPool.hpp"
#include "UserSettings.hpp"



UserSettings::UserSettings(int argc,  char* argv[]) {

    executablePath = argv[0];

    // set up a vector of arguments
    std::vector<std::string> settings;
    for (int i=1; i<argc; i++)
        {
        std::string arg = argv[i];
        settings.push_back(arg);
        }

    // check the number of arguments
    if (argc == 1)
        {
        usage();
        std::cout << "   * Executable path                                       = \"" << executablePath << "\"" << std::endl;
        Msg::error("Expecting command line arguments");
        }
    if (argc % 2 == 0)
        {
        usage();
        std::cout << "   * Executable path                                       = \"" << executablePath << "\"" << std::endl;
        Msg::error("Expecting an odd number of arguments");
        }

    // initialize settings variables
    inputFile                       = "";
    treeFile                        = "";
    outputFile                      = "";
    simFile                         = "";
    numSims                         = 50;
    numChains                       = 1;
    seed                            = (uint32_t)time(NULL);
    numThreads                      = ThreadPool::defaultNumThreads();
    numAuxiliaryTables              = 10;
    temperature                     = 0.2;
    tuneTemperature                 = true;
    numMcmcCycles                   = 400000;
    burnIn                          = 0;
    tuningFrequency                 = 100;
    numGammaCategories              = 1;
    treeLengthMean                  = 1.0;
    treeLengthSD                    = 1.0;
    shapeLambda                     = 1.0;
    printFrequency                  = 1000;
    sampleFrequency                 = 100;
    isConcentrationParameterFixed   = true;
	priorConcMean                   = 3.0;
    priorMeanTables                 = 1.2;
	priorConcVariance               = 1.0;
    etAlpha                         = 1.2;
    etPi                            = 1.2;
    etTheta                         = 1.2;
    etLength                        = 1.2;
    tuningLocal                     = 3200.0;
    tuningTreeLength                = log(4.0);
    tuningBaseFrequencies           = 600.0;
    tuningBaseFrequenciesSingle     = 600.0;
    tuningExchangabilityRates       = 600.0;
    tuningExchangabilityRatesSingle = 600.0;
    tuningGammaShape                = log(4.0);
    tuningHeat                      = 0.1;
    tuningBrlen                     = 200.0;
    
    // interpret the arguments
    std::string currentArg = "";
    for (int i=0; i<settings.size(); i++)
        {
        if (currentArg == "")
            currentArg = settings[i];
        else
            {
            if (currentArg == "-i")
                inputFile = settings[i];
            else if (currentArg == "-t")
                treeFile = settings[i];
            else if (currentArg == "-o")
                outputFile = settings[i];
            else if (currentArg == "-si")
                simFile = settings[i];
            else if (currentArg == "-nSi")
                numSims = stoi(settings[i]);
            else if (currentArg == "-tu")
                tuningFrequency = stoi(settings[i]);
            else if (currentArg == "-aux")
                numAuxiliaryTables = stoi(settings[i]);
            else if (currentArg == "-threads")
                numThreads = stoi(settings[i]);
            else if (currentArg == "-seed")
                seed = (uint32_t)stoul(settings[i]);
            else if (currentArg == "-nc")
                numChains = stoi(settings[i]);
            else if (currentArg == "-temp")
                temperature = stof(settings[i]);
            else if (currentArg == "-tuneTemp")
                tuneTemperature = yesNo(settings[i]);
            else if (currentArg == "-n")
                numMcmcCycles = stoi(settings[i]);
            else if (currentArg == "-p")
                printFrequency = stoi(settings[i]);
            else if (currentArg == "-s")
                sampleFrequency = stoi(settings[i]);
            else if (currentArg == "-b")
                burnIn = stoi(settings[i]);
            else if (currentArg == "-g")
                numGammaCategories = stoi(settings[i]);
            else if (currentArg == "-lenMean")
                treeLengthMean = stof(settings[i]);
            else if (currentArg == "-lenSD")
                treeLengthSD = stof(settings[i]);
            else if (currentArg == "-e")
                shapeLambda = stof(settings[i]);
             else if (currentArg == "-c")
                isConcentrationParameterFixed = yesNo(settings[i]);
            else if (currentArg == "-m")
                priorConcMean = stof(settings[i]);
            else if(currentArg == "-mT")
                priorMeanTables = stof(settings[i]);
            else if (currentArg == "-k")
                priorConcVariance = stof(settings[i]);
            else if (currentArg == "-ePi")
                etPi = stof(settings[i]);
            else if (currentArg == "-eTheta")
                etTheta = stof(settings[i]);
            else if (currentArg == "-eShape")
                etAlpha = stof(settings[i]);
            else if (currentArg == "-eLength")
                etLength = stof(settings[i]);
            else if (currentArg == "-tLength")
                tuningTreeLength = stof(settings[i]);
            else if (currentArg == "-tLocal")
                tuningLocal = stof(settings[i]);
            else if (currentArg == "-tShape")
                tuningGammaShape = stof(settings[i]);
            else if (currentArg == "-tFreqs")
                tuningBaseFrequencies = stof(settings[i]);
            else if (currentArg == "-tFreqsSingle")
                tuningBaseFrequenciesSingle = stof(settings[i]);
            else if (currentArg == "-tRates")
                tuningExchangabilityRates = stof(settings[i]);
            else if (currentArg == "-tRatesSingle")
                tuningExchangabilityRatesSingle = stof(settings[i]);
            else if (currentArg == "-tTBR")
                tuningHeat = stof(settings[i]);
            else if (currentArg == "-tBrlen")
                tuningBrlen = stof(settings[i]);
            else
                Msg::error("Could not interpret argument " + settings[i]);
            currentArg = "";
            }
        }
    if (numAuxiliaryTables < 1)
        Msg::error("The number of auxiliary tables (-aux) must be at least 1");
}

double UserSettings::getAlphaT(void) {

    return treeLengthMean * treeLengthMean / (treeLengthSD * treeLengthSD);
}

double UserSettings::getBetaT(void) {

    return treeLengthMean / (treeLengthSD * treeLengthSD);
}

void UserSettings::print(void) {

    std::cout << std::fixed << std::setprecision(3);
    
    std::cout << "   User Settings" << std::endl;
    std::cout << "   * Executable path                                       = \"" << executablePath << "\"" << std::endl;
    std::cout << "   * Input file path and name                              = \"" << inputFile << "\"" << std::endl;
    std::cout << "   * Tree file path and name                               = \"" << treeFile << "\"" << std::endl;
    std::cout << "   * Output file path and name                             = \"" << outputFile << "\"" << std::endl;
    
    std::cout << "   * Tree length mean                                      = " << treeLengthMean << std::endl;
    std::cout << "   * Tree length standard deviation                        = " << treeLengthSD << std::endl;
    std::cout << "   * Exponential parameter for gamma shape                 = " << shapeLambda << std::endl;
    std::cout << "   * Number of gamma categories                            = " << numGammaCategories << std::endl;
    
    std::cout << "   * Concentration parameter is fixed                      = " << ((isConcentrationParameterFixed == true) ? "yes" : "no") << std::endl;
    std::cout << "   * Concentration parameter mean (when r.v.)              = " << priorConcMean << std::endl;
    std::cout << "   * Expected number of tables (when conc. is r.v.)        = " << priorMeanTables << std::endl;
    std::cout << "   * Concentration parameter variance (when r.v.)          = " << priorConcVariance << std::endl;
    std::cout << "   * Expected number of tables for tree length             = " << etLength << std::endl;
    std::cout << "   * Expected number of tables for gamma shape             = " << etAlpha << std::endl;
    std::cout << "   * Expected number of tables for base frequencies        = " << etPi << std::endl;
    std::cout << "   * Expected number of tables for exchangability rates    = " << etTheta << std::endl;
    
    std::cout << "   * Number of MCMC chains                                 = " << numChains << std::endl;
    std::cout << "   * Random number seed                                    = " << seed << std::endl;
    std::cout << "   * Number of threads                                     = " << numThreads << std::endl;
    std::cout << "   * Number of auxiliary tables (Neal's m)                 = " << numAuxiliaryTables << std::endl;
    std::cout << "   * MC^3 temperature (starting value, if tuned)           = " << temperature << std::endl;
    std::cout << "   * MC^3 temperature tuned during the burn-in             = " << ((tuneTemperature == true) ? "yes" : "no") << std::endl;
    std::cout << "   * Number of MCMC cycles                                 = " << numMcmcCycles << std::endl;
    std::cout << "   * Sample Frequency                                      = " << sampleFrequency << std::endl;
    std::cout << "   * Print Frequency                                       = " << printFrequency << std::endl;
    std::cout << "   * Burn in                                               = " << burnIn << std::endl;
    std::cout << "   * Tuning Frequency                                      = " << tuningFrequency << std::endl;
    std::cout << "   * MCMC tuning parameter for the tree topology parameter = " << tuningLocal << std::endl;
    std::cout << "   * MCMC tuning parameter for the gamma shape parameter   = " << tuningGammaShape << std::endl;
    std::cout << "   * MCMC tuning parameter for the base frequencies (D)    = " << tuningBaseFrequencies << std::endl;
    std::cout << "   * MCMC tuning parameter for the substitution rates (D)  = " << tuningExchangabilityRates << std::endl;
    std::cout << "   * MCMC tuning parameter for the base frequencies (B)    = " << tuningBaseFrequenciesSingle << std::endl;
    std::cout << "   * MCMC tuning parameter for the substitution rates (B)  = " << tuningExchangabilityRatesSingle << std::endl;
    std::cout << "   * MCMC tuning parameter for the tree length parameter   = " << tuningTreeLength << std::endl;
    std::cout << "   * MCMC tuning parameter for the TBR heat parameter      = " << tuningHeat << std::endl;
    std::cout << "   * MCMC tuning parameter for the branch proportions      = " << tuningBrlen << std::endl;
    std::cout << std::endl;
}

void UserSettings::usage(void) {

    std::cout << "   Program options for data input/output:" << std::endl;
    std::cout << "   * -i               : Input file name" << std::endl;
    std::cout << "   * -t               : Tree file name (for constraining the analysis to a fixed tree)" << std::endl;
    std::cout << "   * -o               : Output file name" << std::endl;
    std::cout << "   * -seed            : Random number seed (default: from the clock; printed at start-up)" << std::endl;
    std::cout << "   * -threads         : Number of threads (default: the number of performance cores)" << std::endl;
    std::cout << std::endl;
    
    std::cout << "   Program options for model parameters:" << std::endl;
    std::cout << "   * -lenMean         : Tree length mean" << std::endl;
    std::cout << "   * -lenSD           : Tree length standard deviation" << std::endl;
    std::cout << "   * -lenLam          : Tree length exponential parameter" << std::endl;
    std::cout << "   * -e               : Exponential parameter for shape parameter describing ASRV" << std::endl;
    std::cout << "   * -g               : Number of gamma categories" << std::endl;
    std::cout << std::endl;
    
    std::cout << "   Program options for DPP:" << std::endl;
    std::cout << "   * -c               : Concentration parameter is fixed (yes) or a random variable (no)" << std::endl;
    std::cout << "   * -aux             : Number of auxiliary tables when reseating (m in Neal's Algorithm 8; default 10)" << std::endl;
    std::cout << "   * -k               : Prior mean of the number of categories when the concentration parameter is fixed" << std::endl;
    std::cout << "   * -m               : Prior mean of the concentration parameter when it is a random variable" << std::endl;
    std::cout << "   * -mT              : Expected number of tables when the concentration parameter is a random variable" << std::endl;
    std::cout << "   * -v               : Prior variance of the concentration parameter when it is a random variable" << std::endl;
    std::cout << "   * -eLength         : Expected number of tables for tree length when the concentration parameter is fixed" << std::endl;
    std::cout << "   * -eShape          : Expected number of tables for gamma shape when the concentration parameter is fixed" << std::endl;
    std::cout << "   * -ePi             : Expected number of tables for base frequencies when the concentration parameter is fixed" << std::endl;
    std::cout << "   * -eTheta          : Expected number of tables for exchangability rates when the concentration parameter is fixed" << std::endl;
    std::cout << std::endl;
    std::cout << "   Program options for MCMC:" << std::endl;
    std::cout << "   * -n               : Number of MCMC cycles" << std::endl;
    std::cout << "   * -p               : Print-to-screen frequency" << std::endl;
    std::cout << "   * -s               : Sample-to-file frequency" << std::endl;
    std::cout << "   * -b               : The number of burn-in iterations" << std::endl;
    std::cout << "   * -tu              : How frequently to tune the tunable parameters during the burn-in" << std::endl;
    std::cout << "   * -nc              : Number of Metropolis-coupled chains (MC^3)" << std::endl;
    std::cout << "   * -temp            : MC^3 temperature (chain i is heated to power 1/(1+i*temp))" << std::endl;
    std::cout << "   * -tuneTemp        : Tune the MC^3 temperature during the burn-in (yes/no; default yes)" << std::endl;
    std::cout << "   * -tLocal          : MCMC tuning parameter for the tree topology parameter" << std::endl;
    std::cout << "   * -tTBR            : MCMC tuning parameter for heat parameter for biased TBR proposals" << std::endl;
    std::cout << "   * -tBrlen          : MCMC tuning parameter for the branch proportions update" << std::endl;
    std::cout << "   * -tShape          : MCMC tuning parameter for the gamma shape parameter" << std::endl;
    std::cout << "   * -tFreqs          : MCMC tuning parameter for the base frequencies (Dirichlet Simplex)" << std::endl;
    std::cout << "   * -tRates          : MCMC tuning parameter for the substitution rates (Dirichlet Simplex)" << std::endl;
    std::cout << "   * -tFreqsSingle    : MCMC tuning parameter for the base frequencies (Beta Simplex)" << std::endl;
    std::cout << "   * -tRatesSingle    : MCMC tuning parameter for the substitution rates (Beta Simplex)" << std::endl;
    std::cout << "   * -tLength         : MCMC tuning parameter for the tree length parameter" << std::endl;
}

bool UserSettings::yesNo(std::string str) {

    if ( (str[0] == 'Y' || str[0] == 'y') && (str[1] == 'E' || str[1] == 'e') && (str[2] == 'S' || str[2] == 's') )
        return true;
    else if ( (str[0] == 'N' || str[0] == 'n') && (str[1] == 'O' || str[1] == 'o') )
        return false;
    Msg::error("Unknown string " + str);
    return false;
}

