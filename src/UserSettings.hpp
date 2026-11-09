#ifndef UserSettings_hpp
#define UserSettings_hpp

#include <cstdint>
#include <string>



class UserSettings {

    public:
                        UserSettings(void) = delete;
                        UserSettings(int argc,  char* argv[]);
        double          getAlphaT(void);
        double          getBetaT(void);
        int             getBurnIn(void) { return burnIn; }
        double          getGammaShapeLambda(void) { return shapeLambda; }
        double          getExpectedNumberTreeLengthTables(void) { return etLength; }
        double          getExpectedNumberAlphaTables(void) { return etAlpha; }
        double          getExpectedNumberPiTables(void) { return etPi; }
        double          getExpectedNumberThetaTables(void) { return etTheta; }
        bool            getIsConcentrationParameterFixed(void) { return isConcentrationParameterFixed; }
        double          getPriorConcMean(void) { return priorConcMean; }
        double          getPriorMeanTables(void) { return priorMeanTables; }
        double          getPriorConcVariance(void) { return priorConcVariance; }
        std::string     getInputFile(void) { return inputFile; }
        std::string     getOutputFile(void) { return outputFile; }
        int             getNumGammaCategories(void) { return numGammaCategories; }
        int             getNumMcmcCycles(void) { return numMcmcCycles; }
        int             getPrintFrequency(void) { return printFrequency; }
        int             getSampleFrequency(void) { return sampleFrequency; }
        int             getTuningFrequency(void) { return tuningFrequency; }
        int             getNumChains(void) { return numChains; }
        uint32_t        getSeed(void) { return seed; }
        int             getNumSims(void) { return numSims; }
        std::string     getSimFile(void) { return simFile; }
        double          getTemperature(void) { return temperature; }
        std::string     getTreeFile(void) { return treeFile; }
        double          getTuningBrlen(void) { return tuningBrlen; }
        double          getTuningLocal(void) { return tuningLocal; }
        double          getTuningTreeLength(void) { return tuningTreeLength; }
        double          getTuningBaseFrequencies(void) { return tuningBaseFrequencies; }
        double          getTuningBaseFrequenciesSingle(void) { return tuningBaseFrequenciesSingle; }
        double          getTuningExchangabilityRates(void) { return tuningExchangabilityRates; }
        double          getTuningExchangabilityRatesSingle(void) { return tuningExchangabilityRatesSingle; }
        double          getTuningGammaShape(void) { return tuningGammaShape; }
        double          getTuningHeat(void) { return tuningHeat; }
        void            print(void);
        void            setInputFile(std::string s) { inputFile = s; }
        void            setNumMcmcCycles(int x) { numMcmcCycles = x; }
        void            setOutputFile(std::string s) { outputFile = s; }
        void            usage(void);
        
    private:
        uint32_t        seed;
        bool            yesNo(std::string str);
        std::string     executablePath;
        std::string     inputFile;
        std::string     treeFile;
        std::string     outputFile;
        std::string     simFile;
        int             numSims;
        int             numChains;
        int             numMcmcCycles;
        int             printFrequency;
        int             sampleFrequency;
        int             burnIn;
        int             tuningFrequency;
        int             numGammaCategories;
        double          temperature;
        double          treeLengthMean;
        double          treeLengthSD;
        double          shapeLambda;
        bool            isConcentrationParameterFixed;
        double          priorConcMean;
        double          priorMeanTables;
        double          priorConcVariance;
        double          etLength;
        double          etAlpha;
        double          etPi;
        double          etTheta;
        double          tuningLocal;
        double          tuningHeat;
        double          tuningBrlen;
        double          tuningTreeLength;
        double          tuningBaseFrequencies;
        double          tuningBaseFrequenciesSingle;
        double          tuningExchangabilityRates;
        double          tuningExchangabilityRatesSingle;
        double          tuningGammaShape;
};

#endif
