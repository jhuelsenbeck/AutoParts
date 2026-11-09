#ifndef ParameterBaseFrequencies_hpp
#define ParameterBaseFrequencies_hpp

#include <vector>
#include "Parameter.hpp"
#include "SubsetParameters.hpp"
class Model;
class UserSettings;



class ParameterBaseFrequencies : public Parameter {

    public:
                                ParameterBaseFrequencies(void) = delete;
                                ParameterBaseFrequencies(const ParameterBaseFrequencies& parm);
                                ParameterBaseFrequencies(Model* m, UserSettings* s);
        void                    accept(void);
        std::vector<double>&    getFreqs(void) { return freqs[0]; }
        std::string             getHeader(void);
        std::vector<double>     getValue(void);
        std::string             getValuesAsString(int precision);
        double                  lnProbability(void);
        Parameter*              newRandomInstance(void);
        void                    print(void);
        void                    reject(void);
        std::string             type(void);
        void                    setSubsetParameters(SubsetParameters& p) { p.baseFrequencies = &freqs[0]; }
        double                  update(void);

    private:
        std::vector<double>     freqs[2];
};

#endif
