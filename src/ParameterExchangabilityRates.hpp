#ifndef ParameterExchangabilityRates_hpp
#define ParameterExchangabilityRates_hpp

#include <vector>
#include "Parameter.hpp"
#include "SubsetParameters.hpp"
class Model;
class UserSettings;



class ParameterExchangabilityRates : public Parameter {

    public:
                                ParameterExchangabilityRates(void) = delete;
                                ParameterExchangabilityRates(const ParameterExchangabilityRates& parm);
                                ParameterExchangabilityRates(Model* m, UserSettings* s);
        void                    accept(void);
        std::string             getHeader(void);
        std::vector<double>     getValue(void);
        std::string             getValuesAsString(int precision);
        double                  lnProbability(void);
        Parameter*              newRandomInstance(void);
        void                    print(void);
        void                    reject(void);
        std::string             type(void);
        void                    setSubsetParameters(SubsetParameters& p) { p.exchangeabilityRates = &rates[0]; }
        double                  update(void);
        std::vector<double>&    getRates(void) { return rates[0]; }

    private:
        std::vector<double>     rates[2];
};

#endif
