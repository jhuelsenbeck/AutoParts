#ifndef ParameterGammaShape_hpp
#define ParameterGammaShape_hpp

#include <vector>
#include "Parameter.hpp"
#include "SubsetParameters.hpp"
class Model;
class UserSettings;



class ParameterGammaShape : public Parameter {

    public:
                                ParameterGammaShape(void) = delete;
                                ParameterGammaShape(const ParameterGammaShape& parm);
                                ParameterGammaShape(Model* m, UserSettings* s, double lam, int nc);
        void                    accept(void);
        std::string             getHeader(void);
        std::vector<double>&    getRates(void) { return gammaRates[0]; }
        std::vector<double>     getValue(void);
        std::string             getValuesAsString(int precision);
        double                  lnProbability(void);
        Parameter*              newRandomInstance(void);
        void                    print(void);
        void                    reject(void);
        std::string             type(void);
        void                    setSubsetParameters(SubsetParameters& p) { p.gammaCategoryRates = &gammaRates[0]; }
        double                  update(void);

    private:
        double                  lambda;
        double                  shape[2];
        int                     numGammaCategories;
        std::vector<double>     gammaRates[2];
};

#endif
