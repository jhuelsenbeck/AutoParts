#ifndef ParameterTree_hpp
#define ParameterTree_hpp

#include "Parameter.hpp"
#include "SubsetParameters.hpp"
#include "Tree.hpp"
class Model;
class UserSettings;



class ParameterTree : public Parameter {

    public:
                                ParameterTree(void) = delete;
                                ParameterTree(const ParameterTree& parm);
                                ParameterTree(Model* m, UserSettings* s, std::vector<std::string> tn);
        void                    accept(void);
        Tree*                   getActiveTree(void) { return trees[0]; }
        std::string             getHeader(void);
        std::vector<double>     getValue(void);
        std::string             getValuesAsString(int precision);
        double                  lnProbability(void);
        Parameter*              newRandomInstance(void);
        void                    print(void);
        void                    reject(void);
        std::string             type(void);
        void                    setSubsetParameters(SubsetParameters& p) { p.tree = trees[0]; }
        double                  update(void);
        void                    setModel(Model* m){trees[0]->setModel(m); trees[1]->setModel(m);}

    private:
        double                  updateBrlen(void);
        double                  updateLocal(void);
        double                  updateTbr(void);
        Tree*                   trees[2];
};

#endif
