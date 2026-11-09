#ifndef Parameter_hpp
#define Parameter_hpp

#include <limits>
#include <string>
#include <vector>
class Model;
struct SubsetParameters;
class Table;
class UserSettings;



// Returned by update() as the log proposal ratio when the proposed value falls outside
// the parameter's support: the move is then rejected, leaving the value unchanged.
// (Drawing again until the value is inside the support would make the proposal
// asymmetric near the bounds, which the Hastings ratio does not allow for.)
inline double rejectProposal(void) { return -std::numeric_limits<double>::infinity(); }

class Parameter {

    public:
                                    Parameter(void) = delete;
                                    Parameter(const Parameter& parm);
                                    Parameter(Model* m, UserSettings* s, std::string n) { modelPtr = m; userSettings = s; name = n; }
        virtual                    ~Parameter(void) { }
        virtual void                accept(void) = 0;
        virtual std::string         getHeader(void) = 0;
        virtual std::string         getName(void) { return name; }
        virtual double              getProposalProbability(void) { return proposalProbability; }
        virtual Table*              getTable(void) { return myTable; }
        virtual bool                getUpdateModifiesEigens(void) { return updateModifiesEigens; }
        virtual std::string         getValuesAsString(int precision) = 0;
        virtual std::vector<double> getValue(void) = 0;
        virtual double              lnProbability(void) = 0;
        virtual Parameter*          newRandomInstance(void) = 0;
        virtual void                print(void) = 0;
        virtual void                reject(void) = 0;
        virtual void                setName(std::string s) { name = s; }
        virtual void                setProposalProbability(double x) { proposalProbability = x; }
        virtual void                setSubsetParameters(SubsetParameters& p) = 0;   // supply this parameter's value to p
        virtual void                setTable(Table* t) { myTable = t; }
        virtual void                setModel(Model* m) {modelPtr = m;}
        virtual void                setUpdateModifiesEigens(bool tf) { updateModifiesEigens = tf; }
        virtual std::string         type(void) = 0;
        virtual double              update(void) = 0;
    
    protected:
        bool                        anyLessThanMin(std::vector<double>& f, double minVal);
        std::string                 name;
        double                      proposalProbability;
        bool                        updateModifiesEigens;
        Model*                      modelPtr;
        Table*                      myTable;
        UserSettings*               userSettings;
};

#endif
