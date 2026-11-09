#ifndef Restaurant_hpp
#define Restaurant_hpp

#include <map>
#include <set>
#include <string>
#include "Table.hpp"

class Model;
class Parameter;
class UserSettings;



class Restaurant {

    public:
                            Restaurant(void) = delete;
                            Restaurant(const Restaurant& r);
                            Restaurant(Model* mp, UserSettings*, bool sf, double a, int np, Parameter* parm);
                           ~Restaurant(void);
        Table*              addTable(void);
        static double       calculateAlphaFromExpectedNumberOfTables(double expT, int np);
        static double       expectedNumberOfTables(double a, int np);
        Table*              findTableWithPatron(int idx);
        double              getConcentrationParameter(void) { return alpha; }
        int                 getNumPatrons(void) { return numPatrons; }
        int                 getNumTables(void) { return (int)tables.size(); }
        Parameter*          getParameter(void) { return parameter; }
        double              lnProbability(void);
        void                print(void);
        void                removeTable(Table* tab);
        std::string         rgf(void);
        void                setModel(Model* m) {modelPtr = m; parameter->setModel(m); for(Table* t : tables) {t->setModel(m);}}
        double              update(double power);
        
    private:
        Table*              addAuxiliaryTable(void);
        Table*              chooseTable(std::map<Table*,double,TableLess>& lnProbs);
        Parameter*          copyParameter(Parameter* parmToCopy);
        void                normalize(std::map<Table*,double,TableLess>& lnProbs);
        double              sampleAlpha(int k, int n, double oldAlpha, double a, double b);
        Model*              modelPtr;
        UserSettings*       settingsPtr;
        double              alpha;
        int                 numPatrons;
        TableSet            tables;
        Parameter*          parameter;
        bool                isSeatingRv;
        double              gammaAlpha;
        double              gammaBeta;
};

#endif
