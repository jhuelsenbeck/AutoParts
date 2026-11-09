#ifndef Table_hpp
#define Table_hpp

#include <set>
#include "Model.hpp"
#include "Parameter.hpp"

class Restaurant;


class Table {

    public:
                        Table(void);
        void            addParameter(Parameter* p) { parm = p; }
        void            addPatron(int idx) { patrons.insert(idx); }
        void            clean(void);
        int             getNumPatrons(void) { return (int)patrons.size(); }
        Parameter*      getParameter(void) { return parm; }
        std::set<int>&  getPatrons(void) { return patrons; }
        Restaurant*     getRestaurant(void) { return myRestaurant; }
        bool            hasPatron(int idx);
        void            removeAllPatrons(void) { patrons.clear(); }
        void            removePatron(int idx) { patrons.erase(idx); }
        void            setRestaurant(Restaurant* r) { myRestaurant = r; }
        void            setModel(Model* m) {parm->setModel(m);}
        long            getSerial(void) const { return serial; }
        void            setSerial(long x) { serial = x; }
    
    private:
        std::set<int>   patrons;
        Parameter*      parm;
        Restaurant*     myRestaurant;
        long            serial;
};

// Orders tables by the order in which they were created rather than by memory address,
// so that iterating over tables gives the same order in every run.
struct TableLess {
    bool operator()(const Table* a, const Table* b) const { return a->getSerial() < b->getSerial(); }
};
typedef std::set<Table*,TableLess> TableSet;

#endif
