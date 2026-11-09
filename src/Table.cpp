#include "Parameter.hpp"
#include "Table.hpp"



Table::Table(void) {

    parm = NULL;
    myRestaurant = NULL;
    serial = 0;
}

void Table::clean(void) {

    myRestaurant = NULL;
    delete parm;
    parm = NULL;
    patrons.clear();
}

bool Table::hasPatron(int idx) {

    std::set<int>::iterator it = patrons.find(idx);
    if (it != patrons.end())
        return true;
    return false;
}
