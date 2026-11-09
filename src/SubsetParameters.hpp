#ifndef SubsetParameters_hpp
#define SubsetParameters_hpp

#include <vector>
class Tree;



// The parameter values used to calculate the likelihood of one data subset. The
// pointers refer to values held by the parameters themselves; nothing is copied.
// Parameter::setSubsetParameters() replaces the value that a parameter supplies, which
// is how a subset's likelihood is calculated as if it were seated at another table.

struct SubsetParameters {

    Tree*                   tree;
    double                  treeLength;
    std::vector<double>*    baseFrequencies;
    std::vector<double>*    exchangeabilityRates;
    std::vector<double>*    gammaCategoryRates;
};

#endif
