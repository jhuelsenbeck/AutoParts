#ifndef StateSets_H
#define StateSets_H

#include <vector>

class Alignment;
class UserSettings;
class Node;

class StateSets {

	public:
                        StateSets(Alignment* ap, UserSettings* sp);
                       ~StateSets(void);
        unsigned*       getStsPtr(int space, int node) { return stsPtr[space][node]; }
        const int*      getPatternWeights(void) const { return patternWeights.data(); }
        int             getNumChar(void) const { return numChar; }
        int             calcStateSetFor(int s1, int n1, int s2, int n2, int s3, int n3);
        int             initializeStateSets(std::vector<Node*>& dp);

    private:
        Alignment*      alignmentPtr;
        UserSettings*   settingsPtr;
        unsigned*       sts;
        unsigned***     stsPtr;
        int**           lengths;
        unsigned**      codedMatrix;
        int             numTaxa;
        int             numChar;
        int             numNodes;
        int             oneStateSetSize;
        std::vector<int> patternWeights;    // number of sites with each pattern
};

#endif
