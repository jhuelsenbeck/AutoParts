#ifndef NodePair_hpp
#define NodePair_hpp

#include "Node.hpp"



class NodePair {

    public:
                NodePair(void);
                NodePair(Node* n1, Node* n2);
                NodePair(const Node* n1, const Node* n2);
        Node*   getEnd1(void) { return e1; }
        Node*   getEnd2(void) { return e2; }
        Node*   getEnd1(void) const { return e1; }
        Node*   getEnd2(void) const { return e2; }
        void    print(void);
        void    setEnds(Node* n1, Node* n2);
        
    private:
        Node*   e1;
        Node*   e2;
};

// Node pairs are ordered by the ends' positions in the tree (offsets), not by memory
// address, so that iterating over branches gives the same order in every run.
inline int nodePairKey(const Node* n) { return (n == NULL) ? -1 : n->getOffset(); }

struct CompNodePair {

    bool operator()(const NodePair& np1, const NodePair& np2) const {
        
        int a1 = nodePairKey(np1.getEnd1()), b1 = nodePairKey(np2.getEnd1());
        if (a1 != b1)
            return a1 < b1;
        return nodePairKey(np1.getEnd2()) < nodePairKey(np2.getEnd2());
        }
};

#endif
