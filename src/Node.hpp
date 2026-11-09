#ifndef Node_hpp
#define Node_hpp

#include <set>
#include <string>
#include <vector>
class Node;
class RandomVariable;

// Orders nodes by their position in the tree's node list (offset) rather than by memory
// address, so that iterating over a node's neighbors gives the same order in every run.
struct NodeLess {
    bool operator()(const Node* a, const Node* b) const;
};
typedef std::set<Node*,NodeLess> NodeSet;



class Node {

    public:
                            Node(void);
        void                addNeighbor(Node* p) { neighbors.insert(p); }
        Node*               chooseNeighborAtRandom(RandomVariable* rng, Node* excludingNode);
        void                flipActiveCl(void);
        void                flipActiveTp(void);
        int                 getActiveCl(void) { return activeCl; }
        int                 getActiveTp(void) { return activeTp; }
        Node*               getAncestor(void) { return ancestor; }
        bool                getClNeedsUpdate(void) { return clNeedsUpdate; }
        std::vector<Node*>  getDescendants(void);
        Node*               getFirstNeighbor(void) { return *(neighbors.begin()); }
        int                 getIndex(void) { return index; }
        bool                getIsLeaf(void) { return isLeaf; }
        std::string         getName(void) { return name; }
        NodeSet&            getNeighbors(void) { return neighbors; }
        void                getNeighbors(NodeSet& n);
        void                getNeighbors(std::vector<Node*>& n);
        int                 getNumNeighbors(void) { return (int)neighbors.size(); }
        int                 getOffset(void) const { return offset; }
        bool                getTpNeedsUpdate(void) { return tpNeedsUpdate; }
        bool                isDescendant(Node* p);
        void                removeNeighbor(Node* p) { neighbors.erase(p); }
        void                removeNeighbors(void) { neighbors.clear(); }
        void                setActiveCl(int x) { activeCl = x; }
        void                setActiveTp(int x) { activeTp = x; }
        void                setAncestor(Node* p) { ancestor = p; }
        void                setIndex(int x) { index = x; }
        void                setIsLeaf(bool tf) { isLeaf = tf; }
        void                setName(std::string s) { name = s; }
        void                setOffset(int x) { offset = x; }
        void                setClNeedsUpdate(bool tf) { clNeedsUpdate = tf; }
        void                setTpNeedsUpdate(bool tf) { tpNeedsUpdate = tf; }
    
    private:
        NodeSet             neighbors;
        Node*               ancestor;
        int                 index;
        bool                isLeaf;
        std::string         name;
        int                 offset;
        int                 activeCl;
        int                 activeTp;
        bool                clNeedsUpdate;
        bool                tpNeedsUpdate;
};

inline bool NodeLess::operator()(const Node* a, const Node* b) const {

    return a->getOffset() < b->getOffset();
}

#endif
