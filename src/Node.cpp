#include "Msg.hpp"
#include "Node.hpp"
#include "RandomVariable.hpp"



Node::Node(void) {

    ancestor = NULL;
    index = 0;
    isLeaf = false;
    offset = 0;
    name = "";
    activeCl = 0;
    activeTp = 0;
    clNeedsUpdate = false;
    tpNeedsUpdate = false;
}

Node* Node::chooseNeighborAtRandom(RandomVariable* rng, Node* excludingNode) {

    Node* p = NULL;
    
    do {
        int whichNode = rng->uniformRv() * neighbors.size();
        int i = 0;
        for (Node* n : neighbors)
            {
            if (i == whichNode)
                {
                p = n;
                break;
                }
            i++;
            }
        } while (p == excludingNode);
    
    return p;
}

std::vector<Node*> Node::getDescendants(void) {

    std::vector<Node*> des;
    for (Node* n : neighbors)
        {
        if (n != ancestor)
            {
            des.push_back(n);
            }
        }
    return des;
}

void Node::getNeighbors(NodeSet& n) {

    for (Node* p : neighbors)
        n.insert(p);
}

void Node::getNeighbors(std::vector<Node*>& n) {

    for (Node* p : neighbors)
        n.push_back(p);
}

void Node::flipActiveCl(void) {

    if (activeCl == 0)
        activeCl = 1;
    else
        activeCl = 0;
}

void Node::flipActiveTp(void) {

    if (activeTp == 0)
        activeTp = 1;
    else
        activeTp = 0;
}

bool Node::isDescendant(Node* p) {

    bool is_in = (neighbors.find(p) != neighbors.end()) && (p != ancestor);
    return is_in;
}
