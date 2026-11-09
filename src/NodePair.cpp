#include <iostream>
#include "NodePair.hpp"



NodePair::NodePair(void) {

    e1 = NULL;
    e2 = NULL;
}

NodePair::NodePair(Node* n1, Node* n2) {

    if (nodePairKey(n1) < nodePairKey(n2))
        {
        e1 = n1;
        e2 = n2;
        }
    else
        {
        e1 = n2;
        e2 = n1;
        }
}

NodePair::NodePair(const Node* n1, const Node* n2) {

    if (nodePairKey(n1) < nodePairKey(n2))
        {
        e1 = (Node*)n1;
        e2 = (Node*)n2;
        }
    else
        {
        e1 = (Node*)n2;
        e2 = (Node*)n1;
        }
}

void NodePair::print(void) {

    std::cout << "(" << e1 << "," << e2 << ")" << std::endl;
}

void NodePair::setEnds(Node* n1, Node* n2) {

    if (nodePairKey(n1) < nodePairKey(n2))
        {
        e1 = n1;
        e2 = n2;
        }
    else
        {
        e1 = n2;
        e2 = n1;
        }
}


