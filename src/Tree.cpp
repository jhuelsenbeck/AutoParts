#include <cmath>
#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "Branch.hpp"
#include "BranchFactory.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "StateSets.hpp"
#include "Tree.hpp"



Tree::Tree(Tree& t) {

    clone(t);
}

Tree::Tree(Model* m, std::vector<std::string> tn) {

    modelPtr = m;

    // initialize the tree randomly with a three-way split at the root

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    // check size
    if (tn.size() < 3)
        Msg::error("Tree is too small");
        
    // initialize some important variables
    setNumLeaves( (int)tn.size() );
    setTaxonNames( tn );

    // set up the initial tree
    Node* r = addNode();
    setRoot(r);
    std::vector<Node*> possibleAttachmentNodes;
    for (int i=0; i<3; i++)
        {
        Node* p = addNode();
        possibleAttachmentNodes.push_back(p);
        p->setName(tn[i]);
        p->setIndex(i);
        p->setIsLeaf(true);
        p->addNeighbor(r);
        r->addNeighbor(p);
        p->setAncestor(r);
        }

    for (int i=3; i<tn.size(); i++)
        {
        Node* p = possibleAttachmentNodes[ (int)(rng.uniformRv()*possibleAttachmentNodes.size()) ];
        
        Node* newInt = addNode();
        Node* newTip = addNode();
        newTip->setIsLeaf(true);
        newTip->setName(tn[i]);
        newTip->setIndex(i);

        possibleAttachmentNodes.push_back(newInt);
        possibleAttachmentNodes.push_back(newTip);

        if (p == getRoot())
            {
            Msg::error("Selected root node for random attachment point");
            }
        else
            {
            // add the branch on another branch
            Node* pAnc = p->getAncestor();
            p->removeNeighbor(pAnc);
            pAnc->removeNeighbor(p);
            newInt->addNeighbor(p);
            newInt->addNeighbor(pAnc);
            newInt->addNeighbor(newTip);
            newInt->setAncestor(pAnc);
            p->addNeighbor(newInt);
            p->setAncestor(newInt);
            pAnc->addNeighbor(newInt);
            newTip->addNeighbor(newInt);
            newTip->setAncestor(newInt);
            }
        }
        
    // get the down pass sequence
    initalizeDownPassSequence();

    // re-index interior nodes
    int intIdx = (int)tn.size();
    for (Node* p : downPassSequence)
        {
        if (p->getIsLeaf() == false)
            p->setIndex(intIdx++);
        }
    
    // add the branches to the tree and initialize them from flat Dirichlet
    double treeLength = 0.0;
    for (Node* p : downPassSequence)
        {
        Node* pAnc = p->getAncestor();
        if (pAnc != NULL)
            {
            Branch* b = addBranch(p, pAnc, Probability::Exponential::rv(&rng, 1.0));
            treeLength += b->getProportion();
            }
        }
    for (std::map<NodePair,Branch*,CompNodePair>::iterator it = branches.begin(); it != branches.end(); it++)
        {
        double x = it->second->getProportion() / treeLength;
        it->second->setProportion(x);
        }
    if ( fabs(1.0 - branchProportionSum()) > 0.0001 )
        Msg::error("Problem setting branch lengths in random tree constructor (" + std::to_string(branchProportionSum()) + ")");

    rootOnTaxon(0);
}

Tree::Tree(Model* m, std::vector<std::string> tn, std::string newickStr) {

    modelPtr = m;
    
    taxonNames = tn;
    // break the newick string into individual tokens
    std::vector<std::string> newickTokens = parseNewickString(newickStr);

    // read the vector of Newick tokens one at a time to construct the tree
    Node* p = NULL;
    bool readingBranchLength = false;
    numLeaves = 0;
    for (int i=0; i<newickTokens.size(); i++)
        {
        std::string token = newickTokens[i];
        
        if (token == "(")
            {
            // add a node
            if (p == NULL)
                {
                // this is the first node of the tree
                p = addNode();
                root = p;
                }
            else
                {
                Node* newNode = addNode();
                p->addNeighbor(newNode);
                newNode->addNeighbor(p);
                newNode->setAncestor(p);
                p = newNode;
                addBranch(p, p->getAncestor());
                }
            readingBranchLength = false;
            }
        else if (token == ")" || token == ",")
            {
            p = p->getAncestor();
            readingBranchLength = false;
            }
        else if (token == ":")
            {
            readingBranchLength = true;
            }
        else if (token == ";")
            {
            // check that the tree is complete
            if (p != root)
                Msg::error("Ill-formatted Newick string");
            }
        else
            {
            // add a tip node or set a branch length
            if (readingBranchLength == false)
                {
                // add tip
                Node* newNode = addNode();
                newNode->setIsLeaf(true);
                p->addNeighbor(newNode);
                newNode->addNeighbor(p);
                newNode->setAncestor(p);
                newNode->setName(token);
                p = newNode;
                addBranch(p, p->getAncestor());
                numLeaves++;
                
                bool foundTaxon = false;
                for (int i=0; i<taxonNames.size(); i++)
                    {
                    if (taxonNames[i] == token)
                        {
                        p->setIndex(i);
                        foundTaxon = true;
                        break;
                        }
                    }
                if (foundTaxon == false)
                    Msg::error("Could not find taxon \"" + token + "\" in list of taxon names");
                }
            else
                {
                // assign branch length
                double x = std::stod(token);
                Branch* b = findBranch(p, p->getAncestor());
                if (b == NULL)
                    Msg::error("Could not find branch when making tree from Newick string");
                b->setProportion(x);
                readingBranchLength = false;
                }
            }
        }

    // initialize the downpass sequence
    initalizeDownPassSequence();

    // re-index interior nodes
    std::vector<Node*>& dps = getDownPassSequence();
    int intIdx = (int)taxonNames.size();
    for (Node* p : dps)
        {
        if (p->getIsLeaf() == false)
            p->setIndex(intIdx++);
        }
        
    // initialize the branch proportions
    double treeLength = 0.0;
    for (std::map<NodePair,Branch*,CompNodePair>::iterator it = branches.begin(); it != branches.end(); it++)
        treeLength += it->second->getProportion();
    for (std::map<NodePair,Branch*,CompNodePair>::iterator it = branches.begin(); it != branches.end(); it++)
        {
        double x =  it->second->getProportion();
        it->second->setProportion(x/treeLength);
        }
        
    //print();
}

Tree::~Tree(void) {

    deleteAllNodes();
    removeAllBranches();
}

Tree& Tree::operator=(Tree& t) {

    if (this != &t)
        clone(t);
    return *this;
}

Branch* Tree::addBranch(Node* e1, Node* e2) {

    BranchFactory& bf = BranchFactory::branchFactoryInstance();
    Branch* branch = bf.getBranch();
    branch->setEnds(e1, e2);
    branches.insert( std::make_pair(NodePair(e1, e2), branch) );
    return branch;
}

Branch* Tree::addBranch(Node* e1, Node* e2, double x) {

    Branch* b = addBranch(e1, e2);
    if (b != NULL)
        b->setProportion(x);
    return b;
}

Node* Tree::addNode(void) {

    Node* n = new Node;
    n->setOffset((int)nodes.size());
    nodes.push_back(n);
    return n;
}

double Tree::branchProportionSum(void) {

    double sum = 0.0;
    for (std::map<NodePair,Branch*,CompNodePair>::iterator it = branches.begin(); it != branches.end(); it++)
        sum += it->second->getProportion();
    return sum;
}

void Tree::calculateReconnectionProbabilities(CutInfo& info, std::map<NodePair,double,CompNodePair>& vals, double heat) {

    StateSets* ssPtr = modelPtr->getStateSetPtr();
    Alignment* alnPtr = modelPtr->getAlignment();
	int length1 = ssPtr->initializeStateSets(info.subtree1);
	int length2 = ssPtr->initializeStateSets(info.subtree2);

    double lambda = (double)(nodes.size()-1) / modelPtr->getTreeLength();

	double z = (lambda / (4.0/3.0 + lambda));
	double a0 = log(0.25 + 0.75 * z);
	double a1 = log(0.25 - 0.25 * z);
	double marginalLikelihoodTerm = alnPtr->getNumSites() * ((2 * getNumTaxa() - 3) * a0 - log(4.0) );
	double tempFactor = a1 - a0;

    // calculate the unscaled probabilities of all reconnection possibilities
    int nNodes1 = (int)info.subtree1.size();
    int nNodes2 = (int)info.subtree2.size();
    double maxP = 0.0;
    for (int i=0; i<info.subtree1.size(); i++)
        {
        Node* p1 = info.subtree1[i];
        Node* pAnc1 = p1->getAncestor();
        if ( (pAnc1 == NULL && nNodes1 == 1) || pAnc1 != NULL )
            {
            for (int j=0; j<info.subtree2.size(); j++)
                {
                Node* p2 = info.subtree2[j];
                Node* pAnc2 = p2->getAncestor();
                if ( (pAnc2 == NULL && nNodes2 == 1) || pAnc2 != NULL )
                    {
                    NodePair nodePair(p1, p2);
                    
                    // parsimony thing here
                    unsigned *ss1 = ssPtr->getStsPtr( 1, p1->getIndex() );
                    unsigned *ss2 = ssPtr->getStsPtr( 1, p2->getIndex() );
                    int length = length1 + length2;
                    for (int c=0; c<alnPtr->getNumSites(); c++)
                        {
                        unsigned x = (*ss1);
                        unsigned y = (*ss2);
                        unsigned zA = (x & y);
                        if ( zA == 0 )
                            length += alnPtr->getNumSitesOfPattern( c );
                        ss1++;
                        ss2++;
                        }
					double lnP = tempFactor * length + marginalLikelihoodTerm;
                    lnP *= heat;
                    
                    if (vals.size() == 0)
                        {
                        maxP = lnP;
                        }
                    else
                        {
                        if (lnP > maxP)
                            maxP = lnP;
                        }
                    vals.insert( std::make_pair(nodePair,lnP) );
                    }
                }
            }
        }
        
    // rescale the probabilities
    double sum = 0.0;
    for (std::map<NodePair,double,CompNodePair>::iterator it = vals.begin(); it != vals.end(); it++)
        {
        if (it->second - maxP < -300.0)
            it->second = 0.0;
        else
            it->second = exp(it->second - maxP);
        sum += it->second;
        }
        
    // normalize
    for (std::map<NodePair,double,CompNodePair>::iterator it = vals.begin(); it != vals.end(); it++)
        it->second /= sum;
        
#   if 0
    int i = 0;
    sum = 0.0;
    for (std::map<NodePair,double,CompNodePair>::iterator it = vals.begin(); it != vals.end(); it++)
        {
        std::cout << std::fixed << std::setprecision(5);
        std::cout << i+1 << "  ("<< it->first.getEnd1() << " - " << it->first.getEnd2() << ") -- ";
        std::cout << it->second << std::endl;
        i++;
        sum += it->second;
        }
    std::cout << "sum = " << sum << std::endl;
#   endif
}

void Tree::chooseBackbone(RandomVariable* rng, std::vector<Node*>& backboneNodes) {

    // choose an interior branch
    Node* p = NULL;
    do {
        p = nodes[(int)(rng->uniformRv()*nodes.size())];
        } while(p->getIsLeaf() == true || p->getAncestor() == NULL);
    Node* pAnc = p->getAncestor();
        
    Node* pUp = p->chooseNeighborAtRandom(rng, pAnc);
    Node* pDn = pAnc->chooseNeighborAtRandom(rng, p);
    
    backboneNodes.push_back(pUp);
    backboneNodes.push_back(p);
    backboneNodes.push_back(pAnc);
    backboneNodes.push_back(pDn);
}

NodePair Tree::choosePair(std::map<NodePair,double,CompNodePair>& vals) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double u = rng.uniformRv();
    double sum = 0.0;
    NodePair chosenPair;
    for (std::map<NodePair,double,CompNodePair>::iterator it = vals.begin(); it != vals.end(); it++)
        {
        sum += it->second;
        if (u < sum)
            {
            chosenPair.setEnds(it->first.getEnd1(), it->first.getEnd2());
            break;
            }
        }
    return chosenPair;
}

void Tree::clone(Tree& t) {

    // copy some instance variables
    modelPtr = t.modelPtr;
    numLeaves = t.numLeaves;
    taxonNames = t.taxonNames;
    
    // make certain we have the saame number of nodes in each tree
    if (nodes.size() != t.nodes.size())
        {
        deleteAllNodes();
        for (int i=0; i<t.nodes.size(); i++)
            addNode();
        }
        
    // copy the nodes
    root = nodes[t.root->getOffset()];
    for (int i=0; i<t.nodes.size(); i++)
        {
        Node* pLft = nodes[i];
        Node* pRht = t.nodes[i];
        
        pLft->setIndex( pRht->getIndex() );
        pLft->setIsLeaf( pRht->getIsLeaf() );
        pLft->setName( pRht->getName() );
        pLft->setActiveCl( pRht->getActiveCl() );
        pLft->setActiveTp( pRht->getActiveTp() );
        pLft->setClNeedsUpdate( pRht->getClNeedsUpdate() );
        pLft->setTpNeedsUpdate( pRht->getTpNeedsUpdate() );
        
        if (pRht->getAncestor() != NULL)
            pLft->setAncestor( nodes[pRht->getAncestor()->getOffset()] );
        else
            pLft->setAncestor(NULL);
        pLft->removeNeighbors();
        NodeSet& rhtNeighbors = pRht->getNeighbors();
        for (Node* n : rhtNeighbors)
            pLft->addNeighbor( nodes[n->getOffset()] );
        }
        
    // copy the branches
    removeAllBranches();
    std::map<NodePair,Branch*,CompNodePair>& rBranches = t.getBranches();
    for (std::map<NodePair,Branch*,CompNodePair>::iterator it = rBranches.begin(); it != rBranches.end(); it++)
        {
        NodePair rPair = it->first;
        Branch* b = it->second;
        Node* e1 = nodes[ rPair.getEnd1()->getOffset() ];
        Node* e2 = nodes[ rPair.getEnd2()->getOffset() ];
        addBranch(e1, e2, b->getProportion());
        }
    
    // copy the down pass sequence
    downPassSequence.clear();
    for (int i=0; i<t.downPassSequence.size(); i++)
        downPassSequence.push_back( nodes[t.downPassSequence[i]->getOffset()] );
}

void Tree::deleteAllNodes(void) {

    for (int i=0; i<nodes.size(); i++)
        delete nodes[i];
    nodes.clear();
}

Branch* Tree::findBranch(Node* e1Idx, Node* e2Idx) {

    // a local key (not the branchKey member) so that threads can look up branches at once
    NodePair key(e1Idx, e2Idx);
    std::map<NodePair, Branch*, CompNodePair>::iterator it = branches.find(key);
    if (it != branches.end())
        return it->second;
    return NULL;
}

Node* Tree::findNodeIndexed(int idx) {

    for (int i=0; i<nodes.size(); i++)
        {
        if (nodes[i]->getIndex() == idx)
            return nodes[i];
        }
    return NULL;
}

NodePair Tree::findOriginalPair(std::map<NodePair,double,CompNodePair>& vals, CutInfo& info) {

    Node* p1 = NULL;
    Node* e1 = info.reconnectionPair1.getEnd1();
    Node* e2 = info.reconnectionPair1.getEnd2();
    if (e1 == NULL && e2 != NULL)
        p1 = e2;
    else if (e1 != NULL && e2 == NULL)
        p1 = e1;
    else if (e1 != NULL && e2 != NULL)
        {
        if (e1->getAncestor() == e2)
            p1 = e1;
        else if (e2->getAncestor() == e1)
            p1 = e2;
        else
            Msg::error("Cannot find e1 or e2 as ancestor for first subtree");
        }
    else
        Msg::error("Both ends are NULL");
        
    Node* p2 = NULL;
    e1 = info.reconnectionPair2.getEnd1();
    e2 = info.reconnectionPair2.getEnd2();
    if (e1 == NULL && e2 != NULL)
        p2 = e2;
    else if (e1 != NULL && e2 == NULL)
        p2 = e1;
    else if (e1 != NULL && e2 != NULL)
        {
        if (e1->getAncestor() == e2)
            p2 = e1;
        else if (e2->getAncestor() == e1)
            p2 = e2;
        else
            Msg::error("Cannot find e1 or e2 as ancestor for second subtree");
        }
    else
        Msg::error("Both ends are NULL");

    NodePair pair(p1, p2);
    
    std::map<NodePair,double,CompNodePair>::iterator it = vals.find(pair);
    if (it == vals.end())
        Msg::error("Could not find original pair in the map");

    return pair;
}

void Tree::flipAllActiveCls(void) {

    for (int i=0; i<nodes.size(); i++)
        nodes[i]->flipActiveCl();
}

void Tree::flipAllActiveCls(Node* p) {

    while (p != NULL)
        {
        p->flipActiveCl();
        p = p->getAncestor();
        }
}

void Tree::flipAllActiveTps(void) {

    for (int i=0; i<nodes.size(); i++)
        nodes[i]->flipActiveTp();
}

std::string Tree::getNewick(void) {

    std::stringstream ss;
    writeTree(root, ss);
    ss << ";";
    std::string newick = ss.str();
    return newick;
}

void Tree::initalizeDownPassSequence(void) {

    downPassSequence.clear();
    passDown(root, root);
}

bool Tree::isBinary(void) {

    for (Node* n : downPassSequence)
        {
        if (n->getIsLeaf() == true)
            {
            if (n->getNumNeighbors() != 1)
                return false;
            }
        else
            {
            if (n->getNumNeighbors() != 3)
                return false;
            }
        }

    return true;
}

void Tree::normalizeBranchProportions(void) {

    double sum = branchProportionSum();
    for (std::map<NodePair,Branch*,CompNodePair>::iterator it = branches.begin(); it != branches.end(); it++)
        it->second->setProportion( it->second->getProportion()/sum );
}

std::vector<std::string> Tree::parseNewickString(std::string ns) {

    std::vector<std::string> vec;
    
    std::string token = "";
    for (int i=0; i<ns.length(); i++)
        {
        char c = ns[i];
        if (c == '(' || c == ')' || c == ',' || c == ':' || c == ';')
            {
            if (token != "")
                {
                vec.push_back(token);
                }
            token = c;
            vec.push_back(token);
            token = "";
            }
        else
            {
            token += c;
            }
        }
        
#   if 0
    for (int i=0; i<vec.size(); i++)
        std::cout << "token[" << i << "] \"" << vec[i] << "\"" << std::endl;
#   endif

    return vec;
}

void Tree::passDown(Node* p, Node* from) {

    if (p != NULL)
        {
        NodeSet& neighbors = p->getNeighbors();
        for (Node* n : neighbors)
            {
            if (n != from)
                passDown(n, p);
            }
        downPassSequence.push_back(p);
        if (p != from)
            p->setAncestor(from);
        else
            p->setAncestor(NULL);
        }
}

void Tree::passDown(Node* p, Node* from, std::vector<Node*>& dp) {

    if (p != NULL)
        {
        NodeSet& neighbors = p->getNeighbors();
        for (Node* n : neighbors)
            {
            if (n != from)
                passDown(n, p, dp);
            }
        dp.push_back(p);
        if (p != from)
            p->setAncestor(from);
        else
            p->setAncestor(NULL);
        }
}

void Tree::print(void) {

    printNode(root, 0);
}

void Tree::print(std::string h) {

    std::cout << h << std::endl;
    print();
}

void Tree::printBranches(std::string header) {

    std::cout << header << std::endl;
    for (std::map<NodePair, Branch*, CompNodePair>::iterator it = branches.begin(); it != branches.end(); it++)
        std::cout << it->first.getEnd1() << " " << it->first.getEnd2() << " " << it->second->getProportion() << std::endl;
}

void Tree::printNode(Node* p, int indent) {

    if (p != NULL)
        {
        NodeSet& neighbors = p->getNeighbors();

        for (int i=0; i<indent; i++)
            std::cout << " ";
        std::cout << p->getIndex();
        std::cout << " ( ";
        for (Node* n : neighbors)
            {
            if (n == p->getAncestor())
                std::cout << "a.";
            std::cout << n->getIndex() << " ";
            }
        std::cout << ") ";
        
        Branch* b = findBranch(p, p->getAncestor());
        if (b != NULL)
            std::cout << std::fixed << std::setprecision(6) << b->getProportion() << " ";
        else
            std::cout << std::fixed << std::setprecision(6) << 0.0 << " ";

        std::cout << p->getName() << " ";
        std::cout << p->getClNeedsUpdate() << p->getActiveCl() << " " << p->getTpNeedsUpdate() << p->getActiveTp() << " ";
            
        std::cout << " " << p;
        if (p == root)
            std::cout << " <- Root";

        std::cout << std::endl;

        for (Node* n : neighbors)
            {
            if (n != p->getAncestor())
                printNode(n, indent + 2);
            }
        }
}

Branch* Tree::randomBranch(RandomVariable* rng) {

    int whichBranch = (int)(rng->uniformRv() * branches.size());
    int i = 0;
    for ( std::map<NodePair, Branch*, CompNodePair>::iterator it = branches.begin(); it != branches.end(); it++)
        {
        if (i == whichBranch)
            return it->second;
        i++;
        }
    return NULL;
}

Branch* Tree::randomInteriorBranch(RandomVariable* rng) {

    Branch* b = NULL;
    do
        {
        b = randomBranch(rng);
        } while (b->isTip() == true);
    return b;
}

double Tree::randomlyCutTree(RandomVariable* rng, CutInfo& info) {

    info.orphanedNodes.clear();
    info.root1 = NULL;
    info.root2 = NULL;
    info.subtree1.clear();
    info.subtree2.clear();
    info.missingBranchLength = 0.0;
    
    double lnP = 0.0;
        
    // pick a random branch on the tree
    Branch* b = randomBranch(rng);
    Node* p = b->getDescendantNode();
    Node* pAnc = b->getAncestralNode();
    if (pAnc != p->getAncestor())
        Msg::error("problem cutting tree");
        
    // break between p and pAnc
    p->removeNeighbor(pAnc);
    pAnc->removeNeighbor(p);
    p->setAncestor(NULL);
    info.missingBranchLength = b->getProportion();
    removeBranch(p, pAnc);
    
    // remember original connection points
    NodeSet& pNeighbors = p->getNeighbors();
    std::vector<Node*> recon;
    for (Node* n : pNeighbors)
        {
        if (n != pAnc)
            recon.push_back(n);
        }
    if (recon.size() == 0)
        info.reconnectionPair1.setEnds(p, NULL);
    else if (recon.size() == 2)
        {
        info.reconnectionPair1.setEnds(recon[0], recon[1]);
        lnP -= log(findBranch(recon[0], p)->getProportion() + findBranch(recon[1], p)->getProportion());
        }
    else
        Msg::error("Expecting 0 or 2 nodes for reconnection pair from p");
    recon.clear();
    NodeSet& pAncNeighbors = pAnc->getNeighbors();
    for (Node* n : pAncNeighbors)
        {
        if (n != p)
            recon.push_back(n);
        }
    if (recon.size() == 0)
        info.reconnectionPair2.setEnds(pAnc, NULL);
    else if (recon.size() == 2)
        {
        info.reconnectionPair2.setEnds(recon[0], recon[1]);
        lnP -= log(findBranch(recon[0], pAnc)->getProportion() + findBranch(recon[1], pAnc)->getProportion());
        }
    else
        Msg::error("Expecting 0 or 2 nodes for reconnection pair from pAnc");

    // find tips above p
    std::set<int> pTips;
    tipDescendantsFrom(p, p, pTips);
    int smallestTip1 = -1, smallestTip2 = -1;
    for (int i=0; i<getNumTaxa(); i++)
        {
        std::set<int>::iterator it = pTips.find(i);
        if (it == pTips.end())
            {
            if (smallestTip2 == -1)
                smallestTip2 = i;
            }
        else
            {
            if (smallestTip1 == -1)
                smallestTip1 = i;
            }
        }
        
    // initialize subtree 1
    info.root1 = findNodeIndexed(smallestTip1);
    passDown(info.root1, info.root1, info.subtree1);
    
    // initialize subtree 2
    info.root2 = findNodeIndexed(smallestTip2);
    passDown(info.root2, info.root2, info.subtree2);
    
    // remove superfluous nodes
    std::vector<Node*> superfluousNodes;
    for (Node* n : info.subtree1)
        {
        if (n->getNumNeighbors() == 2)
            superfluousNodes.push_back(n);
        }
    for (Node* n : info.subtree2)
        {
        if (n->getNumNeighbors() == 2)
            superfluousNodes.push_back(n);
        }
    if (superfluousNodes.size() > 2)
        Msg::error("Too many superfluous nodes in randomlyCutTree");
    for (Node* n : superfluousNodes)
        {
        if (n == info.root1 || n == info.root2)
            Msg::error("Superfluous node cannot be a root of the subtree");
        std::vector<Node*> neighbors;
        n->getNeighbors(neighbors);
        Branch* b0 = findBranch(neighbors[0], n);
        Branch* b1 = findBranch(neighbors[1], n);
        if (b0 == NULL || b1 == NULL)
            Msg::error("Could not find branch(es) for superfluous node");
        double x = b0->getProportion() + b1->getProportion();
        removeBranch(neighbors[0], n);
        removeBranch(neighbors[1], n);
        n->removeNeighbors();
        info.orphanedNodes.push_back(n);
        
        neighbors[0]->removeNeighbor(n);
        neighbors[1]->removeNeighbor(n);
        neighbors[0]->addNeighbor(neighbors[1]);
        neighbors[1]->addNeighbor(neighbors[0]);
        addBranch(neighbors[0], neighbors[1], x);
        }

    info.subtree1.clear();
    info.subtree2.clear();
    passDown(info.root1, info.root1, info.subtree1);
    passDown(info.root2, info.root2, info.subtree2);
    
    return lnP;
    
#   if 0
    std::cout << "Final subtree 1:" << std::endl;
    printNode(info.root1, 3);
    std::cout << "Final subtree 2:" << std::endl;
    printNode(info.root2, 3);
    std::cout << "Cut info:" << std::endl;
    std::cout << "   orphantedNodes = ";
    for (int i=0; i<info.orphanedNodes.size(); i++)
        std::cout << info.orphanedNodes[i]->getIndex() << " ";
    std::cout << std::endl;
    std::cout << "   root1 = " << info.root1->getIndex() << std::endl;
    std::cout << "   subtree1 = ";
    for (int i=0; i<info.subtree1.size(); i++)
        std::cout << info.subtree1[i]->getIndex() << " ";
    std::cout << std::endl;
    std::cout << "   root2 = " << info.root2->getIndex() << std::endl;
    std::cout << "   subtree2 = ";
    for (int i=0; i<info.subtree2.size(); i++)
        std::cout << info.subtree2[i]->getIndex() << " ";
    std::cout << std::endl;
    std::cout << "   missingBranchLength = " << info.missingBranchLength << std::endl;
    
    double sum1 = 0.0, sum2 = 0.0;
    for (Node* n : info.subtree1)
        {
        Node* nAnc = n->getAncestor();
        if (nAnc != NULL)
            {
            Branch* b = findBranch(n, nAnc);
            sum1 += b->getProportion();
            }
        }
    for (Node* n : info.subtree2)
        {
        Node* nAnc = n->getAncestor();
        if (nAnc != NULL)
            {
            Branch* b = findBranch(n, nAnc);
            sum2 += b->getProportion();
            }
        }
    std::cout << sum1 << " + " << sum2 << " + " << info.missingBranchLength << " = " << sum1 + sum2 + info.missingBranchLength << std::endl;
#   endif
}

double Tree::reconnect(Node* p1, Node* p2, CutInfo& info) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

    Node* pAnc1 = p1->getAncestor();
    Node* pAnc2 = p2->getAncestor();
    
    double lnP = 0.0;
    if (pAnc1 == NULL && pAnc2 != NULL)
        {
        // case 1
        if (info.orphanedNodes.size() != 1)
            Msg::error("Expecting one orphaned node");
        Node* n = info.orphanedNodes[0];
        
        Branch* b2 = findBranch(p2, pAnc2);
        if (b2 == NULL)
            Msg::error("Could not find branch b2");
        double prop2 = b2->getProportion();
        double x = prop2 * rng.uniformRv();
        double y = prop2 - x;
        double z = info.missingBranchLength;
        removeBranch(p2, pAnc2);
        
        p2->removeNeighbor(pAnc2);
        pAnc2->removeNeighbor(p2);
        n->addNeighbor(p1);
        n->addNeighbor(p2);
        n->addNeighbor(pAnc2);
        p1->addNeighbor(n);
        p2->addNeighbor(n);
        pAnc2->addNeighbor(n);
        
        addBranch(p1, n, z);
        addBranch(p2, n, x);
        addBranch(pAnc2, n, y);
        
        lnP = -log(prop2);
        }
    else if (pAnc1 != NULL && pAnc2 == NULL)
        {
        // case 2
        if (info.orphanedNodes.size() != 1)
            Msg::error("Expecting one orphaned node");
        Node* n = info.orphanedNodes[0];

        Branch* b1 = findBranch(p1, pAnc1);
        if (b1 == NULL)
            Msg::error("Could not find branch b1");
        double prop1 = b1->getProportion();
        double x = prop1 * rng.uniformRv();
        double y = prop1 - x;
        double z = info.missingBranchLength;
        removeBranch(p1, pAnc1);

        p1->removeNeighbor(pAnc1);
        pAnc1->removeNeighbor(p1);
        n->addNeighbor(p1);
        n->addNeighbor(p2);
        n->addNeighbor(pAnc1);
        p1->addNeighbor(n);
        p2->addNeighbor(n);
        pAnc1->addNeighbor(n);
        
        addBranch(p2, n, z);
        addBranch(p1, n, x);
        addBranch(pAnc1, n, y);

        lnP = -log(prop1);
        }
    else if (pAnc1 != NULL && pAnc2 != NULL)
        {
        // case 3
        if (info.orphanedNodes.size() != 2)
            Msg::error("Expecting two orphaned nodes");
        Node* n1 = info.orphanedNodes[0];
        Node* n2 = info.orphanedNodes[1];
        
        Branch* b1 = findBranch(p1, pAnc1);
        if (b1 == NULL)
            Msg::error("Could not find branch b1");
        Branch* b2 = findBranch(p2, pAnc2);
        if (b2 == NULL)
            Msg::error("Could not find branch b2");
        double prop1 = b1->getProportion();
        double prop2 = b2->getProportion();
        double x1 = prop1 * rng.uniformRv();
        double y1 = prop1 - x1;
        double x2 = prop2 * rng.uniformRv();
        double y2 = prop2 - x2;
        double z = info.missingBranchLength;
        removeBranch(p1, pAnc1);
        removeBranch(p2, pAnc2);

        p1->removeNeighbor(pAnc1);
        pAnc1->removeNeighbor(p1);
        n1->addNeighbor(p1);
        n1->addNeighbor(pAnc1);
        p1->addNeighbor(n1);
        pAnc1->addNeighbor(n1);
        addBranch(p1, n1, x1);
        addBranch(pAnc1, n1, y1);
        
        p2->removeNeighbor(pAnc2);
        pAnc2->removeNeighbor(p2);
        n2->addNeighbor(p2);
        n2->addNeighbor(pAnc2);
        p2->addNeighbor(n2);
        pAnc2->addNeighbor(n2);
        addBranch(p2, n2, x2);
        addBranch(pAnc2, n2, y2);
       
        n1->addNeighbor(n2);
        n2->addNeighbor(n1);
        addBranch(n1, n2, z);
        lnP  = -log(prop1);
        lnP += -log(prop2);
        }
    else
        Msg::error("Problem reconnecting trees");
                
    return lnP;
}

void Tree::removeBranch(Node* e1, Node* e2) {

    BranchFactory& bf = BranchFactory::branchFactoryInstance();
    branchKey.setEnds(e1, e2);
    std::map<NodePair, Branch*, CompNodePair>::iterator it = branches.find(branchKey);
    if (it == branches.end())
        Msg::error("Could not find branch in branch table");
    bf.returnToPool(it->second);
    branches.erase(branchKey);
}

void Tree::removeAllBranches(void) {

    BranchFactory& bf = BranchFactory::branchFactoryInstance();
    for (std::map<NodePair, Branch*, CompNodePair>::iterator it=branches.begin(); it != branches.end(); it++)
        bf.returnToPool(it->second);
    branches.clear();
}

void Tree::rootOnTaxon(int idx) {

    // find the tip node
    Node* p = NULL;
    for (int i=0; i<nodes.size(); i++)
        {
        if (nodes[i]->getIndex() == idx)
            p = nodes[i];
        }
    if (p == NULL)
        Msg::error("Could not find tip node with index " + std::to_string(idx));
    if (p->getIsLeaf() == false)
        Msg::error("Node indexed " + std::to_string(idx) + " should be a tip node");
        
    // find the branch
    Node* pAnc = p->getAncestor();
    if (pAnc == NULL)
        {
        // we may have a tree already root on node idx
        if (p->getNumNeighbors() != 1)
            Msg::error("Could not find ancestor of node indexed " + std::to_string(idx));
        pAnc = p->getFirstNeighbor();
        if (pAnc == NULL)
            Msg::error("Could not find ancestor of node indexed " + std::to_string(idx) + " with pAnc still null");
        }
    //std::cout << "panc = " << pAnc->getIndex() << std::endl;
    Branch * b = findBranch(p, pAnc);
    if (b == NULL)
        {
        print();
        Msg::error("Could not find branch for node indexed " + std::to_string(idx));
        }
        
    // reroot
    root = pAnc;
    initalizeDownPassSequence();
}

void Tree::setAllClUpdateFlags(bool tf) {

    for (int n=0; n<nodes.size(); n++)
        nodes[n]->setClNeedsUpdate(tf);
}

void Tree::setAllClUpdateFlags(bool tf, Node* p) {

    while (p != NULL)
        {
        p->setClNeedsUpdate(tf);
        p = p->getAncestor();
        }
}

void Tree::setAllTpUpdateFlags(bool tf) {

    for (int n=0; n<nodes.size(); n++)
        nodes[n]->setTpNeedsUpdate(tf);
}

void Tree::setAllTpUpdateFlags(bool tf, Node* p) {

    while (p != NULL)
        {
        p->setTpNeedsUpdate(tf);
        p = p->getAncestor();
        }
}

void Tree::tipDescendantsFrom(Node* p, Node* from, std::set<int>& tips) {

    if (p != NULL)
        {
        NodeSet& neighbors = p->getNeighbors();
        for (Node* n : neighbors)
            {
            if (n != from)
                tipDescendantsFrom(n, p, tips);
            }
        if (p->getIsLeaf() == true)
            tips.insert(p->getIndex());
        }
}

void Tree::writeTree(Node* p, std::stringstream& ss) {

    if (p != NULL)
        {
        if (p->getIsLeaf() == true)
            {
            //ss << p->getName();
            ss << p->getIndex() + 1;
            Branch* b = findBranch(p, p->getAncestor());
            if (b != NULL)
                ss << ":" << b->getProportion();
            }
        else
            {
            ss << "(";
            }
        NodeSet& neighbors = p->getNeighbors();
        std::vector<Node*> myDescendants;
        for (Node* n : neighbors)
            {
            if (n != p->getAncestor())
                myDescendants.push_back(n);
            }
        for (int i=0; i<(int)myDescendants.size(); i++)
            {
            writeTree(myDescendants[i], ss);
            if ( (i + 1) != (int)myDescendants.size() )
                ss << ",";
            }
        if (p->getIsLeaf() == false)
            {
            ss << ")";
            Branch* b = findBranch(p, p->getAncestor());
            if (b != NULL)
                ss << ":" << b->getProportion();
            }
        }
}
