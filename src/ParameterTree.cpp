#include <cmath>
#include <iostream>
#include "Branch.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "ParameterTree.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Tree.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"



ParameterTree::ParameterTree(const ParameterTree& parm) : Parameter(parm) {

    std::cout << "ParameterTree copy constructor" << std::endl;

    trees[0] = new Tree(*parm.trees[0]);
    trees[1] = new Tree(*parm.trees[1]);
}

ParameterTree::ParameterTree(Model* m, UserSettings* s, std::vector<std::string> tn) : Parameter(m, s, "Tree") {

    updateModifiesEigens = false;
    trees[0] = new Tree(modelPtr, tn);
    trees[1] = new Tree(*trees[0]);
}

void ParameterTree::accept(void) {

    *trees[1] = *trees[0];
}

std::string ParameterTree::getHeader(void) {

    return "";
}

std::vector<double> ParameterTree::getValue(void) {

    std::vector<double> vals;
    return vals;
}

std::string ParameterTree::getValuesAsString(int precision) {

    std::string newickStr = trees[0]->getNewick();
    return newickStr;
}

double ParameterTree::lnProbability(void) {

    return 2 * trees[0]->getNumTaxa() - 2;
}

Parameter* ParameterTree::newRandomInstance(void) {

    return new ParameterTree(modelPtr, userSettings, trees[0]->getTaxonNames());
}

void ParameterTree::print(void) {

    trees[0]->print();
}

void ParameterTree::reject(void) {

    *trees[0] = *trees[1];
}

std::string ParameterTree::type(void) {

    return "tree";
}

double ParameterTree::update(void) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double u = rng.uniformRv();
    if (u <= 0.33)
        return updateTbr();
    else if (u > 0.33 && u <= 0.67)
        return updateBrlen();
    else
        return updateLocal();
}

double ParameterTree::updateBrlen(void) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double alpha0 = UpdateInfo::updateInfo().attempt(UpdateType::BRANCH_PROPORTIONS);
    Tree* t = trees[0];

    // choose a branch at random
    Branch* b = t->randomBranch(&rng);

	// propose new proportions
	double oldP = b->getProportion();
	std::vector<double> alphaForward(2);
	std::vector<double> alphaReverse(2);
	std::vector<double> oldProportions(2);
	std::vector<double> newProportions(2);
    oldProportions[0] = oldP;
    oldProportions[1] = 1.0 - oldP;
    for (int i=0; i<2; i++)
        alphaForward[i] = oldProportions[i] * alpha0;
    bool err = Probability::Dirichlet::rv(&rng, alphaForward, newProportions);
    if (err == true || anyLessThanMin(newProportions, MIN_PROPORTION) == true)
        return rejectProposal();     // outside the support (see Parameter.hpp)
    double newP = newProportions[0];
    double lnForwardProb = Probability::Dirichlet::lnPdf(alphaForward, newProportions);

    for (int i=0; i<2; i++)
        alphaReverse[i] = newProportions[i] * alpha0;
    double lnReverseProb = Probability::Dirichlet::lnPdf(alphaReverse, oldProportions);
	
	// update the branch lengths
    double factor = (1.0 - newP) / (1.0 - oldP);
    double sum = 0.0;
    std::map<NodePair,Branch*,CompNodePair>& brlens = t->getBranches();
    for (std::map<NodePair,Branch*,CompNodePair>::iterator it = brlens.begin(); it != brlens.end(); it++)
        {
        if (it->second != b)
            {
            double p = it->second->getProportion();
            p *= factor;
            it->second->setProportion(p);
            }
        else
            {
            it->second->setProportion(newP);
            }
        sum += it->second->getProportion();
        }
    
    if ( fabs(1.0-sum) > 0.0001 )
        Msg::error("Problem with branch proportions in updateBrlen " + std::to_string(sum));
    
	return (lnReverseProb - lnForwardProb) + (t->getNumNodes()-2) * log(factor);
}

double ParameterTree::updateLocal(void) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    double alpha0 = UpdateInfo::updateInfo().attempt(UpdateType::LOCAL);
    Tree* t = trees[0];

    // choose the backbone and modify the branch lengths randomly
    std::vector<Node*> backboneNodes;
    t->chooseBackbone(&rng, backboneNodes);
    std::vector<Branch*> backboneBranches;
    for (int i=0; i<3; i++)
        {
        Branch* b = t->findBranch( backboneNodes[i], backboneNodes[i+1] );
        if (b == NULL)
            Msg::error("Can't find branch in LOCAL");
        backboneBranches.push_back(b);
        }
        
    // change the length of the backbone
    std::vector<double> oldP(2);
    std::vector<double> alphaForward(2);
    oldP[0] = 0.0;
    for (int i=0; i<3; i++)
        {
        Branch* b = backboneBranches[i];
        oldP[0] += b->getProportion();
        }
    oldP[1] = 1.0 - oldP[0];
    alphaForward[0] = oldP[0] * alpha0;
    alphaForward[1] = oldP[1] * alpha0;
    std::vector<double> newP(2);
    bool err = Probability::Dirichlet::rv(&rng, alphaForward, newP);
    if (err == true || anyLessThanMin(newP, 0.0001) == true)
        return rejectProposal();     // outside the support (see Parameter.hpp)
    double newV[3];
    for (int i=0; i<3; i++)
        newV[i] = backboneBranches[i]->getProportion() * (newP[0]/oldP[0]);
    std::vector<double> alphaReverse(2);
    alphaReverse[0] = newP[0] * alpha0;
    alphaReverse[1] = newP[1] * alpha0;
    std::map<NodePair,Branch*,CompNodePair>& branches = t->getBranches();
    for (std::map<NodePair,Branch*,CompNodePair>::iterator it = branches.begin(); it != branches.end(); it++)
        {
        Branch* b = it->second;
        double x = b->getProportion();
        b->setProportion(x * (newP[1]/oldP[1]));
        }
    for (int i=0; i<3; i++)
        backboneBranches[i]->setProportion(newV[i]);
    //t->normalizeBranchProportions();

	double lnForwardProb = Probability::Dirichlet::lnPdf(alphaForward, newP);
	double lnReverseProb = Probability::Dirichlet::lnPdf(alphaReverse, oldP);
	double lnProposalProb = lnReverseProb - lnForwardProb;
	lnProposalProb += 2.0 * (log(newP[0]) - log(oldP[0])) + (branches.size()-3) * (log(newP[1]) - log(oldP[1]));

    // update tree
    if (rng.uniformRv() < 0.5)
        {
        double newPt = rng.uniformRv()*newP[0];
        if (newPt < backboneBranches[0]->getProportion() + backboneBranches[1]->getProportion())
            {
            // no change in topology
            double x = backboneBranches[0]->getProportion() + backboneBranches[1]->getProportion();
            backboneBranches[0]->setProportion(newPt);
            backboneBranches[1]->setProportion(x-newPt);
            }
        else
            {
            // topology changes
            Node* q = backboneNodes[1];
            q->removeNeighbor(backboneNodes[0]);
            q->removeNeighbor(backboneNodes[2]);
            backboneNodes[0]->removeNeighbor(q);
            backboneNodes[2]->removeNeighbor(q);
            backboneNodes[0]->addNeighbor(backboneNodes[2]);
            backboneNodes[2]->addNeighbor(backboneNodes[0]);
            
            backboneNodes[2]->removeNeighbor(backboneNodes[3]);
            backboneNodes[3]->removeNeighbor(backboneNodes[2]);
            backboneNodes[2]->addNeighbor(q);
            backboneNodes[3]->addNeighbor(q);
            q->addNeighbor(backboneNodes[2]);
            q->addNeighbor(backboneNodes[3]);
            
            double x = backboneBranches[0]->getProportion() + backboneBranches[1]->getProportion();
            //double y = backboneBranches[2]->getLength();
            t->removeBranch(backboneNodes[0], backboneNodes[1]);
            t->removeBranch(backboneNodes[1], backboneNodes[2]);
            t->removeBranch(backboneNodes[2], backboneNodes[3]);
            t->addBranch(backboneNodes[0], backboneNodes[2], x);
            t->addBranch(backboneNodes[2], backboneNodes[1], newP[0]-newPt);
            t->addBranch(backboneNodes[1], backboneNodes[3], newPt-x);
            
            t->initalizeDownPassSequence();
            }
        }
    else
        {
        double newPt = rng.uniformRv()*newP[0];
        if (newPt > backboneBranches[0]->getProportion())
            {
            // no change in topology
            double x = backboneBranches[0]->getProportion();
            backboneBranches[1]->setProportion(newPt-x);
            backboneBranches[2]->setProportion(newP[0]-newPt);
            }
        else
            {
            // topology changes
            Node* q = backboneNodes[2];
            q->removeNeighbor(backboneNodes[1]);
            q->removeNeighbor(backboneNodes[3]);
            backboneNodes[1]->removeNeighbor(q);
            backboneNodes[3]->removeNeighbor(q);
            backboneNodes[1]->addNeighbor(backboneNodes[3]);
            backboneNodes[3]->addNeighbor(backboneNodes[1]);
            
            backboneNodes[0]->removeNeighbor(backboneNodes[1]);
            backboneNodes[1]->removeNeighbor(backboneNodes[0]);
            backboneNodes[0]->addNeighbor(q);
            backboneNodes[1]->addNeighbor(q);
            q->addNeighbor(backboneNodes[0]);
            q->addNeighbor(backboneNodes[1]);
            
            double x = backboneBranches[1]->getProportion() + backboneBranches[2]->getProportion();
            double y = backboneBranches[0]->getProportion();
            t->removeBranch(backboneNodes[0], backboneNodes[1]);
            t->removeBranch(backboneNodes[1], backboneNodes[2]);
            t->removeBranch(backboneNodes[2], backboneNodes[3]);
            t->addBranch(backboneNodes[0], backboneNodes[2], newPt);
            t->addBranch(backboneNodes[2], backboneNodes[1], y-newPt);
            t->addBranch(backboneNodes[1], backboneNodes[3], x);

            //t->rootOnTaxon(0);
            t->initalizeDownPassSequence();
            }
        }

    t->rootOnTaxon(0);
    t->normalizeBranchProportions();
    
    // perform a few correctness checks
    if (t->isBinary() == false)
        Msg::error("Proposed TBR tree is not binary");
    if ( fabs(1.0 - t->branchProportionSum()) > 0.0001 )
        Msg::error("Branch proportions are wrong in LOCAL " + std::to_string(t->branchProportionSum()));
    
    return lnProposalProb;
}

double ParameterTree::updateTbr(void) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();
    Tree* t = trees[0];
    double lnProposalProbability = 0.0;
    double heat = UpdateInfo::updateInfo().attempt(UpdateType::TBR);
    
    // cut the tree
    CutInfo info;
    lnProposalProbability -= t->randomlyCutTree(&rng, info);

    // calculate the scores of all potential reconnection points
    std::map<NodePair,double,CompNodePair> vals;
    t->calculateReconnectionProbabilities(info, vals, heat);
    
    // choose a reconnection point
    NodePair newReconnectionPair = t->choosePair(vals);
    
    // find the old reconnection point
    NodePair oldReconnectionPair = t->findOriginalPair(vals, info);
    
    // reconnect tree
    lnProposalProbability += t->reconnect(newReconnectionPair.getEnd1(), newReconnectionPair.getEnd2(), info);
    
    // modify the proposal probability by the connection probabilities
    std::map<NodePair,double,CompNodePair>::iterator it = vals.find(oldReconnectionPair);
    if (it == vals.end())
        Msg::error("Could not find old connection point in map");
    lnProposalProbability += log(it->second);
    it = vals.find(newReconnectionPair);
    if (it == vals.end())
        Msg::error("Could not find old connection point in map");
    lnProposalProbability -= log(it->second);
        
    // reroot tree in standard way
    t->rootOnTaxon(0);
    
    // perform a few correctness checks
    if (t->isBinary() == false)
        Msg::error("Proposed TBR tree is not binary");
    if ( fabs(1.0 - t->branchProportionSum()) > 0.0001 )
        Msg::error("Branch proportions are wrong " + std::to_string(t->branchProportionSum()));
    
    return lnProposalProbability;
}
