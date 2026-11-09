#include <iomanip>
#include <iostream>
#include "Alignment.hpp"
#include "Chunk.hpp"
#include "ConditionalLikelihoods.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "TransitionProbabilities.hpp"
#include "Tree.hpp"
#include <unordered_map>

#undef NO_DATA

#define MAX_GAMMA_CATEGORIES 64


Chunk::Chunk(const Chunk& c) {

    alignment = c.alignment;
    model = c.model;
    id = c.id;
    numSites = c.numSites;
    numGammaCategories = c.numGammaCategories;

    condLikes = new ConditionalLikelihoods(alignment, numGammaCategories);
    transProbs = new TransitionProbabilities(model, this, 2*alignment->getNumTaxa()-1, numGammaCategories);
    initializeTipStates();
}

Chunk::Chunk(int myId, Alignment* a, Model* m, int nc) {

    alignment = a;
    model = m;
    id = myId;
    
    numSites = alignment->getNumSites();
    numGammaCategories = nc;
    
    condLikes = new ConditionalLikelihoods(alignment, numGammaCategories);
    transProbs = new TransitionProbabilities(model, this, 2*alignment->getNumTaxa()-1, numGammaCategories);
    initializeTipStates();
}

Chunk::~Chunk(void) {

    delete alignment;
    delete condLikes;
    delete transProbs;
}

int Chunk::getNumSitesForChunk(void) {

    return alignment->getNumSites();
}

void Chunk::flipActiveTransitionProbabilities(void) {

    transProbs->flipActiveProbs();
}

void Chunk::flipActiveEigens(void) {
    
    transProbs->flipActiveEigens();
}

namespace {

    // The conditional likelihoods are stored as [site][gamma category][state] and the
    // transition probability matrices are column major, so column j of matrix T is
    // T[4*j], ..., T[4*j+3]. Each function below forms, for one child of node p,
    //    x_i = sum_j P(i,j) * cl_child[j]
    // and either stores it in (FIRST == true) or multiplies it into the conditional
    // likelihoods of p. The arithmetic is the same, term for term, as computing all of
    // the children's sums and multiplying them together in one pass.

    template<bool FIRST>
    inline void applyInternalChild(double* clP, const double* cl, const double* const* T, int numSites, int numCats) {

        for (int c=0; c<numSites; c++)
            {
            for (int k=0; k<numCats; k++)
                {
                const double* M = T[k];
                double x0 = cl[0], x1 = cl[1], x2 = cl[2], x3 = cl[3];
                for (int i=0; i<4; i++)
                    {
                    double sum = M[i]*x0 + M[4+i]*x1 + M[8+i]*x2 + M[12+i]*x3;
                    if (FIRST == true)
                        clP[i] = sum;
                    else
                        clP[i] *= sum;
                    }
                clP += 4;
                cl  += 4;
                }
            }
    }

    template<bool FIRST>
    inline void applyTipChild(double* clP, const double* cl, const int* states, const double* const* T, int numSites, int numCats) {

        for (int c=0; c<numSites; c++)
            {
            int s = states[c];
            if (s >= 0)
                {
                // unambiguous tip: the sum over j picks out column s of the matrix
                for (int k=0; k<numCats; k++)
                    {
                    const double* col = T[k] + 4*s;
                    for (int i=0; i<4; i++)
                        {
                        if (FIRST == true)
                            clP[i] = col[i];
                        else
                            clP[i] *= col[i];
                        }
                    clP += 4;
                    }
                cl += 4 * numCats;
                }
            else
                {
                applyInternalChild<FIRST>(clP, cl, T, 1, numCats);
                clP += 4 * numCats;
                cl  += 4 * numCats;
                }
            }
    }
}

void Chunk::initializeTipStates(void) {

    int numTaxa = alignment->getNumTaxa();
    tipStates.assign(numTaxa * numSites, -1);
    for (int i=0; i<numTaxa; i++)
        {
        for (int j=0; j<numSites; j++)
            {
            int possibleNucs[4];
            alignment->getPossibleNucs(alignment->matrixEntry(i, j), possibleNucs);
            int numPossible = 0, state = -1;
            for (int s=0; s<4; s++)
                {
                if (possibleNucs[s] == 1)
                    {
                    numPossible++;
                    state = s;
                    }
                }
            if (numPossible == 1)
                tipStates[i * numSites + j] = state;
            }
        }
}

double Chunk::lnLikelihood(void) {

#   if defined(NO_DATA)
    return 0.0;
#   endif

    // get the tree and a vector containing all of the transition probabilities
    Tree* t = model->getTree(id);
    std::vector<std::vector<NucleotideSquareMatrix_t> >& tp = transProbs->getTransitionProbabilities();
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    if (t->getRoot()->getDescendants().size() != 3)
        Msg::error("Expecting three descendants of the root");
    if (numGammaCategories > MAX_GAMMA_CATEGORIES)
        Msg::error("Too many gamma categories for the likelihood calculation");
    
    const double* T[MAX_GAMMA_CATEGORIES];
    for (int n=0; n<dpSeq.size(); n++)
        {
        Node* p = dpSeq[n];
        if (p->getIsLeaf() == false)
            {
            std::vector<Node*> desc = p->getDescendants();
            if (desc.size() != 2 && desc.size() != 3)
                Msg::error("Incorrect number of descendants!");
            if (desc.size() == 3 && p != t->getRoot())
                Msg::error("Only expecting a three-way split at the root of the tree");
                
            double* clP = condLikes->getConditionalLikelihoods(0, p->getIndex());
            for (int d=0; d<desc.size(); d++)
                {
                int idx = desc[d]->getIndex();
                const double* clD = condLikes->getConditionalLikelihoods(0, idx);
                for (int k=0; k<numGammaCategories; k++)
                    T[k] = tp[k][idx].data();
                if (desc[d]->getIsLeaf() == true)
                    {
                    const int* states = &tipStates[idx * numSites];
                    if (d == 0)
                        applyTipChild<true>(clP, clD, states, T, numSites, numGammaCategories);
                    else
                        applyTipChild<false>(clP, clD, states, T, numSites, numGammaCategories);
                    }
                else
                    {
                    if (d == 0)
                        applyInternalChild<true>(clP, clD, T, numSites, numGammaCategories);
                    else
                        applyInternalChild<false>(clP, clD, T, numSites, numGammaCategories);
                    }
                }
            }
        }
        
    // calculate average probability at root
    std::vector<double>& f = model->getBaseFrequencies(id);
    double lnL = 0.0;
    double catProb = 1.0 / numGammaCategories;
    double* clRoot = condLikes->getConditionalLikelihoods(0, t->getRoot()->getIndex());
    for (int c=0; c<numSites; c++)
        {
        double prob = 0.0;
        for (int k=0; k<numGammaCategories; k++)
            {
            for (int i=0; i<4; i++)
                prob += clRoot[i] * f[i] * catProb;
            clRoot += 4;
            }
        int numSitesOfPattern = alignment->getNumSitesOfPattern(c);
        lnL += numSitesOfPattern * log(prob);
        }
        
    return lnL;
}

void Chunk::print(void) {

    std::cout << "Data Subset ID: " << id << std::endl;
    alignment->print();
}

void Chunk::printTransitionProbabilities(void) {

    transProbs->print();
}

void Chunk::updateRateMatrix(void) {

    std::vector<double> freqs = model->getBaseFrequencies(id);
    std::vector<double> rates = model->getExchangeabilityRates(id);
    transProbs->updateRateMatrix(rates, freqs);
}

void Chunk::updateTransitionProbabilities(void) {

    Tree* t = model->getTree(id);
    std::vector<double> gammaRates = model->getGammaCategoryRates(id);
    transProbs->setTransitionProbabilities(t, gammaRates);
}
