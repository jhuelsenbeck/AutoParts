#include <cmath>
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

// Conditional likelihoods are rescaled to avoid underflow on large trees. The likelihood is
// first computed without rescaling. Every factor in the calculation is at most one, so if
// each site's probability at the root is at least MIN_UNSCALED_SITE_PROB, nothing that
// contributes to it can have underflowed and that result is used. Otherwise it is computed
// again with rescaling: when the largest value for a site at a node drops below
// 2^-SCALE_EXPONENT, the site's values at that node are multiplied by 2^SCALE_EXPONENT,
// which is exact.
#define SCALE_EXPONENT 256
#define MIN_UNSCALED_SITE_PROB 1e-290


Chunk::Chunk(const Chunk& c) {

    alignment = c.alignment;
    model = c.model;
    id = c.id;
    numSites = c.numSites;
    numGammaCategories = c.numGammaCategories;

    condLikes = new ConditionalLikelihoods(alignment, numGammaCategories);
    transProbs = new TransitionProbabilities(model, this, 2*alignment->getNumTaxa()-1, numGammaCategories);
    initializeTipStates();
    ownCondLikePtrs.resize(condLikes->getNumNodes());
    for (int i=0; i<condLikes->getNumNodes(); i++)
        ownCondLikePtrs[i] = condLikes->getConditionalLikelihoods(0, i);
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
    ownCondLikePtrs.resize(condLikes->getNumNodes());
    for (int i=0; i<condLikes->getNumNodes(); i++)
        ownCondLikePtrs[i] = condLikes->getConditionalLikelihoods(0, i);
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

    // scale up any site whose conditional likelihoods at this node are close to underflowing
    inline void rescaleNode(double* clP, int* scaleCounts, int numSites, int numCats) {

        const double threshold = std::ldexp(1.0, -SCALE_EXPONENT);
        const double factor    = std::ldexp(1.0,  SCALE_EXPONENT);
        const int n = 4 * numCats;
        for (int c=0; c<numSites; c++)
            {
            double maxVal = 0.0;
            for (int i=0; i<n; i++)
                maxVal = (clP[i] > maxVal) ? clP[i] : maxVal;
            while (maxVal < threshold && maxVal > 0.0)
                {
                for (int i=0; i<n; i++)
                    clP[i] *= factor;
                maxVal *= factor;
                scaleCounts[c]++;
                }
            clP += n;
            }
    }

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

    Tree* t = model->getTree(id);
    std::vector<std::vector<NucleotideSquareMatrix_t> >& tp = transProbs->getTransitionProbabilities();
    std::vector<double>& f = model->getBaseFrequencies(id);
    bool tooSmall = false;
    double lnL = computeLnLikelihood(false, tooSmall, t, tp, ownCondLikePtrs.data(), f, scaleCounts);
    if (tooSmall == true)
        lnL = computeLnLikelihood(true, tooSmall, t, tp, ownCondLikePtrs.data(), f, scaleCounts);
    return lnL;
}

// The current parameter values for this subset (those of the tables it is seated at).
SubsetParameters Chunk::currentParameters(void) {

    SubsetParameters p;
    p.tree                 = model->getTree(id);
    p.treeLength           = model->getTreeLength(id);
    p.baseFrequencies      = &model->getBaseFrequencies(id);
    p.exchangeabilityRates = &model->getExchangeabilityRates(id);
    p.gammaCategoryRates   = &model->getGammaCategoryRates(id);
    return p;
}

// The likelihood of this subset under the parameter values p, calculated in the
// workspace ws; this subset's own buffers are only read. If newEigens is false, the
// subset's current eigensystem is used (p must then have its current base frequencies
// and exchangeability rates). Different threads can call this at once for the same
// subset, each with its own workspace.
double Chunk::lnLikelihood(const SubsetParameters& p, bool newEigens, LikelihoodWorkspace& ws) {

#   if defined(NO_DATA)
    return 0.0;
#   endif

    EigenSystem* e = &transProbs->getEigens();
    if (newEigens == true)
        {
        ws.eigens.updateReversibleRateMatrix(*p.exchangeabilityRates, *p.baseFrequencies);
        e = &ws.eigens;
        }
        
    int numNodes = condLikes->getNumNodes();
    int numTaxa = alignment->getNumTaxa();
    if (ws.transitionProbabilities.size() != numGammaCategories || ws.transitionProbabilities[0].size() < numNodes)
        {
        ws.transitionProbabilities.resize(numGammaCategories);
        for (int k=0; k<numGammaCategories; k++)
            ws.transitionProbabilities[k].resize(numNodes);
        }
    TransitionProbabilities::calculate(*e, p.tree, p.treeLength, *p.gammaCategoryRates, ws.transitionProbabilities);
    
    // tips use this subset's own (fixed) conditional likelihoods; internal nodes use the workspace
    size_t nodeSize = (size_t)numSites * 4 * numGammaCategories;
    if (ws.condLikes.size() < numNodes * nodeSize)
        ws.condLikes.resize(numNodes * nodeSize);
    ws.condLikePtrs.resize(numNodes);
    for (int i=0; i<numNodes; i++)
        ws.condLikePtrs[i] = (i < numTaxa) ? ownCondLikePtrs[i] : ws.condLikes.data() + i * nodeSize;

    bool tooSmall = false;
    double lnL = computeLnLikelihood(false, tooSmall, p.tree, ws.transitionProbabilities, ws.condLikePtrs.data(), *p.baseFrequencies, ws.scaleCounts);
    if (tooSmall == true)
        lnL = computeLnLikelihood(true, tooSmall, p.tree, ws.transitionProbabilities, ws.condLikePtrs.data(), *p.baseFrequencies, ws.scaleCounts);
    return lnL;
}

double Chunk::computeLnLikelihood(bool rescale, bool& tooSmall, Tree* t,
                                  std::vector<std::vector<NucleotideSquareMatrix_t> >& tp,
                                  double* const* clPtrs, const std::vector<double>& f,
                                  std::vector<int>& scaleCountsForSites) {

    tooSmall = false;
    std::vector<Node*>& dpSeq = t->getDownPassSequence();
    if (t->getRoot()->getDescendants().size() != 3)
        Msg::error("Expecting three descendants of the root");
    if (numGammaCategories > MAX_GAMMA_CATEGORIES)
        Msg::error("Too many gamma categories for the likelihood calculation");
    
    const double* T[MAX_GAMMA_CATEGORIES];
    if (rescale == true)
        scaleCountsForSites.assign(numSites, 0);
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
                
            double* clP = clPtrs[p->getIndex()];
            for (int d=0; d<desc.size(); d++)
                {
                int idx = desc[d]->getIndex();
                const double* clD = clPtrs[idx];
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
            if (rescale == true)
                rescaleNode(clP, scaleCountsForSites.data(), numSites, numGammaCategories);
            }
        }
        
    // calculate average probability at root
    double lnL = 0.0;
    double catProb = 1.0 / numGammaCategories;
    const double* clRoot = clPtrs[t->getRoot()->getIndex()];
    for (int c=0; c<numSites; c++)
        {
        double prob = 0.0;
        for (int k=0; k<numGammaCategories; k++)
            {
            for (int i=0; i<4; i++)
                prob += clRoot[i] * f[i] * catProb;
            clRoot += 4;
            }
        if (rescale == false && !(prob >= MIN_UNSCALED_SITE_PROB))
            {
            tooSmall = true;    // possible underflow: the caller recomputes with rescaling
            return 0.0;
            }
        int numSitesOfPattern = alignment->getNumSitesOfPattern(c);
        double lnProb = log(prob);
        if (rescale == true && scaleCountsForSites[c] > 0)
            lnProb -= scaleCountsForSites[c] * SCALE_EXPONENT * std::log(2.0);
        lnL += numSitesOfPattern * lnProb;
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
