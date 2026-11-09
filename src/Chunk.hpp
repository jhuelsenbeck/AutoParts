#ifndef Chunk_hpp
#define Chunk_hpp

#include <vector>
#include "EigenSystem.hpp"
#include "SubsetParameters.hpp"

class Alignment;
class ConditionalLikelihoods;
class Model;
class TransitionProbabilities;
class Tree;

// Scratch space for calculating a subset's likelihood under parameter values other than
// its current ones, without touching the subset's own buffers. Each thread has its own.
struct LikelihoodWorkspace {

    EigenSystem                                             eigens;
    std::vector<std::vector<NucleotideSquareMatrix_t> >     transitionProbabilities;    // [gamma category][node index]
    std::vector<double>                                     condLikes;                  // internal nodes
    std::vector<double*>                                    condLikePtrs;               // [node index]
    std::vector<int>                                        scaleCounts;
};


class Chunk {

    public:
                                    Chunk(void) = delete;
                                    Chunk(const Chunk& c);
                                    Chunk(int myId, Alignment* a, Model* m, int nc);
                                   ~Chunk(void);
        void                        flipActiveTransitionProbabilities(void);
        void                        flipActiveEigens(void);
        int                         getId(void) { return id; }
        int                         getNumSitesForChunk(void);
        void                        print(void);
        SubsetParameters            currentParameters(void);
        double                      lnLikelihood(void);
        double                      lnLikelihood(const SubsetParameters& p, bool newEigens, LikelihoodWorkspace& ws);
        void                        printTransitionProbabilities(void);
        void                        setId(int i) { id = i; }
        void                        setModel(Model* m) {model = m;}
        void                        updateRateMatrix(void);
        void                        updateTransitionProbabilities(void);
        
    private:
        double                      computeLnLikelihood(bool rescale, bool& tooSmall, Tree* t,
                                                        std::vector<std::vector<NucleotideSquareMatrix_t> >& tp,
                                                        double* const* clPtrs, const std::vector<double>& f,
                                                        std::vector<int>& scaleCountsForSites);
        void                        initializeTipStates(void);
        int                         id;
        int                         numSites;
        int                         numGammaCategories;
        Alignment*                  alignment;
        Model*                      model;
        ConditionalLikelihoods*     condLikes;
        TransitionProbabilities*    transProbs;
        std::vector<int>            tipStates;          // [taxon * numSites + site]: nucleotide 0-3, or -1 if ambiguous
        std::vector<int>            scaleCounts;
        std::vector<double*>        ownCondLikePtrs;    // [node index]: this subset's own conditional likelihoods        // [site]: times the site's conditional likelihoods were scaled up by 2^SCALE_EXPONENT
};

#endif
