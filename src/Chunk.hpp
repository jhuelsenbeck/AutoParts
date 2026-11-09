#ifndef Chunk_hpp
#define Chunk_hpp

#include <vector>

class Alignment;
class ConditionalLikelihoods;
class Model;
class TransitionProbabilities;


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
        double                      lnLikelihood(void);
        void                        printTransitionProbabilities(void);
        void                        setId(int i) { id = i; }
        void                        setModel(Model* m) {model = m;}
        void                        updateRateMatrix(void);
        void                        updateTransitionProbabilities(void);
        
    private:
        void                        initializeTipStates(void);
        int                         id;
        int                         numSites;
        int                         numGammaCategories;
        Alignment*                  alignment;
        Model*                      model;
        ConditionalLikelihoods*     condLikes;
        TransitionProbabilities*    transProbs;
        std::vector<int>            tipStates;          // [taxon * numSites + site]: nucleotide 0-3, or -1 if ambiguous
};

#endif
