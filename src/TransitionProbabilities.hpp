#ifndef TransitionProbabilities_hpp
#define TransitionProbabilities_hpp

#include "EigenSystem.hpp"
#include <vector>
class Chunk;
class Model;
class Tree;



class TransitionProbabilities {

    public:
                                                                TransitionProbabilities(void) = delete;
                                                                TransitionProbabilities(Model* m, Chunk* c, int nb, int nc);
                                                               ~TransitionProbabilities(void);
        static void                                             calculate(EigenSystem& e, Tree* t, double treeLength, const std::vector<double>& gammaCategoryRates, std::vector<std::vector<NucleotideSquareMatrix_t> >& tp);
        void                                                    flipActiveEigens(void);
        EigenSystem&                                            getEigens(void) { return *eigens; }
        void                                                    flipActiveProbs(void);
        std::vector<std::vector<NucleotideSquareMatrix_t> >&    getTransitionProbabilities(void) { return transitionProbabilities[activeProbs]; }
        void                                                    print(void);
        void                                                    setTransitionProbabilities(Tree* t, std::vector<double> gammaCategoryRates);
        void                                                    updateRateMatrix(std::vector<double>& rates, std::vector<double>& f);
    
    private:
        int                                                     numBranches;
        int                                                     numGammaCategories;
        int                                                     activeProbs;
        Chunk*                                                  myDataChunk;
        EigenSystem*                                            eigens;
        Model*                                                  model;
        std::vector<std::vector<NucleotideSquareMatrix_t> >     transitionProbabilities[2];
};

#endif
