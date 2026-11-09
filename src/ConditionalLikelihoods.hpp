#ifndef ConditionalLikelihoods_hpp
#define ConditionalLikelihoods_hpp

class Alignment;



class ConditionalLikelihoods {

    public:
                                    ConditionalLikelihoods(void) = delete;
                                    ConditionalLikelihoods(Alignment* aln, int nc);
                                   ~ConditionalLikelihoods(void);
        double*                     getConditionalLikelihoods(int space, int nodeIdx) { return cls[space][nodeIdx]; }
        int                         getNumNodes(void) { return numNodes; }
        void                        print(void);
        
    private:
        std::vector<std::string>    checkConditionalLikelihoods(Alignment* aln, double tolerance);
        int                         numTaxa;
        int                         numSites;
        int                         numNodes;
        int                         numStates;
        int                         numRateCategories;
        double*                     clsRaw[2];
        double**                    cls[2];
};

#endif
