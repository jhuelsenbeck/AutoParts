#include <fstream>
#include <iostream>
#include "Alignment.hpp"
#include "Branch.hpp"
#include "Chunk.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Node.hpp"
#include "Parameter.hpp"
#include "ParameterBaseFrequencies.hpp"
#include "ParameterExchangabilityRates.hpp"
#include "ParameterGammaShape.hpp"
#include "ParameterTree.hpp"
#include "ParameterTreeLength.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Restaurant.hpp"
#include "StateSets.hpp"
#include "Table.hpp"
#include "Tree.hpp"
#include "UserSettings.hpp"

#undef DEBUG_MODEL

Model::Model(UserSettings* s, int nt, int nss) {

    settings = s;
    stateSets = NULL;
    numSubsets =  nss;
    int numTaxa = nt;
    std::vector<std::string> taxonNames;
    for (int i=0; i<numTaxa; i++)
        taxonNames.push_back("Taxon_" + std::to_string(i+1));
        
    initializeParameters(settings, taxonNames);
}

Model::Model(Alignment* aln, UserSettings* s) {

    std::cout << "   Initializing phylogenetic model" << std::endl;

    settings = s;
    stateSets = NULL;
    alignment = aln;
    initializeStateSets(settings);
    initializeDataChunks(settings);
    initializeParameters(settings, alignment->getTaxonNames());
    initializeProposalProbabilities();
    for (int i=0; i<chunks.size(); i++)
        {
        chunks[i]->updateRateMatrix();
        chunks[i]->updateTransitionProbabilities();
        }
    std::cout << std::endl;
}

Model::Model(Alignment* aln, UserSettings* s, Model* tm) {

    std::cout << "   Initializing phylogenetic model" << std::endl;

    settings = s;
    stateSets = NULL;
    alignment = aln;
    initializeStateSets(settings);
    initializeDataChunks(settings);

    treeRestaurant       = new Restaurant(*tm->treeRestaurant);
    treeRestaurant->setModel(this);
    treeLengthRestaurant = new Restaurant(*tm->treeLengthRestaurant);
    treeLengthRestaurant->setModel(this);
    freqsRestaurant      = new Restaurant(*tm->freqsRestaurant);
    freqsRestaurant->setModel(this);
    ratesRestaurant      = new Restaurant(*tm->ratesRestaurant);
    ratesRestaurant->setModel(this);
    shapeRestaurant      = new Restaurant(*tm->shapeRestaurant);
    shapeRestaurant->setModel(this);
    initializeProposalProbabilities();
    
    treeRestaurant->print();
    treeLengthRestaurant->print();
    freqsRestaurant->print();
    ratesRestaurant->print();
    shapeRestaurant->print();

    for (int i=0; i<chunks.size(); i++)
        {
        chunks[i]->updateRateMatrix();
        chunks[i]->updateTransitionProbabilities();
        }
    std::cout << std::endl;
}

Model::Model(Model& m) {

    settings = m.settings;
    alignment = m.alignment;
    numSubsets = m.numSubsets;
    
    if (m.stateSets != NULL)
        initializeStateSets(settings);

    treeRestaurant       = new Restaurant(*m.treeRestaurant);
    treeRestaurant->setModel(this);
    treeLengthRestaurant = new Restaurant(*m.treeLengthRestaurant);
    treeLengthRestaurant->setModel(this);
    freqsRestaurant      = new Restaurant(*m.freqsRestaurant);
    freqsRestaurant->setModel(this);
    ratesRestaurant      = new Restaurant(*m.ratesRestaurant);
    ratesRestaurant->setModel(this);
    shapeRestaurant      = new Restaurant(*m.shapeRestaurant);
    shapeRestaurant->setModel(this);
    
    initializeProposalProbabilities();
    
    for (int i=0; i<m.chunks.size(); i++)
        {
        Chunk* c = new Chunk(*m.chunks[i]);
        c->setModel(this);
        c->updateRateMatrix();
        c->updateTransitionProbabilities();
        chunks.push_back(c);
        }
}

Model::~Model(void) {

    for (int i=0; i<chunks.size(); i++)
        delete chunks[i];
        
    delete treeRestaurant;
    delete treeLengthRestaurant;
    delete freqsRestaurant;
    delete ratesRestaurant;
    delete shapeRestaurant;
    
    if (stateSets != NULL)
        delete stateSets;
}

char Model::charCode(int x) {

    if (x == 0)
        return 'A';
    else if (x == 1)
        return 'C';
    else if (x == 2)
        return 'G';
    else
        return 'T';
}

int Model::getDegreeTree(void) {

    return treeRestaurant->getNumTables();
}

int Model::getDegreeTreeLength(void) {

    return treeLengthRestaurant->getNumTables();
}

int Model::getDegreeGammaShape(void) {

    return shapeRestaurant->getNumTables();
}

int Model::getDegreeBaseFrequencies(void) {

    return freqsRestaurant->getNumTables();
}

int Model::getDegreeRates(void) {

    return ratesRestaurant->getNumTables();
}

std::vector<double>& Model::getBaseFrequencies(int subsetId) {

    Table* tab = freqsRestaurant->findTableWithPatron(subsetId);
    if (tab == NULL)
        Msg::error("Could not find base frequencies table for subset " + std::to_string(subsetId));
    Parameter* parm = tab->getParameter();
    ParameterBaseFrequencies* freqsParm = dynamic_cast<ParameterBaseFrequencies*>(parm);
    if (freqsParm == NULL)
        Msg::error("Could not find base frequencies for subset " + std::to_string(subsetId));
    return freqsParm->getFreqs();
}

std::vector<double>& Model::getExchangeabilityRates(int subsetId) {

    Table* tab = ratesRestaurant->findTableWithPatron(subsetId);
    if (tab == NULL)
        Msg::error("Could not find base frequencies table for subset " + std::to_string(subsetId));
    Parameter* parm = tab->getParameter();
    ParameterExchangabilityRates* rateParm = dynamic_cast<ParameterExchangabilityRates*>(parm);
    if (rateParm == NULL)
        Msg::error("Could not find exchangability rates for subset " + std::to_string(subsetId));
    return rateParm->getRates();
}

std::vector<double>& Model::getGammaCategoryRates(int subsetId) {

    Table* tab = shapeRestaurant->findTableWithPatron(subsetId);
    Parameter* parm = tab->getParameter();
    ParameterGammaShape* shapeParm = dynamic_cast<ParameterGammaShape*>(parm);
    return shapeParm->getRates();
}

double Model::getGammaShape(int subsetId) {

    Table* tab = shapeRestaurant->findTableWithPatron(subsetId);
    Parameter* parm = tab->getParameter();
    ParameterGammaShape* shapeParm = dynamic_cast<ParameterGammaShape*>(parm);
    return shapeParm->getValue()[0];
}

void Model::getHeader(std::string& v) {

    v.clear();
    
    // tree length
    v += "RGF(Length)";
    v += '\t';
    v += "L(Alpha)";
    v += '\t';
    for (int n=0; n<treeLengthRestaurant->getNumPatrons(); n++)
        {
        std::vector<double> s = treeLengthRestaurant->findTableWithPatron(n)->getParameter()->getValue();
        for (int i=0; i<s.size(); i++)
            {
            v += "L(" + std::to_string(n+1) + "," + std::to_string(i+1) + ")";
            v += '\t';
            }
        }
    
    // frequencies
    std::string pis[4] = { "A", "C", "G", "T" };
    v += "RGF(Pi)";
    v += '\t';
    v += "Pi(Alpha)";
    v += '\t';
    for (int n=0; n<freqsRestaurant->getNumPatrons(); n++)
        {
        std::vector<double> f = freqsRestaurant->findTableWithPatron(n)->getParameter()->getValue();
        for (int i=0; i<f.size(); i++)
            {
            v += "Pi(" + std::to_string(n+1) + "," + pis[i] + ")";
            v += '\t';
            }
        }

    // exchangability rates
    std::string rs[6] = { "AC", "AG", "AT", "CG", "CT", "GT" };
    v += "RGF(R)";
    v += '\t';
    v += "R(Alpha)";
    v += '\t';
    for (int n=0; n<ratesRestaurant->getNumPatrons(); n++)
        {
        std::vector<double> r = ratesRestaurant->findTableWithPatron(n)->getParameter()->getValue();
        for (int i=0; i<r.size(); i++)
            {
            v += "R(" + std::to_string(n+1) + "," + rs[i] + ")";
            v += '\t';
            }
        }

    // gamma shape
    v += "RGF(Alpha)";
    v += '\t';
    v += "Alpha(Alpha)";
    v += '\t';
    for (int n=0; n<shapeRestaurant->getNumPatrons(); n++)
        {
        std::vector<double> s = shapeRestaurant->findTableWithPatron(n)->getParameter()->getValue();
        for (int i=0; i<s.size(); i++)
            {
            v += "Alpha(" + std::to_string(n+1) + "," + std::to_string(i+1) + ")";
            v += '\t';
            }
        }
}

void Model::getParameterValues(std::string& v) {

    v.clear();
    v.reserve(1000);
    
    // tree length
    v += treeLengthRestaurant->rgf();
    v += '\t';
    v += std::to_string(treeLengthRestaurant->getConcentrationParameter());
    v += '\t';
    for (int n=0; n<treeLengthRestaurant->getNumPatrons(); n++)
        {
        std::vector<double> s = treeLengthRestaurant->findTableWithPatron(n)->getParameter()->getValue();
        for (int i=0; i<s.size(); i++)
            {
            v += std::to_string(s[i]);
            v += '\t';
            }
        }
    
    // frequencies
    v += freqsRestaurant->rgf();
    v += '\t';
    v += std::to_string(freqsRestaurant->getConcentrationParameter());
    v += '\t';
    for (int n=0; n<freqsRestaurant->getNumPatrons(); n++)
        {
        std::vector<double> f = freqsRestaurant->findTableWithPatron(n)->getParameter()->getValue();
        for (int i=0; i<f.size(); i++)
            {
            v += std::to_string(f[i]);
            v += '\t';
            }
        }

    // exchangability rates
    v += ratesRestaurant->rgf();
    v += '\t';
    v += std::to_string(ratesRestaurant->getConcentrationParameter());
    v += '\t';
    for (int n=0; n<ratesRestaurant->getNumPatrons(); n++)
        {
        std::vector<double> r = ratesRestaurant->findTableWithPatron(n)->getParameter()->getValue();
        for (int i=0; i<r.size(); i++)
            {
            v += std::to_string(r[i]);
            v += '\t';
            }
        }

    // gamma shape
    v += shapeRestaurant->rgf();
    v += '\t';
    v += std::to_string(shapeRestaurant->getConcentrationParameter());
    v += '\t';
    for (int n=0; n<shapeRestaurant->getNumPatrons(); n++)
        {
        std::vector<double> s = shapeRestaurant->findTableWithPatron(n)->getParameter()->getValue();
        for (int i=0; i<s.size(); i++)
            {
            v += std::to_string(s[i]);
            v += '\t';
            }
        }
}

double Model::getAlphaTree(void) {

    return treeRestaurant->getConcentrationParameter();
}

double Model::getAlphaTreeLength(void) {

    return treeLengthRestaurant->getConcentrationParameter();
}

double Model::getAlphaShape(void) {

    return shapeRestaurant->getConcentrationParameter();
}

double Model::getAlphaFrequencies(void) {

    return freqsRestaurant->getConcentrationParameter();
}

double Model::getAlphaRates(void) {

    return ratesRestaurant->getConcentrationParameter();
}

std::string Model::getPartitionTree(void) {

    return treeRestaurant->rgf();
}

std::string Model::getPartitionTreeLength(void) {

    return treeLengthRestaurant->rgf();
}

std::string Model::getPartitionGammaShape(void) {

    return shapeRestaurant->rgf();
}

std::string Model::getPartitionBaseFrequencies(void) {

    return freqsRestaurant->rgf();
}

std::string Model::getPartitionRates(void) {

    return ratesRestaurant->rgf();
}

std::vector<std::string> Model::getTaxonNames(void) {

    return alignment->getTaxonNames();
}

Tree* Model::getTree(int id) {

    Parameter* p = treeRestaurant->findTableWithPatron(id)->getParameter();
    ParameterTree* pt = dynamic_cast<ParameterTree*>(p);
    if (pt == NULL)
        Msg::error("Failed to get tree parameter");
    Tree* t = pt->getActiveTree();
    return t;
}

double Model::getTreeLength(void) {

    double len = 0.0;
    int totalNumSites = 0;
    for (int i=0; i<numSubsets; i++)
        {
        double iLen = getTreeLength(i);
        int ns = getChunk(i)->getNumSitesForChunk();
        totalNumSites += ns;
        len += iLen * ns;
        }
    len /= totalNumSites;
    return len;
}

double Model::getTreeLength(int subsetId) {

    Table* tab = treeLengthRestaurant->findTableWithPatron(subsetId);
    Parameter* parm = tab->getParameter();
    ParameterTreeLength* treeLengthParm = dynamic_cast<ParameterTreeLength*>(parm);
    return treeLengthParm->getValue()[0];
}

void Model::initializeDataChunks(UserSettings* s) {

    std::set<int> pids = alignment->getPartitionIds();
    numSubsets = (int)pids.size();
            
    for (int i : pids)
        {
        Alignment* a = alignment->dataForPartition(i);
        a->compress();
        Chunk* c = new Chunk(i-1, a, this, s->getNumGammaCategories());
        chunks.push_back(c);
        }
           
#   if defined(DEBUG_MODEL)
    for (Chunk* c : chunks)
        c->print();
#   endif
}

void Model::initializeParameters(UserSettings* s, std::vector<std::string> tn) {

    // set up the tree restaurant with only one table and a fixed seating
    treeRestaurant = new Restaurant(this, s, false, 0.0, numSubsets, new ParameterTree(this, s, tn));
    
    // set up the tree-length restaurant
    double alphaLength = Restaurant::calculateAlphaFromExpectedNumberOfTables(s->getExpectedNumberTreeLengthTables(), numSubsets);
    treeLengthRestaurant = new Restaurant(this, s, true, alphaLength, numSubsets, new ParameterTreeLength(this, s, s->getAlphaT(), s->getBetaT()) );
    
    // set up the base frequencies restaurant
    double alphaPi = Restaurant::calculateAlphaFromExpectedNumberOfTables(s->getExpectedNumberPiTables(), numSubsets);
    freqsRestaurant = new Restaurant(this, s, true, alphaPi, numSubsets, new ParameterBaseFrequencies(this, s));
    
    // set up the exchangeability rates restaurant
    double alphaTheta = Restaurant::calculateAlphaFromExpectedNumberOfTables(s->getExpectedNumberThetaTables(), numSubsets);
    ratesRestaurant = new Restaurant(this, s, true, alphaTheta, numSubsets, new ParameterExchangabilityRates(this, s));

    // set up the gamma shape restaurant
    double alphaAlpha = Restaurant::calculateAlphaFromExpectedNumberOfTables(s->getExpectedNumberAlphaTables(), numSubsets);
    shapeRestaurant = new Restaurant(this, s, true, alphaAlpha, numSubsets, new ParameterGammaShape(this, s, s->getGammaShapeLambda(), s->getNumGammaCategories()));
    
#   if defined(DEBUG_MODEL)
    treeRestaurant->print();
    treeLengthRestaurant->print();
    freqsRestaurant->print();
    ratesRestaurant->print();
    shapeRestaurant->print();
#   endif
}

void Model::initializeProposalProbabilities(void) {

    proposalProbabilities.push_back( std::make_pair(treeRestaurant,       10.0) );
    proposalProbabilities.push_back( std::make_pair(treeLengthRestaurant,  2.0) );
    proposalProbabilities.push_back( std::make_pair(freqsRestaurant,       1.0) );
    proposalProbabilities.push_back( std::make_pair(ratesRestaurant,       1.0) );
    proposalProbabilities.push_back( std::make_pair(shapeRestaurant,       1.0) );
    
    double sum = 0.0;
    for (const auto& r : proposalProbabilities)
        sum += r.second;
    for (auto& r : proposalProbabilities)
        r.second /= sum;
        
#   if defined(DEBUG_MODEL)
    for (const auto& r : proposalProbabilities)
        std::cout << r.first->getParameter()->getName() << " " << r.second << std::endl;
#   endif
}

void Model::initializeStateSets(UserSettings* s) {

	stateSets = new StateSets(alignment, s);
}

// Returns the cached log-likelihood of each data subset, computing it on the first call.
// Restaurant::update keeps these values current, so they are not recomputed afterwards.
std::vector<double>& Model::getChunkLnLs(void) {

    if (chunkLnLValid == false)
        {
        chunkLnL.resize(chunks.size());
        for (int i=0; i<chunks.size(); i++)
            {
            chunks[i]->updateRateMatrix();
            chunks[i]->updateTransitionProbabilities();
            chunkLnL[i] = chunks[i]->lnLikelihood();
            }
        chunkLnLValid = true;
        }
    return chunkLnL;
}

double Model::lnLikelihood(void) {

    double lnL = 0.0;
    for (int i=0; i<chunks.size(); i++)
        lnL += chunks[i]->lnLikelihood();
    return lnL;
}

double Model::lnLikelihood(int idx) {

    return chunks[idx]->lnLikelihood();
}

// Return the prior probability
double Model::lnProbability(void) {

    double lnP = 0.0;
    lnP += treeRestaurant->lnProbability();
    lnP += treeLengthRestaurant->lnProbability();
    lnP += freqsRestaurant->lnProbability();
    lnP += ratesRestaurant->lnProbability();
    lnP += shapeRestaurant->lnProbability();
    return lnP;
}

Alignment* Model::simulate(std::string fn, int ns) {

    // dynamically allocate the matrix
    int nt = getTree(0)->getNumTaxa();
    int** mat = new int*[nt];
    mat[0] = new int[nt * ns * numSubsets];
    for (int i=1; i<nt; i++)
        mat[i] = mat[i-1] + (ns*numSubsets);
    for (int i=0; i<nt; i++)
        for (int j=0; j<ns*numSubsets; j++)
            mat[i][j] = 0;
    
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    // simulate
    for (int s=0; s<numSubsets; s++)
        {
        // get all the parameter for chunk s
        Tree* t = getTree(s);
        double treeLength = getTreeLength(s);
        double alpha = getGammaShape(s);
        std::vector<double>& theta = getExchangeabilityRates(s);
        std::vector<double>& bf = getBaseFrequencies(s);

        std::vector<double> gammaCategories;
        Probability::Gamma::discretization(gammaCategories, alpha, alpha, settings->getNumGammaCategories(), false);
        
        // set up the rate matrix
        double q[4][4];
        for (int i=0, k=0; i<4; i++)
            {
            for (int j=i+1; j<4; j++)
                {
                q[i][j] = theta[k] * bf[j];
                q[j][i] = theta[k] * bf[i];
                k++;
                }
            }
        double averageRate = 0.0;
        for (int i=0; i<4; i++)
            {
            double sum = 0.0;
            for (int j=0; j<4; j++)
                {
                if (i != j)
                    sum += q[i][j];
                }
            q[i][i] = -sum;
            averageRate += bf[i] * sum;
            }
        double scaleFactor = 1.0 / averageRate;
        for (int i=0; i<4; i++)
            for (int j=0; j<4; j++)
                q[i][j] *= scaleFactor;
            
        // simulate the data for this chunk
        std::vector<Node*>& dpSeq = t->getDownPassSequence();
        int* siteInfo = new int[t->getNumNodes()];
        
        for (int c=0; c<ns; c++)
            {
            int randCat = (int)(rng.uniformRv() * gammaCategories.size());
            double r = gammaCategories[randCat];
            for (int n=(int)dpSeq.size()-1; n>=0; n--)
                {
                Node* p = dpSeq[n];
                if (p->getAncestor() == NULL)
                    {
                    // draw sequences form stationary distribution
                    double u = rng.uniformRv();
                    double sum = 0.0;
                    for (int i=0; i<4; i++)
                        {
                        sum += bf[i];
                        if (u < sum)
                            {
                            siteInfo[p->getIndex()] = i;
                            break;
                            }
                        }
                    }
                else
                    {
                    // simulate up the branch
                    Branch* b = t->findBranch(p, p->getAncestor());
                    double v = b->getProportion() * treeLength * r;
                    int curState = siteInfo[p->getAncestor()->getIndex()];
                    double x = 0.0;           
                    while (x < v)
                        {
                        x += Probability::Exponential::rv(&rng, -q[curState][curState]);
                        if (x < v)
                            {
                            double u = rng.uniformRv();
                            double sum = 0.0;
                            for (int i=0; i<4; i++)
                                {
                                if (i != curState)
                                    {
                                    sum += -q[curState][i] / q[curState][curState];
                                    if (u < sum)
                                        {
                                        curState = i;
                                        break;
                                        }
                                    }
                                }
                            }
                        }
                    siteInfo[p->getIndex()] = curState;
                        
                    }
                }
            
            // add the site to the alignment
            for (int i=0; i<nt; i++)
                mat[i][s*ns+c] = siteInfo[i];
            }
            
        delete [] siteInfo;
        }
        
        
        
    // print to a file
    std::ofstream simStrm( fn.c_str(), std::ios::out );
    if (!simStrm)
        Msg::error("Cannot open file \"" + fn + "\"");
    std::vector<std::string> tn = getTree(0)->getTaxonNames();
    int longestLen = 0;
    for (int i=0; i<tn.size(); i++)
        {
        if (tn[i].length() > longestLen)
            longestLen = (int)tn[i].length();
        }
    simStrm << nt << " " << ns*numSubsets << std::endl;
    for (int i=0; i<nt; i++)
        {
        simStrm << tn[i] << " ";
        for (int j=0; j<longestLen-tn[i].length(); j++)
            simStrm << " ";
        for (int j=0; j<ns*numSubsets; j++)
            simStrm << charCode( mat[i][j] );
        simStrm << std::endl;
       }
    for (int i=0; i<numSubsets; i++)
        {
        int firstSite = i * ns + 1;
        int lastSite = (i+1) * ns;
        simStrm << "charset " << "subset_" << i+1 << " = " << firstSite << "-" << lastSite << ";" << std::endl;
        }
    simStrm.close();

    Alignment* aln = new Alignment(fn);

    return aln;
}

// Picks a restaurant to update and returns the likelihood
double Model::update(double power) {

    // pick a parameter to update
    Restaurant* restaurantToUpdate = NULL;
    RandomVariable& rv = RandomVariable::randomVariableInstance();
    double u = rv.uniformRv();
    double sum = 0.0;
    for (auto& r : proposalProbabilities)
        {
        sum += r.second;
        if (u < sum)
            {
            restaurantToUpdate = r.first;
            break;
            }
        }
    if (restaurantToUpdate == NULL)
        Msg::error("Could not find a restaurant to update (" + std::to_string(u) + ")");
        
    // call the update function in the restaurant
    return restaurantToUpdate->update(power);
}
