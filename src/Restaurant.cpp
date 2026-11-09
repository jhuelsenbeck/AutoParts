#include <cmath>
#include <iomanip>
#include <iostream>
#include "Chunk.hpp"
#include "Model.hpp"
#include "Msg.hpp"
#include "Parameter.hpp"
#include "ParameterBaseFrequencies.hpp"
#include "ParameterExchangabilityRates.hpp"
#include "ParameterGammaShape.hpp"
#include "ParameterTree.hpp"
#include "ParameterTreeLength.hpp"
#include "Probability.hpp"
#include "RandomVariable.hpp"
#include "Restaurant.hpp"
#include "Table.hpp"
#include "TableFactory.hpp"
#include "ThreadPool.hpp"
#include "UpdateInfo.hpp"
#include "UserSettings.hpp"



Restaurant::Restaurant(const Restaurant& r) {

    modelPtr = r.modelPtr;
    settingsPtr = r.settingsPtr;
    alpha = r.alpha;
    numPatrons = r.numPatrons;
    isSeatingRv = r.isSeatingRv;
    gammaAlpha = r.gammaAlpha;
    gammaBeta = r.gammaBeta;
        
    parameter = copyParameter(r.parameter);
    if (parameter == NULL)
        Msg::error("Could not find parameter to copy in restaurant copy constructor");
        
    TableFactory& tf = TableFactory::tableFactory();
    for (Table* t : r.tables)
        {
        Table* newT = tf.getTable();
        newT->setRestaurant(this);
        std::set<int>& tPatrons = t->getPatrons();
        for (int i : tPatrons)
            newT->addPatron(i);
        newT->addParameter( copyParameter(t->getParameter()) );
        tables.insert(newT);
        }
}

Restaurant::Restaurant(Model* mp, UserSettings* s, bool sf, double a, int np, Parameter* parm) {

    isSeatingRv = sf;
    modelPtr    = mp;
    settingsPtr = s;
    numPatrons  = np;
    parameter   = parm;
    
    if (isSeatingRv == false)
        {
        // seating is not a random variable, so we seat all patrons at one table
        std::cout << "   * Parameter " << parameter->getName() << " subset arrangement is fixed (all subsets grouped together)" << std::endl;
        alpha = 0.0;
        Table* t = addTable();
        for (int i=0; i<numPatrons; i++)
            t->addPatron(i);
        std::cout << "   *    Fixed RGF = " << rgf() << std::endl;
        }
    else
        {
        // seating is a random variable following a DPP
        RandomVariable& rng = RandomVariable::randomVariableInstance();
        
        if (settingsPtr->getIsConcentrationParameterFixed() == true)
            {
            std::cout << "   * Parameter " << parameter->getName() << " subset arrangement is a random variable following a DPP with fixed alpha" << std::endl;
            alpha = a; // a is passed into the function based on E(#Tables)
            std::cout << "   *    Fixed Alpha = " << alpha << std::endl;
            }
        else
            {
            std::cout << "   * Parameter " << parameter->getName() << " subset arrangement is a random variable following a DPP with a gamma hyperprior on alpha" << std::endl;
            double m = settingsPtr->getPriorConcMean();
            if(settingsPtr->getPriorMeanTables() != 0){
                double eT = settingsPtr->getPriorMeanTables();
                m = calculateAlphaFromExpectedNumberOfTables(eT, np);
            }
            double v = settingsPtr->getPriorConcVariance();
            gammaBeta  = m / v;
            gammaAlpha = gammaBeta * m;

            // If the prior is for a small number of clusters, there is a chance you will initialize it to 0.0 (or something else unstable), throwing an error
            int initialization_attempt = 0;
            do {
                alpha = Probability::Gamma::rv(&rng, gammaAlpha, gammaBeta);
                initialization_attempt++;
            }
            while(alpha < 1e-6 && initialization_attempt <= 100);
            if(initialization_attempt > 100){
                Msg::error("Failed to initialize a valid alpha under the prior Gamma(" + std::to_string(gammaAlpha) + "," + std::to_string(gammaBeta) + ") for " + parameter->getName());
            }

            std::cout << "   *    Initial Alpha = " << alpha << std::endl;
            }
        
        // This can help with instantiation issues we are seeing with extreme alphas
        Table* initTable = addTable();
        initTable->addPatron(0);

        for (int i=1; i<numPatrons; i++)
            {
            double newTableProb = alpha / (alpha + i);
            double u1 = rng.uniformRv();
            if (u1 < newTableProb)
                {
                // add a new table
                Table* t = addTable();
                t->addPatron(i);
                }
            else
                {
                // add the patron to one of the existing tables
                double u2 = rng.uniformRv();
                double sum = 0.0;
                for (Table* t : tables)
                    {
                    sum += (double)t->getNumPatrons() / i;
                    if (u2 < sum)
                        {
                        t->addPatron(i);
                        break;
                        }
                    }
                }
            }
        std::cout << "   *    Initial RGF = " << rgf() << std::endl;

        }
}

Restaurant::~Restaurant(void) {

    delete parameter;
}

Table* Restaurant::addTable(void) {

    Table* t = TableFactory::tableFactory().getTable();
    Parameter* newParm = parameter->newRandomInstance();
    newParm->setTable(t);
    t->addParameter( newParm );
    t->setRestaurant(this);
    tables.insert(t);
    return t;
}

Table* Restaurant::addAuxiliaryTable(void) {

    Table* t = TableFactory::tableFactory().getTable();
    Parameter* newParm = parameter->newRandomInstance();
    newParm->setTable(t);
    t->addParameter( newParm );
    t->setRestaurant(this);
    return t;
}

double Restaurant::calculateAlphaFromExpectedNumberOfTables(double expT, int np) {

    if (expT > np)
        Msg::error("The expected number of tables cannot be larger than the number of patrons (" + std::to_string(np) + "<" + std::to_string(expT) + ")");
     if (expT < 1.0)
        Msg::error("The expected number of tables cannot be less than one");
       
    double a = 0.000001;
    double ea = expectedNumberOfTables(a, np);
    bool goUp;
    if (ea < expT)
        goUp = true;
    else
        goUp = false;
    double increment = 0.1;
    while ( fabs(ea - expT) > 0.000001 )
        {
        if (ea < expT && goUp == true)
            {
            a += increment;
            }
        else if (ea > expT && goUp == false)
            {
            a -= increment;
            }
        else if (ea < expT && goUp == false)
            {
            increment /= 2.0;
            goUp = true;
            a += increment;
            }
        else
            {
            increment /= 2.0;
            goUp = false;
            a -= increment;
            }
        ea = expectedNumberOfTables(a, np);
        }
    //cout << ea << " <-> " << expT << " " << "alpha=" << a << endl;
    return a;
}

Table* Restaurant::chooseTable(std::map<Table*,double,TableLess>& lnProbs) {

    RandomVariable& rv = RandomVariable::randomVariableInstance();
    double u = rv.uniformRv();
    double sum = 0.0;
    for (std::map<Table*,double,TableLess>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        {
        sum += it->second;
        if (u < sum)
            return it->first;
        }
    // the probabilities can sum to slightly less than one because of rounding
    if (lnProbs.empty() == false)
        return lnProbs.rbegin()->first;
    return NULL;
}

Parameter* Restaurant::copyParameter(Parameter* parmToCopy) {

    Parameter* parm = NULL;
    if (dynamic_cast<ParameterTree*>(parmToCopy) != NULL)
        parm = new ParameterTree(*dynamic_cast<ParameterTree*>(parmToCopy));
    else if (dynamic_cast<ParameterTreeLength*>(parmToCopy) != NULL)
        parm = new ParameterTreeLength(*dynamic_cast<ParameterTreeLength*>(parmToCopy));
    else if (dynamic_cast<ParameterBaseFrequencies*>(parmToCopy) != NULL)
        parm = new ParameterBaseFrequencies(*dynamic_cast<ParameterBaseFrequencies*>(parmToCopy));
    else if (dynamic_cast<ParameterExchangabilityRates*>(parmToCopy) != NULL)
        parm = new ParameterExchangabilityRates(*dynamic_cast<ParameterExchangabilityRates*>(parmToCopy));
    else if (dynamic_cast<ParameterGammaShape*>(parmToCopy) != NULL)
        parm = new ParameterGammaShape(*dynamic_cast<ParameterGammaShape*>(parmToCopy));
    return parm;
}

double Restaurant::expectedNumberOfTables(double a, int np) {

    double expectedNum = 0.0;
    for (int i=1; i<=np; i++)
        expectedNum += ( 1.0 / (i - 1.0 + a) );
    expectedNum *= a;
    return expectedNum;
}

Table* Restaurant::findTableWithPatron(int idx) {

    for (Table* t : tables)
        {
        if (t->hasPatron(idx) == true)
            return t;
        }
    return NULL;
}

double Restaurant::lnProbability(void) {

    double lnP = 0.0;
    for (Table* t : tables)
        {
        Parameter* parm = t->getParameter();
        lnP += parm->lnProbability();
        }
    return lnP;
}

void Restaurant::normalize(std::map<Table*,double,TableLess>& lnProbs) {

    // find maximum value
    bool first = true;
    double maxVal = 0.0;
    for (std::map<Table*,double,TableLess>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        {
        if (first == true)
            {
            maxVal = it->second;
            first = false;
            }
        else
            {
            if (it->second > maxVal)
                maxVal = it->second;
            }
        }
        
    // rescale
    for (std::map<Table*,double,TableLess>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        it->second -= maxVal;
        
    // exponentiate
    double sum = 0.0;
    for (std::map<Table*,double,TableLess>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        {
        double lnP = it->second;
        double p = exp(lnP);
        it->second = p;
        sum += p;
        }
        
    // rescale to probs
    for (std::map<Table*,double,TableLess>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        {
        double p = it->second / sum;
        it->second = p;
        }

#   if 0
    for (std::map<Table*,double,TableLess>::iterator it = lnProbs.begin(); it != lnProbs.end(); it++)
        std::cout << it->first << " " << it->second << std::endl;
#   endif
}

void Restaurant::print(void) {

    std::cout << "Restaurant " << parameter->getName() << " (" << rgf() << ")" << std::endl;
    int i = 0;
    for (Table* t : tables)
        {
        Parameter* parm = t->getParameter();
        std::cout << "   Table " << i << ": ";
        std::cout << parm->getValuesAsString(5);
        i++;

        std::cout << " [ ";
        std::set<int>& patrons = t->getPatrons();
        for (int p : patrons)
            std::cout << p << " ";
        std::cout << "]" << std::endl;
        }
}

std::string Restaurant::rgf(void) {

    std::vector<int> rgf(numPatrons);
    for (int i=0; i<numPatrons; i++)
        rgf[i] = 0;
    int idx = 0;
    for (int i=0; i<numPatrons; i++)
        {
        if (rgf[i] == 0)
            {
            idx++;
            Table* t = findTableWithPatron(i);
            if (t == NULL)
                Msg::error("Could not find patron " + std::to_string(i) + " in the restaurant (" + parameter->getName() + ")");
            std::set<int> tPatrons = t->getPatrons();
            for (int p : tPatrons)
                {
                rgf[p] = idx;
                }
            }
        }
    
    std::string rgfStr = "";
    for (int i=0; i<numPatrons; i++)
        {
        if (i != 0)
            rgfStr += ",";
        rgfStr += std::to_string(rgf[i]);
        }
    return rgfStr;
}

void Restaurant::removeTable(Table* tab) {

    TableFactory& tf = TableFactory::tableFactory();
    tf.returnToPool(tab);
    tables.erase(tab);
}

double Restaurant::sampleAlpha(int k, int n, double oldAlpha, double a, double b) {

    RandomVariable& rng = RandomVariable::randomVariableInstance();

	/* Step 1: Draw a Beta(oldAlpha+1, n) distribution to get eta */
	std::vector<double> z(2);
	std::vector<double> f(2);
	z[0] = oldAlpha + 1.0;
	z[1] = (double)n;
    Probability::Dirichlet::rv(&rng, z, f);
	double eta = f[0];
	
	/* Step 2: Draw a new value for alpha, based on k and eta */
	double u = rng.uniformRv();
	double x = ( a + (double)k - 1.0 ) / ( (double)n * (b - log(eta)) );
	double newAlpha;
	if ( (u / (1.0 - u)) < x)
		newAlpha = Probability::Gamma::rv(&rng, a + k, b - log(eta));
	else
		newAlpha = Probability::Gamma::rv(&rng, a + k - 1.0, b - log(eta));
		
	return newAlpha;
}

double Restaurant::update(double power) {

    // The update of a restaurant does two things. First, we visit all of the
    // current tables in the restaurant, and update the parameter on each. We
    // then reseat patrons at tables according to Algorithm 8 of Neal (2000).
        
    RandomVariable& rng = RandomVariable::randomVariableInstance();
    
    // the current likelihood of each patron (data subset), cached in the model and kept
    // up to date below, so it does not need to be recomputed at the start of each update
    std::vector<double>& curLnL = modelPtr->getChunkLnLs();
    std::vector<double> newLnL(numPatrons);
    
    // Likelihoods are calculated in parallel. Random numbers are drawn, and results are
    // combined, only outside of the parallel loops and in a fixed order, so a run gives
    // the same results whatever the number of threads.
    ThreadPool& pool = ThreadPool::threadPool();
    modelPtr->prepareWorkspaces(pool.getNumThreads());
        
    // 1. update the parameters on each table. The tables have no patrons in common, so
    // the proposals for all of the tables are made first (in table order, as they draw
    // random numbers), the likelihoods of all of the patrons are then calculated in one
    // parallel loop, and finally each proposal is accepted or rejected, in table order.
    // Done this way, tables with few patrons do not leave most of the threads idle.
    UpdateInfo& updateInfo = UpdateInfo::updateInfo();
    std::vector<Table*> tableList(tables.begin(), tables.end());
    int numTables = (int)tableList.size();
    std::vector<double> tableLnPriorRatio(numTables), tableLnProposalRatio(numTables);
    std::vector<UpdateType> tableUpdateType(numTables);
    std::vector<int> patronTable(numPatrons, -1);
    std::vector<int> patronList;
    patronList.reserve(numPatrons);
    for (int k=0; k<numTables; k++)
        {
        // change the parameter
        Parameter* parm = tableList[k]->getParameter();
        updateInfo.setAttemptedUpdate(UpdateType::NO_MOVE);
        int numSitesAtTable = 0;     // the proposal's step size is scaled by this (see UpdateInfo)
        for (int i : tableList[k]->getPatrons())
            numSitesAtTable += modelPtr->getChunk(i)->getNumObservedSites();
        updateInfo.setNumSites(numSitesAtTable);
        double lnPriorRatio = -parm->lnProbability();
        tableLnProposalRatio[k] = parm->update();
        lnPriorRatio += parm->lnProbability();
        tableLnPriorRatio[k] = lnPriorRatio;
        tableUpdateType[k] = updateInfo.getAttemptedUpdate();
        for (int i : tableList[k]->getPatrons())
            {
            patronTable[i] = k;
            patronList.push_back(i);
            }
        }

    updateInfo.setNumSites(0);

    // calculate the likelihoods (each patron has its own buffers, so the patrons are done in parallel)
    pool.parallelFor((int)patronList.size(), [&](int j, int) {
        int i = patronList[j];
        Chunk* chunk = modelPtr->getChunk(i);
        if (tableList[patronTable[i]]->getParameter()->getUpdateModifiesEigens() == true)
            {
            chunk->flipActiveEigens();
            chunk->updateRateMatrix();
            }
        chunk->flipActiveTransitionProbabilities();
        chunk->updateTransitionProbabilities();
        newLnL[i] = chunk->lnLikelihood();
        });

    for (int k=0; k<numTables; k++)
        {
        Table* tab = tableList[k];
        Parameter* parm = tab->getParameter();
        std::set<int>& tablePatrons = tab->getPatrons();
        double lnLikelihoodRatio = 0.0;
        for (int i : tablePatrons)
            lnLikelihoodRatio += (newLnL[i] - curLnL[i]);
            
        // accept/reject
        double lnR = (lnLikelihoodRatio * power) + tableLnPriorRatio[k] + tableLnProposalRatio[k]; // only the likelihood is heated
        bool accept = false;
        if (log(rng.uniformRv()) < lnR)
            accept = true;
            
        // update
        updateInfo.setAttemptedUpdate(tableUpdateType[k]);
        if (accept == true)
            {
            // proposed update is accepted
            parm->accept();
            for (int i : tablePatrons)
                curLnL[i] = newLnL[i];
            updateInfo.accept();
            }
        else
            {
            // proposed update is rejected
            parm->reject();
            for (int i : tablePatrons)
                {
                Chunk* chunk = modelPtr->getChunk(i);
                if (parm->getUpdateModifiesEigens() == true)
                    chunk->flipActiveEigens();
                chunk->flipActiveTransitionProbabilities();
                }
            updateInfo.reject();
            }
        }

    // 2. loop over all of the patrons, reseating each (alg. 8)
    if (isSeatingRv == true)
        {
        TableSet auxiliaryTables;
        int numAuxiliaryTables = settingsPtr->getNumAuxiliaryTables();   // m in Neal (2000)
        for (int n=0; n<numPatrons; n++)
            {
            // get the chunk for the patron, and its current parameter values (which must be
            // read before the patron leaves its table)
            Chunk* chunk = modelPtr->getChunk(n);
            SubsetParameters patronParameters = chunk->currentParameters();
            
            // remove the patron from its current table
            Table* tab = findTableWithPatron(n);
            if (tab == NULL)
                Msg::error("Could not find table with patron " + std::to_string(n));
            tab->removePatron(n);
            
            // if the patron was alone at its table, Algorithm 8 keeps that table (and its
            // parameter value) as the first auxiliary table, so only m-1 fresh tables are
            // drawn; if it is not chosen, it is removed below along with the other empty tables
            int numFreshTables = numAuxiliaryTables;
            if (tab->getNumPatrons() == 0)
                {
                auxiliaryTables.insert(tab);
                numFreshTables--;
                }

            // make the remaining auxiliary tables, drawing the parameter values from the prior
            // distribution; the auxiliary tables are added to the tables set
            for (int i=0; i<numFreshTables; i++)
                auxiliaryTables.insert( addTable() );
                            
            // calculate the likelihood when the patron is seated at each of the tables
            // (including the auxiliary tables), in parallel, each thread in its own workspace
            std::vector<Table*> candidates(tables.begin(), tables.end());
            std::vector<double> candidateLnL(candidates.size());
            bool newEigens = parameter->getUpdateModifiesEigens();
            pool.parallelFor((int)candidates.size(), [&](int j, int w) {
                SubsetParameters p = patronParameters;
                candidates[j]->getParameter()->setSubsetParameters(p);
                candidateLnL[j] = chunk->lnLikelihood(p, newEigens, modelPtr->getWorkspace(w));
                });

            std::map<Table*,double,TableLess> tableLikes;
            std::map<Table*,double,TableLess> tableLnL;     // unheated log-likelihood of the patron at each table
            for (int j=0; j<candidates.size(); j++)
                {
                Table* t = candidates[j];
                double rawLnL = candidateLnL[j];
                double lnL = rawLnL * power;
                TableSet::iterator it = auxiliaryTables.find(t);
                double lnP = 0.0;
                if (it == auxiliaryTables.end())
                    lnP = log( (double)t->getNumPatrons() );  // the patron is not seated, so this is n_{-i,c}
                else
                    lnP = log(alpha/numAuxiliaryTables);
                tableLikes.insert( std::make_pair(t, lnL+lnP) );
                tableLnL.insert( std::make_pair(t, rawLnL) );
                }

            // choose a table in proportion to the likelihood to reseat
            normalize(tableLikes);
            tab = chooseTable(tableLikes);
            if (tab == NULL)
                Msg::error("Could not find a new table for reseating patron " + std::to_string(n));
            tab->addPatron(n);
            
            // sort out the likelihood for the patron
            if (tab->getParameter()->getUpdateModifiesEigens() == true)
                chunk->updateRateMatrix();
            chunk->updateTransitionProbabilities();
            curLnL[n] = tableLnL[tab];
            
            // remove any unused tables
            std::vector<Table*> emptyTables;
            for (Table* t : tables)
                {
                if (t->getNumPatrons() == 0)
                    emptyTables.push_back(t);
                }
            for (int i=0; i<emptyTables.size(); i++)
                removeTable(emptyTables[i]);
            auxiliaryTables.clear();
            }
        }
 
    // 3. the total likelihood, from the cached per-patron values
    double lnL = 0.0;
    for (int i=0; i<numPatrons; i++)
        lnL += curLnL[i];

#   if defined(CHECK_LNL_CACHE)
    for (int i=0; i<numPatrons; i++)
        {
        Chunk* chunk = modelPtr->getChunk(i);
        chunk->updateRateMatrix();
        chunk->updateTransitionProbabilities();
        double fresh = chunk->lnLikelihood();
        if (fresh != curLnL[i])
            Msg::error("Cached likelihood for subset " + std::to_string(i) + " is stale (" + std::to_string(curLnL[i]) + " vs " + std::to_string(fresh) + ")");
        }
#   endif
        
    // 4. update the concentration parameter (Gibb's sampler)
	if ( settingsPtr->getIsConcentrationParameterFixed() == false && isSeatingRv == true )
		alpha = sampleAlpha( getNumTables(), numPatrons, alpha, gammaAlpha, gammaBeta );

    return lnL;
}
