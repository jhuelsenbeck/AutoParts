#include <iostream>
#include <limits>
#include "Hungarian.hpp"
#include "Msg.hpp"
#include "Partition.hpp"
#include "Subset.hpp"



Partition::Partition(std::vector<int> x) {

    // we assume the vector, x, is RGF format with the first index being 1
    rgf = x;
    numElements = (int)rgf.size();
    
    int numSubsets = 0;
    for (int i=0; i<rgf.size(); i++)
        {
        if (rgf[i] > numSubsets)
            numSubsets = rgf[i];
        }
    
    for (int i=0; i<numSubsets; i++)
        {
        Subset* ss = new Subset(numElements);
        for (int j=0; j<numElements; j++)
            {
            if (rgf[j] == i+1)
                ss->addElement(j);
            }
        subsets.insert(ss);
        }
}

Partition::Partition(const Partition& p) {

    rgf = p.rgf;
    numElements = p.numElements;
    for (std::set<Subset*>::iterator it = p.subsets.begin(); it != p.subsets.end(); it++)
        subsets.insert( new Subset(*(*it)) );
}

Partition::Partition(std::map<Partition*,int>& partList, double& score) {

    // check that all of the partitions in the list have the same size
    int partitionSize = partList.begin()->first->getNumElements();
    std::map<Partition*,int>::iterator mostFrequentPartitionIt;
    int cnt = 0;
    for (std::map<Partition*,int>::iterator p=partList.begin(); p != partList.end(); p++)
        {
        if ( p->first->getNumElements() != partitionSize )
            Msg::error("Partition list has an inconsistent number of elements");
        if (p->second > cnt)
            {
            cnt = p->second;
            mostFrequentPartitionIt = p;
            }
        }
    //mostFrequentPartitionIt = partList.begin();
    
    // set this partition to be equal to the most frequent partition
    rgf = mostFrequentPartitionIt->first->rgf;
    numElements = mostFrequentPartitionIt->first->numElements;
    for (std::set<Subset*>::iterator it = mostFrequentPartitionIt->first->subsets.begin(); it != mostFrequentPartitionIt->first->subsets.end(); it++)
        subsets.insert( new Subset(*(*it)) );

    /* heuristic search looking for better partitions (i.e., partitions that
       minimize the squared distance to all of the partitions in the
       list of partitions, prtList) */
    double bestDist = ssDistance(partList);
    bool foundBetter = false;
    int passNum = 1;
    //std::cout << "   Searching for mean partition" << std::endl << std::endl;
    do
        {
        foundBetter = false;
        for (int i=0; i<numElements; i++)
            {
            std::map<Subset*,double> distanceScores;
            distanceScores.insert( std::make_pair(findSubsetWithElement(i),bestDist) );
            
            // remove the element from its current subset
            Subset* origSs = removeElement(i);
            
            // place the element in each of the other subsets
            for (Subset* ss : subsets)
                {
                if (ss != origSs)
                    {
                    ss->addElement(i);
                    double d = ssDistance(partList);
                    distanceScores.insert( std::make_pair(ss,d) );
                    if ( removeElement(i) != NULL )
                        Msg::error("We shouldn't be removing the element from a set with with no other elements");
                    }
                }
            
            // place the element in a new subset by itself
            Subset* newSs = NULL;
            if (origSs == NULL)
                {
                newSs = new Subset(numElements);
                newSs->addElement(i);
                subsets.insert( newSs );
                double d = ssDistance(partList);
                distanceScores.insert( std::make_pair(newSs,d) );
                if ( removeElement(i) != newSs)
                    Msg::error("The returned subset should have been the same as the one that was just inserted");
                }
                
            // find the best neighbor
            double minD = std::numeric_limits<double>::max();
            Subset* bestSs = NULL;
            for (std::map<Subset*,double>::iterator it = distanceScores.begin(); it != distanceScores.end(); it++)
                {
                if (it->second < minD)
                    {
                    minD = it->second;
                    bestSs = it->first;
                    }
                }
                
            // adjust the partition
            bestSs->addElement(i);
            if (origSs == NULL)
                {
                if (bestSs == newSs)
                    subsets.insert(bestSs);
                else
                    delete newSs;
                }
            else
                {
                if (bestSs == origSs)
                    subsets.insert(bestSs);
                else
                    delete origSs;
                }
            
            // did we find a better partition?
            if ( minD < bestDist )
                {
                bestDist = minD;
                foundBetter = true;
                }
                
            // these two lines are for debugging
            setRgf();
            //std::cout << passNum << " " << i << " -- " << bestDist << " -- " << getRgfString() << std::endl;
            }
        passNum++;
        } while (foundBetter == true);
    score = bestDist;
    
    // put partition into RGF form
    setRgf();
}

Partition::~Partition(void) {

    deleteSubsets();
}

Partition& Partition::operator=(Partition& p) {

    if (this != &p)
        {
        numElements = p.numElements;
        rgf = p.rgf;
        deleteSubsets();
        for (std::set<Subset*>::iterator it = p.subsets.begin(); it != p.subsets.end(); it++)
            subsets.insert( new Subset(*(*it)) );
        }
    return *this;
}

void Partition::deleteSubsets(void) {

    for (Subset* ss : subsets)
        delete ss;
    subsets.clear();
}

int Partition::distance(Partition* p) {

    Partition* p1 = this;
    Partition* p2 = p;
    
    if (p1 == p2)
        return 0;
    
    if ( p1->getNumElements() != p2->getNumElements() )
        return -1;
    
    // partition 1 should have a larger degree than partition 2
    Partition* part1 = p1;
    Partition* part2 = p2;;
    if ( p1->degree() < p2->degree() )
        {
        part1 = p2;
        part2 = p1;
        }
    
    // set the number of rows/columns of the cost matrix
    int ncols = part1->degree();
    int nrows = part2->degree();
    
    // allcoate the cost matrix
    int** r = new int*[nrows];
    r[0] = new int[nrows * ncols];
    for (int i=1; i<nrows; i++)
        r[i] = r[i-1] + ncols;
    for (int i=0; i<nrows; i++)
        for (int j=0; j<ncols; j++)
            r[i][j] = 0;
    
    // initialize the cost matrix
    int i = 0;
    for (Subset* ss2 : part2->subsets)
        {
        std::vector<bool> bf2 = ss2->getBitRep();
        int j = 0;
        for (Subset* ss1 : part1->subsets)
            {
            std::vector<bool> bf1 = ss1->getBitRep();
            for (int k=0; k<numElements; k++)
                {
                if ( (bf1[k] & bf2[k]) == true )
                    r[i][j]++;
                }
            j++;
            }
        i++;
        }
    
#   if 0
    for (int i=0; i<nrows; i++)
        {
        for (int j=0; j<ncols; j++)
            std::cout << r[i][j] << ",";
        std::cout << std::endl;
        }
    std::cout << "done" << std::endl;
    getchar();
#   endif

    // calculate assignment cost
    Hungarian* h = new Hungarian(r, nrows, ncols);
    int d = numElements - h->getAssignmentCost();
    
    // free memory
    delete [] r[0];
    delete [] r;
    delete h;

    return d;
}

Subset* Partition::findSubsetWithElement(int x) {

    for (Subset* ss : subsets)
        {
        if ( ss->isElementPartOfSet(x) == true )
            return ss;
        }
    return NULL;
}

int Partition::flip(int x) {

    if (x == 0)
        return 1;
    return 0;
}

std::string Partition::getRgfString(void) {

    std::string s = "(";
    for (int i=0; i<rgf.size(); i++)
        s += std::to_string(rgf[i]) + ",";
    s += ")";
    return s;
}

std::string Partition::getRgfString(void) const {

    std::string s = "(";
    for (int i=0; i<rgf.size(); i++)
        s += std::to_string(rgf[i]) + ",";
    s += ")";
    return s;
}

void Partition::print(void) {

    for (Subset* ss : subsets)
        {
        std::vector<bool> part = ss->getBitRep();
        std::cout << "{ ";
        for (int i=0; i<numElements; i++)
            {
            if (part[i] == true)
                std::cout << i+1 << " ";
            }
        std::cout << "} ";
        }
    std::cout << std::endl;
}

void Partition::printRgf(void) const {

    for (int i=0; i<rgf.size(); i++)
        std::cout << rgf[i] << " ";
    std::cout << std::endl;
}

void Partition::printRgf(void) {

    for (int i=0; i<rgf.size(); i++)
        std::cout << rgf[i] << " ";
    std::cout << std::endl;
}

Subset* Partition::removeElement(int x) {

    Subset* ss = findSubsetWithElement(x);
    ss->removeElement(x);
    if (ss->numAssigned() == 0)
        {
        subsets.erase(ss);
        return ss;
        }
    return NULL;
}

void Partition::removeEmptySubsets(void) {

    std::vector<Subset*> subsetsToRemove;
    for (Subset* ss : subsets)
        {
        if (ss->numAssigned() == 0)
            subsetsToRemove.push_back(ss);
        }
    for (int i=0; i<subsetsToRemove.size(); i++)
        subsets.erase( subsetsToRemove[i] );
}

void Partition::setRgf(void) {

    rgf.resize(numElements);
    for (int i=0; i<numElements; i++)
        rgf[i] = -1;
    int idx = 1;
    for (int i=0; i<numElements; i++)
        {
        if (rgf[i] == -1)
            {
            Subset* ss = findSubsetWithElement(i);
            std::vector<bool> v = ss->getBitRep();
            for (int j=i; j<numElements; j++)
                {
                if (v[j] == true)
                    rgf[j] = idx;
                }
            idx++;
            }
        }
}

double Partition::ssDistance(std::map<Partition*,int>& partList) {

    // find the partition with the largest degree (including this one, which can have a
    // larger degree than any in the list)
    int numberPartitions = 0;
    int largestDegree = this->degree();
    for (std::map<Partition*,int>::iterator p=partList.begin(); p != partList.end(); p++)
        {
        if ( p->first->degree() > largestDegree )
            largestDegree = p->first->degree();
        numberPartitions += p->second;
        }
    
    // allocate a matrix large enough to hold any cost matrix in the list
    int** r = new int*[largestDegree];
    r[0] = new int[largestDegree * largestDegree];
    for (int i=1; i<largestDegree; i++)
        r[i] = r[i-1] + largestDegree;
    
    // calculate the average distance from this partition to all of the paritions in the list
    double distanceSum = 0.0;
    for (std::map<Partition*,int>::iterator p=partList.begin(); p != partList.end(); p++)
        {
        // get pointers to the two partitions
        Partition* part1 = this;
        Partition* part2 = p->first;
        if ( this->degree() < p->first->degree() )
            {
            part1 = p->first;
            part2 = this;
            }

        // set the number of rows/columns of the cost matrix
        int ncols = part1->degree();
        int nrows = part2->degree();
        
        // zero-out the cost matrix
        for (int i=0; i<nrows; i++)
            for (int j=0; j<ncols; j++)
                r[i][j] = 0;

        // initialize the cost matrix
        int i = 0;
        for (Subset* ss2 : part2->subsets)
            {
            std::vector<bool> bf2 = ss2->getBitRep();
            int j = 0;
            for (Subset* ss1 : part1->subsets)
                {
                std::vector<bool> bf1 = ss1->getBitRep();
                for (int k=0; k<numElements; k++)
                    {
                    if ( (bf1[k] & bf2[k]) == true )
                        r[i][j]++;
                    }
                j++;
                }
            i++;
            }

        // calculate assignment cost
        Hungarian h(r, nrows, ncols);
        int d = numElements - h.getAssignmentCost();
        
        // add to the sum
        distanceSum += (d * d) * p->second;
        
        //delete h;
        }
    
    // free memory
    delete [] r[0];
    delete [] r;
    
    // return the average distance
    return distanceSum / numberPartitions;
}

