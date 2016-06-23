/*
 * TSPRSolver.cpp
 *
 *  Created on: Jun 22, 2016
 */

#include "TSPRSolver.h"

using namespace std;
using namespace lemon;

TSPRSolver::TSPRSolver(const TSP_Data& instance, const vector<DNode>& terminais, 
    const vector<DNode>& postos, const DNode source, const int delta, const vector< double >& chromossome) :
    distance(0) {
    
    vector< ValueKeyPair > rawtour;
    
    // Assumes that instance.getNumNodes() == chromossome.size()
    
    // 1) Obtain a permutation out of the chromossome -- this is the raw tour:
    for (unsigned i=0; i<chromossome.size(); ++i) { rawtour[i] = ValueKeyPair(chromossome[i], i); }
    
    // sort 'rank', which will produce a permutation of [n] stored in ValueKeyPair::second:
    sort(rawtour.begin(), rawtour.end());
    
    // 2) Insert the source to create a feasible solution from the rawtour
    this->tour.push_back(source);
    for (unsigned i=0; i<rawtour.size(); ++i) {
        this->tour.push_back(terminais[rawtour[i].second]);
    }
    this->tour.push_back(source);
    
    // 3) Transforms the rawtour in a feasible one
    cook(this->tour, postos, delta);
    
    // 4) Compute the distance of the rawtour given by the permutation
    for (unsigned i=1; i<this->tour.size(); ++i) {
        // Get the nodes
        const DNode u = this->tour[i-1];
        const DNode v = this->tour[i];
        
        // Calculate the distance
        distance += tsp.AdjMatD.Cost(u,v);
    }
    
}

/* Get the raw solution and turns into a feasible tour */
bool TSPRSolver::cook(TSP_Data& instance, vector<DNode>& ftour, const vector<DNode>& postos, const int delta) {
    int fuelUsed = 0.0;
    bool retVal = true;
    list<DNode> listTour;
    
    // Simple copy the tour into a list
    for (auto u : ftour ) {
        listTour.push_back(u);
    }
    
    // Starting from the second position in the list
    for (list<DNode>::iterator it=listTour.begin()+1; it != listTour.end(); ++it) {
        // check the distance from the previous node
        if (fuelUsed + instance.AdjMatD.Cost(*(it-1), *it) < delta) {
            fuelUsed += instance.AdjMatD.Cost(*(it-1), *it);
        } else {
            // limit would be passed, rerout
            if (rerout(instance, listTour, postos, u, v, (delta-fuelUsed))) {
              // reset the fuel used
              fuelUsed = 0.0;
            } else {
                // some problem happened
                retVal = false;
                break;
            }
        }
    }
    
    return retVal;    
}

TSPRSolver::~TSPRSolver() { }

/* Insert a fuel station to respect the fuel limit */
void rerout(TSP_Data& instance, list<DNode> tourlist, vector<DNode>& postos, DNode u, DNode v, double delta) {
    vector<DNode> pp;
    bool retVal = false;
    DNode x = NULL;
    
    // Take the feasible fuel station from u
    for (auto p : postos) {
        if (instance.AdjMatD.Cost(u,p) <= delta) {
            pp.push_back(p);
        }
    }
    
    // Calculate the minimum triangle to rerout
    double min = DBL_MAX;
    for (auto p : pp) {
        if (instance.AdjMatD.Cost(u,p) + instance.AdjMatD.Cost(p,v) <= min ) {
            x = p;
            min = instance.AdjMatD.Cost(u,p) + instance.AdjMatD.Cost(p,v);
        }
    }
    
    // there are no reachable fuel station
    if (x!=NULL) {        
        // insert in the list
        for (list<DNode>::iterator it=tourlist.begin(); it!=tourlist.end(); ++it) {
            // Find the node
            if (*it == v) {
                tourlist.insert(it, x);
            }
        }
        
        // set the return val
        retVal = true;
    }
        
    return retVal;
}

unsigned TSPRSolver::getTourDistance() const { return distance; }

vector< DNode > TSPRSolver::getTour() const { return tour; }