/**
 * TSPRSolver.cpp
 *
 *  Created on: Jun 22, 2016
 */

#include "TSPRSolver.h"

using namespace std;
using namespace lemon;

TSPRSolver::TSPRSolver(TSP_Data_R& instance, const vector<DNode>& terminais, 
    const vector<DNode>& postos, const DNode source, int delta, double avg, const vector< double >& chromosome) :
    distance(0) {
    
    vector< ValueKeyPair > rawtour;
    
    // Assumes that instance.getNumNodes() == chromosome.size()
    
    // 1) Obtain a permutation out of the chromosome -- this is the raw tour:
    for (unsigned i=0; i<chromosome.size()-1; ++i) { 
        rawtour.push_back(ValueKeyPair(chromosome[i], i)); 
    }
    
    // sort 'rank', which will produce a permutation of [n] stored in ValueKeyPair::second:
    sort(rawtour.begin(), rawtour.end());
    
    // 2) Insert the source to create a feasible solution from the rawtour
    this->tour.push_back(source);
    for (unsigned i=0; i<rawtour.size(); ++i) {
        this->tour.push_back(terminais[rawtour[i].second]);
    }
    this->tour.push_back(source);
    
    // 3) Transforms the rawtour in a feasible one
    if (cook(instance, postos, delta, avg*chromosome.back())) {

        // 4) Compute the distance of the tour given by the permutation and cooking
        for (unsigned i=1; i<this->tour.size(); ++i) {
            // Get the nodes
            const DNode u = this->tour[i-1];
            const DNode v = this->tour[i];
            
            // Calculate the distance
            distance = distance + instance.AdjMatD.Cost(u,v);
        }
    } else {
        // Problems doing the solution
        distance = DBL_MAX;
    }
    
}

/* Get the raw solution and turns into a feasible tour */
bool TSPRSolver::cook(TSP_Data_R& instance, const vector<DNode>& postos, const int delta, double avg) {
    int fuelUsed = 0.0;
    bool retVal = true;
    vector<DNode> cooked;
    double refuel;
    
    // Starting from the second position in the list
    DNode u = this->tour[0];
    DNode v;
    
    // initialize the new vector
    cooked.push_back(u);
    for (unsigned i=1; i<this->tour.size(); ++i) {
        // assign v
        v = this->tour[i];
        
        // check the distance from the previous node
        if (fuelUsed + instance.AdjMatD.Cost(u,v) + avg < delta) {
            fuelUsed += instance.AdjMatD.Cost(u,v);
            cooked.push_back(v);
        } else {
            // limit would be passed, rerout
            refuel = rerout(instance, cooked, postos, u, v, double(delta-fuelUsed));
            if (refuel > 0) {
                // reset the fuel used
                fuelUsed = refuel;
                cooked.push_back(v);
            } else {
                // some problem happened
                retVal = false;
                break;
            }
        }
        // assign u
        u = v;
    }
    
    // refresh the tour
    if (retVal) {
        this->tour = cooked;
    }
    
    return retVal;
}

TSPRSolver::~TSPRSolver() { }

/* Insert a fuel station to respect the fuel limit */
double TSPRSolver::rerout(TSP_Data_R& instance, vector<DNode>& ctour, const vector<DNode>& postos, 
                        DNode u, DNode v, double delta) {
    DNode x;
    vector<DNode> pp;
    double retVal = -1.0;
    
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
    if (min != DBL_MAX) {
        // insert in the tour
        ctour.push_back(x);
        // set the return        
        retVal = instance.AdjMatD.Cost(x,v);
    }
        
    return retVal;
}

double TSPRSolver::getTourDistance() const { return distance; }

vector< DNode > TSPRSolver::getTour() const { return tour; }