/**
 * TSPRSolver.h
 *
 *  Created on: Jun 22, 2016
 */

#ifndef TSPRSOLVER_H
#define TSPRSOLVER_H

#include <list>
#include <limits>
#include <vector>
#include <algorithm>

class TSPRSolver {
public:
    // The constructor "solves" the problem in O (n log n) by transforming
    // the chromossome into a tour (get a permutation from the chromossome)
    TSPSolver(const TSP_Data_R& instance,
        const vector<DNode> &terminais,
        const vector<DNode> &postos,
        const DNode source,
        const int delta,
        const vector< double >& chromossome);
    virtual ~TSPSolver();
    
    // Returns the tour distance (fitness)
    double getTourDistance() const;
    // Returns the tour
    vector< DNode > getTour() const;

private:
    typedef pair< double, unsigned> ValueKeyPair;
    
    double distance;
    vector< DNode > tour;
    
    bool cook(TSP_Data_R& instance, const vector<DNode>& postos, const int delta);
    bool rerout(TSP_Data_R& instance, list<DNode> list, vector<DNode>& postos, DNode u, DNode v, double delta);
};

#endif