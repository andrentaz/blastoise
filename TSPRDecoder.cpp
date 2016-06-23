/*
 * TSPRDecoder.cpp
 *
 *  Created on: 22 Jun, 2016
 */
 
#include "TSPRDecoder.h"

TSPRDecoder::TSPRDecoder(TSP_Data &_instance,
        const vector<DNode> &_terminais, 
        const vector<DNode> &_postos,
        const DNode _source,
        const int _delta) : 
        instance(_instance), terminais(_terminais), postos(_postos), source(_source), delta(_delta) { }

TSPRDecoder::~TSPRDecoder() { }

double TSPRDecoder::decode(const vector< double >& chromossome) const {
    // 1) Create a tour out of the chromossome
    TSPRSolver solver(instance, terminais, postos, source, delta, chromossome);
    
    // 2) Extract the fitness
    const unsigned fitness = solver.getTourDistance();
    
    // 3) Return
    return double(fitness);
}

vector<DNode> getSolution() {
    
}