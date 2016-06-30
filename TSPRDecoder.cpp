/**
 * TSPRDecoder.cpp
 *
 *  Created on: 22 Jun, 2016
 */
 
#include "TSPRDecoder.h"

TSPRDecoder::TSPRDecoder(TSP_Data_R &_instance,
        const vector<DNode> &_terminais, 
        const vector<DNode> &_postos,
        const DNode _source,
        int _delta,
        double _avg) : 
        instance(_instance), terminais(_terminais), postos(_postos), source(_source), delta(_delta), avg(_avg) { }

TSPRDecoder::~TSPRDecoder() { }

double TSPRDecoder::decode(const vector< double >& chromosome) const {
    // 1) Create a tour out of the chromosome
    TSPRSolver solver(instance, terminais, postos, source, delta, avg, chromosome);
    
    // 2) Extract the fitness
    const double fitness = solver.getTourDistance();
    
    // 3) Return
    return fitness;
}