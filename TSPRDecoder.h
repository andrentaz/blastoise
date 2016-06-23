/*
 * TSPR Decoder
 *
 * Created on: Jun 22, 2016
 *
 */
#ifndef TSPRDECODER_H
#define TSPRDECODER_H

#include "TSPRSolver.h"
#include "TSPRInstance.h"
#include <lemon/list_graph.h>

class TSPRDecoder {
public:
    TSPRDecoder(TSP_Data &tsp,
        const vector<DNode> &terminais, 
        const vector<DNode> &postos,
        const DNode source,
        const int delta);
    virtual ~TSPRDecoder();
    
    // Decodes the chromossome into a solution to the TSP
    double decode(const vector< double >& chromossome) const;
    vector<DNode> getSolution();
    
private:
    const TSP_Data& instance;
    const vector<DNode> &terminais;    
    const vector<DNode> &postos;    
    const DNode source;
    const int delta;
};

#endif