/**
 * TSPR Decoder
 *
 * Created on: Jun 22, 2016
 *
 */
#ifndef TSPRDECODER_H
#define TSPRDECODER_H

#include <vector>
#include <lemon/list_graph.h>

#include "mygraphlib.h"
#include "TSPRSolver.h"
#include "algoritmos.h"

using namespace lemon;
using namespace std;

class TSPRDecoder {
public:
    TSPRDecoder(TSP_Data_R &tsp,
        const vector<DNode> &terminais, 
        const vector<DNode> &postos,
        const DNode source,
        int delta,
        double avg);
    virtual ~TSPRDecoder();
    
    // Decodes the chromosome into a solution to the TSP
    double decode(const vector< double >& chromosome) const;
    
private:
    TSP_Data_R& instance;
    const vector<DNode> &terminais;    
    const vector<DNode> &postos;    
    const DNode source;
    const int delta;
    const double avg;
};

#endif