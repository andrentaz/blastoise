// MC658 - Projeto e Análise de Algoritmos III
// ---- LABORATÓRIO TSPR ----
// Prof: Flavio Keidi Miyazawa
// PED: Rafael Arakaki

#include <iostream>
#include <float.h>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "algoritmos.h"

using namespace lemon;
using namespace std;

// ********** ALTERE DAQUI PARA BAIXO (HEURISTICA GENETICA) *****************
// Otimiza o problema TSP-R através de uma heurística baseada em BRKGA

// ATENÇÃO: Não modifique a assinatura deste método.
bool heuristica_hg999999(TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos,
                         const DNode source,
                         int delta, int maxTime, vector<DNode> &sol, double &lbound) {
    return false;
}
// ********** ALTERE DAQUI PARA CIMA (HEURISTICA GENETICA) *****************