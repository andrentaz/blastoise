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

// ********** ALTERE DAQUI PARA BAIXO (MELHOR HEURISTICA) *****************
// Este método deve ter como saída a saída do método heuristica_hv999999() ou heuristica_hg999999(), dependendo de qual for a melhor heuristica.

// ATENÇÃO: Não modifique a assinatura deste método.
bool heuristica_h999999(TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, int maxTime, vector<DNode> &sol, double &lbound) {
    return false;
}
// ********** ALTERE DAQUI PARA CIMA (MELHOR HEURISTICA) *****************