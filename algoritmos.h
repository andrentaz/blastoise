// MC658 - Projeto e Análise de Algoritmos III
// ---- LABORATÓRIO TSPR ----
// Prof: Flavio Keidi Miyazawa
// PED: Rafael Arakaki

#ifndef LAB4_ALGORITMOS_H
#define LAB4_ALGORITMOS_H

#include "tspr.h"

// ATENÇÃO: Não altere este arquivo.

// Rotinas de algoritmos exatos e heurísticos p/ TSPR
bool brach_and_bound999999 (TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, int maxTime, vector<DNode> &sol, double &lbound);
bool heuristica_hv999999   (TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, int maxTime, vector<DNode> &sol, double &lbound);
bool heuristica_hg999999   (TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, int maxTime, vector<DNode> &sol, double &lbound);
bool heuristica_h999999    (TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, int maxTime, vector<DNode> &sol, double &lbound);

#endif //LAB4_ALGORITMOS_H
