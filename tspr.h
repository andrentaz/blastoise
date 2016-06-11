// MC658 - Projeto e Análise de Algoritmos III
// ---- LABORATÓRIO TSPR ----
// Prof: Flavio Keidi Miyazawa
// PED: Rafael Arakaki

#ifndef LAB4_TSPR_H
#define LAB4_TSPR_H

// ATENÇÃO: Não altere este arquivo.

// Estrutura de dados que armazena informações de uma instância do problema TSP with refueling.
class TSP_Data_R {
public:
    TSP_Data_R(ListDigraph &graph,
               DNodeStringMap &nodename,
               DNodePosMap &posicaox,
               DNodePosMap &posy,
               ArcValueMap &eweight);
    ListDigraph &g;
    int NNodes,NEdges;
    DNodeStringMap &vname;
    ArcStringMap ename;
    DNodeColorMap vcolor;
    ArcColorMap ecolor;
    ArcValueMap &weight;
    DNodePosMap &posx;
    DNodePosMap &posy;
    AdjacencyMatrixDirected AdjMatD;
};

// Usage information
void showUsage ();

// Calcula o custo de uma solução
double solutionCost (const TSP_Data_R &tsp, const vector<DNode> &path_sol);

// Escreve a saída do programa
void writeOutputFile (string outputfile, string graphname, vector<DNode> path, int delta, TSP_Data_R &tsp, DNode source, double ub, double lb, double elapsed_time, int max_time);

// Funções para formar caminhos a partir de conjuntos de arestas (usados no B&B).
list<DNode> path_s (DNode u, ListDigraph::NodeMap< list<DNode> > &adj);
vector<DNode> path_search (TSP_Data_R &tsp, const DNode source, const vector<Arc>& setEdges);

// Visualiza um grafo com a solução do TSP-R (um caminho)
void ViewTspRCircuit(TSP_Data_R &tsp, const vector<DNode>& bestRoute, string graphname, double delta, const vector<DNode> &postos);

// Global variables reference
extern bool verbose; //0=not set, 1=print solution after optimization

#endif //LAB4_TSPR_H
