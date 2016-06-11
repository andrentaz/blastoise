// MC658 - Projeto e Análise de Algoritmos III
// ---- LABORATÓRIO TSPR ----
// Prof: Flavio Keidi Miyazawa
// PED: Rafael Arakaki

#include <iostream>
#include <float.h>
#include <gurobi_c++.h>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "tspr.h"
#include "algoritmos.h"

using namespace lemon;
using namespace std;

// ATENÇÃO: Não altere este arquivo.

// Global variables definition
bool verbose = 0; //0=not set, 1=print solution after optimization

// TSP_Data_R initializer
TSP_Data_R::TSP_Data_R(ListDigraph &dgraph,
                       DNodeStringMap &nodename,
                       DNodePosMap &posicaox,
                       DNodePosMap &posicaoy,
                       ArcValueMap &eweight):
        g(dgraph),
        vname(nodename),
        ename(dgraph),
        vcolor(dgraph),
        ecolor(dgraph),
        weight(eweight),
        posx(posicaox),
        posy(posicaoy),
        AdjMatD(dgraph,eweight,MY_INF)
{
    NNodes=countNodes(this->g);
    NEdges=countArcs(this->g);
}

// Usage information
void showUsage (){
    cout << "Usage: ./tspr <modo_operacao> (um dentre: -e branch&bound, -hv heuristica_vizinhanca, -hg heuristica_genetico, -h melhor_heuristica) -t <tempo_max_segundos> {-v: printa solucao ao final} "
    << "-d <valor_delta> -r <nome_vertice_origem> -f <nome_grafo_entrada> -o <nome_arquivo_saida>" << endl;
}

// Main program
int main(int argc, char** argv)
{
    int exec = 0; // 0=not set, 1=B&B/C, 2=heuristic (neighbour), 3=heuristic (genetic), 4=best heuristic.
    int maxTime = 0; //0=not set
    int delta = 0; //0=not set, otherwise=vehicle fuel (capacity)
    string sourceVertex_name; // source vertex name
    string inputFile_name; // Input graph file
    string outputFile_name; // Output sol file

    // Reading program arguments
    for(int i = 1; i < argc; i++){
        const string arg(argv[i]);
        string next;
        if((i+1) < argc)
            next = string(argv[i+1]);
        else
            next = string("");

        if( exec != 0 && (arg.find("-e") == 0 || arg.find("-h") == 0 || arg.find("-hv") == 0 || arg.find("-hg") == 0 ) ){
            cout << "Erro ao ler parametro \"" << arg << "\", somente pode haver um parametro de modo de execucao." << endl;
            showUsage();
            exit(1);
        }
        else if( arg.find("-e") == 0 ){
            exec = 1;
        }
        else if( arg.find("-hv") == 0 ){
            exec = 2;
        }
        else if( arg.find("-hg") == 0 ){
            exec = 3;
        }
        else if( arg.find("-h") == 0 ){
            exec = 4;
        }
        else if( arg.find("-v") == 0 ){
            verbose=1;
        }
        else if( arg.find("-t") == 0 && next.size() > 0){
            maxTime = atoi(next.c_str()); i++; continue;
        }
        else if( arg.find("-d") == 0 && next.size() > 0){
            delta = atoi(next.c_str()); i++; continue;
        }
        else if( arg.find("-r") == 0 && next.size() > 0){
            sourceVertex_name = next; i++; continue;
        }
        else if( arg.find("-f") == 0 && next.size() > 0){
            inputFile_name = next; i++; continue;
        }
        else if( arg.find("-o") == 0 && next.size() > 0){
            outputFile_name = next; i++; continue;
        }
        else{
            cout << "Parametro invalido: \"" << arg << "\"" << " (ou faltando argumento)" << endl;
            showUsage();
            exit(1);
        }
    }

    // Required parameters
    if( exec == 0 ){
        cout << "Nenhum modo de execucao selecionado dentre: -e -h -hg ou -hv" << endl;
        showUsage(); exit(1);
    }
    if (sourceVertex_name.size() < 1 || inputFile_name.size() < 1 || outputFile_name.size() < 1){
        cout << "Argumentos obrigatorios faltando: " << ((sourceVertex_name.size() < 1)? "nome do vertice de origem, ":"") <<
                ((inputFile_name.size() < 1)? "nome do arq de grafo, ":"") << ((outputFile_name.size() < 1)? "nome do arq de saida.":"") << endl;
        showUsage(); exit(1);
    }
    if( maxTime == 0 || delta == 0){
        cout << "Argumentos obrigatorios faltando: " << ((maxTime == 0)?"-t <tempo_limite_segundos> ":"")
        << ((delta == 0)?"-d <valor_delta>":"") << endl;
        showUsage(); exit(1);
    }

    ListDigraph g;
    ArcValueMap weight(g);
    DNodeStringMap vname(g);
    DNodePosMap   posx(g),posy(g);
    DNode source;
    DNodeBoolMap is_terminal(g);

    // Leitura do grafo
    ReadListDigraph(inputFile_name, g, vname, weight, posx, posy, is_terminal, true);
    TSP_Data_R tsp(g,vname,posx,posy,weight);

    // Encontra nó origem
    bool found = false;
    for(DNodeIt v(g); v!=INVALID; ++v){
        if(vname[v].compare(sourceVertex_name) == 0){
            source = v; found=true;
        }
    }
    if(!found){
        cout << "Source vertex not found: " << sourceVertex_name << endl;
        exit(1);
    }

    // Separa nós em terminais e postos
    vector<DNode> terminais;
    vector<DNode> postos; int p = 0; int t = 0;
    for(DNodeIt v(g); v!=INVALID; ++v) {
        if(is_terminal[v] == true){
            terminais.push_back(v); ++t;
        }
        else{
            postos.push_back(v); ++p;
        }
    }
    //cout << "#Debug# " << endl << "Lidos " << p << " postos e " << t << " terminais (total nos: " << countNodes(g) << ")" << endl;

    // Otimiza o problema através de alguma rotina
    vector<DNode> path_sol;
    double lbound = 0, ubound = DBL_MAX;
    double elapsed_time;

    // Os algoritmos devolvem: (1) solução em forma de caminho de nodes (path_sol) e (2) limitante inferior, se houver (lbound).
    // São calculados aqui: (1) tempo de processamento (elapsed_time), (2) custo da solução obtida pelos algoritmos (ubound).
    clock_t before = clock();
    bool foundSolution = false;
    switch(exec){
        case 1:
            foundSolution = brach_and_bound999999(tsp, terminais, postos, source, delta, maxTime, path_sol, lbound);
            break;
        case 2:
            foundSolution = heuristica_hv999999(tsp, terminais, postos, source, delta, maxTime, path_sol, lbound);
            break;
        case 3:
            foundSolution = heuristica_hg999999(tsp, terminais, postos, source, delta, maxTime, path_sol, lbound);
            break;
        case 4:
            foundSolution = heuristica_h999999(tsp, terminais, postos, source, delta, maxTime, path_sol, lbound);
            break;
    }
    clock_t after = clock();
    elapsed_time = (double) (after-before) / CLOCKS_PER_SEC;

    if( foundSolution )
        ubound = solutionCost(tsp, path_sol);

    writeOutputFile(outputFile_name, inputFile_name, path_sol, delta, tsp, source, ubound, lbound, elapsed_time, maxTime);

    // Imprime a solução
    if(verbose){
        if( foundSolution ){
            cout << "#####" << endl << "Solucao (caminho):" << endl;
            for(auto v : path_sol){
                cout << tsp.vname[v] << (is_terminal[v]?" ":"* ");
            }
            cout << endl << "*: posto" << endl << "#####" << endl;
            cout << "Solution cost = " << solutionCost(tsp, path_sol) << endl;
            ViewTspRCircuit(tsp, path_sol, inputFile_name, delta, postos);
            cout << endl;
        }
        else{
            cout << "Nenhuma solucao viavel encontrada." << endl;
        }
    }

    return 0;
}

// Calcula o custo de uma solução TSP-R (caminho). (Obs: não verifica se a solução é realmente viável).
double solutionCost (const TSP_Data_R &tsp, const vector<DNode> &path_sol){
    double sum = 0.0;
    for (auto it = path_sol.begin(); it != path_sol.end(); ++it) {
        DNode u = *it; ++it;
        if(it != path_sol.end()){
            DNode v = *it;
            bool found=false;
            for(OutArcIt e(tsp.g,u); e != INVALID; ++e){
                if( tsp.g.source(e) == u && tsp.g.target(e) == v ){
                    // aresta (u,v) do caminho
                    sum += tsp.weight[e];
                    found=true; break;
                }
            }
            if(!found){
                cout << "solutionCost() ERRO: Aresta de vertice " << tsp.vname[u] << " para " << tsp.vname[v] <<" nao encontrada." << endl;
                exit(1);
            }
            --it;
        }
        else
            break;
    }
    return sum;
}

// Search function to find paths
list<DNode> path_s (DNode u, ListDigraph::NodeMap< list<DNode> > &adj){
    list<DNode> myList;
    while(adj[u].size() > 0){
        DNode v = *adj[u].begin();
        adj[u].remove(v); // remove arco (u,v) apenas nessa direção
        //adj[v].remove(u);
        list<DNode> path = path_s (v, adj);
        if( *path.rbegin() == u ){
            myList.splice(myList.begin(), path); // insere path ao começo da lista (forma ciclo!)
        }
        else{
            myList.splice(myList.end(), path); // insere path ao fim da lista (somente um v de \in adj[u] será este caso).
        }
    }
    myList.insert(myList.begin(), u); // myList: u {<cycle #1> <cycle #2> ... <cycle #n>} {<non-cycle path>}
    return myList;
}

// Get vector<DNode> path in respect to a set of edges.
vector<DNode> path_search (TSP_Data_R &tsp, const DNode source, const vector<Arc>& setEdges){
    ListDigraph::NodeMap< list<DNode> > adj(tsp.g);
    for(auto e : setEdges){
        adj[tsp.g.source(e)].push_back(tsp.g.target(e));
    }
    list<DNode> path_list = path_s(source, adj);
    vector<DNode> path { make_move_iterator(begin(path_list)), make_move_iterator(end(path_list)) };
    return path;
}

// Visualize a solution graph from graph data and path
void ViewTspRCircuit(TSP_Data_R &tsp, const vector<DNode> &bestRoute, string graphname, double delta, const vector<DNode> &postos)
{
    ArcColorMap acolor(tsp.g);  // color of edges
    ArcStringMap aname(tsp.g);  // name of edges
    //for (auto v : postos) {
        //tsp.vcolor[v] = CYAN;
    //} // Se tentar colocar cor fica estranho
    int k = 0;
    for (auto it = bestRoute.begin(); it != bestRoute.end(); ++it) {
        DNode u = *it; ++it;
        if(it != bestRoute.end()){
            DNode v = *it;
            Arc a; bool found=false;
            for(OutArcIt e(tsp.g,u); e != INVALID; ++e){
                if( tsp.g.source(e) == u && tsp.g.target(e) == v ){
                    a = e;
                    aname[a] = to_string(k);
                    acolor[a] = BLUE; ++k;
                    //tsp.vcolor[u] = tsp.vcolor[v] = BLUE;
                    //cout << k << ",";
                    found=true; break;
                }
            }
            if(!found){
                cout << "ERRO: Aresta de vertice " << tsp.vname[u] << " para " << tsp.vname[v] <<" nao encontrada." << endl;
                exit(1);
            }
            --it;
        }
        else
            break;
    }
    ostringstream out;
    out << "TSP with Refueling on graph " << graphname << ", d=" << delta <<" cost =" << solutionCost(tsp, bestRoute);
    ViewListDigraph(tsp.g, tsp.vname, tsp.posx, tsp.posy, tsp.vcolor, acolor, out.str());
}

void writeOutputFile (string outputfile, string graphname, vector<DNode> path, int delta, TSP_Data_R &tsp, DNode source, double ub, double lb, double elapsed_time, int max_time){
    ofstream myfile;
    myfile.open (outputfile);
    myfile << graphname << " " << path.size() << " " << delta << " " << tsp.vname[source] << " " << elapsed_time << " " << max_time << " " << lb << " " << ub << endl;
    for( auto u : path){
        myfile << tsp.vname[u] << endl;
    }
    myfile.close();
}