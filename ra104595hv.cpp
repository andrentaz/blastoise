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

// ********** ALTERE DAQUI PARA BAIXO (HEURISTICA DE VIZINHANÇA) *****************

vector<DNode> heuristica_prob (TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, const double tempoRestante);

// Otimiza o problema TSP-R através de uma heurística que utiliza grafo de vizinhanças (busca tabu, grasp, multi-start local search, etc)

// ATENÇÃO: Não modifique a assinatura deste método.
bool heuristica_hv999999(TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, int maxTime, vector<DNode> &sol, double &lbound){

    // Melhor solução inicial: vazia e de custo infinito.
    vector<DNode> S_best;
    double S_best_cost = DBL_MAX;

    // Essa heurística não tem limitante inferior
    lbound = 0.0;

    clock_t before = clock();
    double elapsed_time = (double) (clock()-before) / CLOCKS_PER_SEC;
    // Laço principal: executa até o limite de tempo.
    while ( elapsed_time < maxTime ){
        vector<DNode> S_new = heuristica_prob (tsp, terminais, postos, source, delta, maxTime - elapsed_time);
        if( S_new.size() > 0 ){
            double S_new_cost = solutionCost(tsp, S_new);
            if( S_new_cost +MY_EPS < S_best_cost ){
                S_best = std::move(S_new); // C++11 (equivalente a S_best=S_new; S_new.clear();)
                S_best_cost = S_new_cost;
            }
        }
        elapsed_time = (double) (clock()-before) / CLOCKS_PER_SEC;
    }
    sol = S_best;
    if(sol.size())
        return true;
    return false;
}

vector<DNode> heuristica_prob (TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, const double tempoRestante){
    vector<DNode> terminaisNaoCobertos = terminais;
    vector<DNode> rota; // Inicia a rota com o nó origem (note que nó origem é um posto)
    rota.push_back(source);
    DNode q = source; // último inserido na rota
    double usoAtualCombustivel = 0.0;

    // Descobre pp[i] (posto mais próximo) para cada terminal t \in T
    ListDigraph::NodeMap <DNode> pp(tsp.g);
    for(auto t : terminais){
        double min = DBL_MAX;
        for(auto p : postos){
            if( tsp.AdjMatD.Cost(t,p) <= min){
                pp[t] = p;
                min = tsp.AdjMatD.Cost(t,p);
            }
        }
    }

    double elapsed_time = 0.0;
    clock_t before = clock();
    // Monta uma rota que atende a todos os terminais
    while(terminaisNaoCobertos.size() > 0 && elapsed_time < tempoRestante){
        vector<DNode> A;
        for(auto t : terminaisNaoCobertos){
            if( usoAtualCombustivel + tsp.AdjMatD.Cost(q,t) + tsp.AdjMatD.Cost(t,pp[t]) <= delta -MY_EPS ){
                A.push_back(t);
            }
        }
        if(A.size() > 0){
            //cout << "A.size=" << A.size() << "\t";
            unsigned long randomIndex = rand() % A.size();
            DNode v = A[randomIndex];
            rota.push_back(v);
            usoAtualCombustivel += tsp.AdjMatD.Cost(q,v);
            q = v;
            for(auto it = terminaisNaoCobertos.begin(); it != terminaisNaoCobertos.end(); it++){
                if( *it == v){ terminaisNaoCobertos.erase(it); break; }
            }
        }
        else{
            vector<DNode> B;
            for(auto p : postos){
                if( usoAtualCombustivel + tsp.AdjMatD.Cost(q,p) <= delta -MY_EPS && p != q ){
                    B.push_back(p);
                }
            }

            //cout << "B.size=" << B.size() << " gasto atual: " << usoAtualCombustivel << " node corrente: " << tsp.g.id(q) << " nao cobertos=" << terminaisNaoCobertos.size() << endl;

            if(B.size() > 0){
                unsigned long randomIndex = rand() % B.size();
                DNode v = B[randomIndex];
                rota.push_back(v);
                usoAtualCombustivel = 0.0;
                q = v;
            }
            else{
                // Heuristica falhou
                rota.clear();
                return rota;
            }
        }
        elapsed_time = (double) (clock() - before) / CLOCKS_PER_SEC;
    }
    // A rota deve retornar ao nó source para fechar o ciclo
    while(q != source && elapsed_time < tempoRestante){
        if( usoAtualCombustivel + tsp.AdjMatD.Cost(q,source) <= delta-MY_EPS){
            rota.push_back(source);
            usoAtualCombustivel += tsp.AdjMatD.Cost(q,source);
            q = source;
        }
        else{
            // Se não for possível, vá para um posto aleatório e tente novamente
            vector<DNode> B;
            for(auto p : postos){
                if( usoAtualCombustivel + tsp.AdjMatD.Cost(q,p) <= delta -MY_EPS && p != q ){
                    B.push_back(p);
                }
            }
            if(B.size() > 0){
                unsigned long randomIndex = rand() % B.size();
                DNode v = B[randomIndex];
                rota.push_back(v);
                usoAtualCombustivel = 0.0;
                q = v;
            }
            else{
                // Heuristica falhou
                rota.clear();
                return rota;
            }
        }
        elapsed_time = (double) (clock() - before) / CLOCKS_PER_SEC;
    }

    if( terminaisNaoCobertos.size() > 0 || q != source ){
        // Heuristica falhou por tempo
        rota.clear();
        return rota;
    }
    return rota;
}

// ********** ALTERE DAQUI PARA CIMA (HEURISTICA DE VIZINHANÇA) *****************
