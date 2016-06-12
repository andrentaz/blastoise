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

bool validaRota(vector<DNode> rota, TSP_Data_R tsp, const vector<DNode> terminais, const vector<DNode> postos, int delta);
vector<DNode> heuristica_prob (TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, const double tempoRestante);
vector<DNode> init_solution (TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, const double tempoRestante);
vector<DNode> neighborTwoOpt (vector<DNode> rota, vectorTSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta);

// Otimiza o problema TSP-R através de uma heurística que utiliza grafo de vizinhanças (busca tabu, grasp, multi-start local search, etc)

// ATENÇÃO: Não modifique a assinatura deste método.
bool heuristica_hv999999(TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, int maxTime, vector<DNode> &sol, double &lbound){

    // Solucao inicial
    vector<DNode> S_acc = init_solution(tsp, terminais, postos, source, delta maxTime);

    // Melhor solução inicial: vazia e de custo infinito.
    vector<DNode> S_best = S_acc;
    double S_best_cost = DBL_MAX;

    // Essa heurística não tem limitante inferior
    lbound = 0.0;
    
    // Valores para funcionamento do Metropolis
    const double NEULER = exp(1.0);
    const double ALPHA = 0.95;
    const int K = 7;

    clock_t before = clock();
    double elapsed_time = (double) (clock()-before) / CLOCKS_PER_SEC;
    // Laço principal: executa até o limite de tempo.
    for ( double T=100.0; T > 0.001; T *= ALPHA ) {
        for (int i=0; i<100; ++i) {
            // Pega uma solucao da vizinhanca
            vector<DNode> S_new = neighborTwoOpt(S_acc, tsp, terminais, postos, delta);
            
            // Calcula a diferenca de custo
            double S_new_cost = solutionCost(tsp, S_new);
            double S_acc_cost = solutionCost(tsp, S_acc);
            double delta_S =  S_new_cost - S_acc_cost;
            
            // Atualiza a solucao se:
            // - a solucao nova tem um valor melhor que a solucao atual
            // - c.c. atualiza com uma probabilidade e^-((c(S_new)-c(S_acc))/kT)
            if( delta_S < 0 || (rand()/(double)RAND_MAX <= pow(NEULER, -delta_S/(K*T))){
                S_acc = std::move(S_new);
                if( S_new_cost  < S_best_cost ){
                    S_best = std::move(S_new); // C++11 (equivalente a S_best=S_new; S_new.clear();)
                    S_best_cost = S_new_cost;
                }
            }
        }
    }
    sol = S_best;
    if(sol.size())
        return true;
    return false;
}

vector<DNode> init_solution (TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, const DNode source, int delta, const double tempoRestante){
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
        vector<DNode> A; // lista de temminais alcancaveis do ultimo node do caminho que nao esgotam o combustivel
        for(auto t : terminaisNaoCobertos){
            if( usoAtualCombustivel + tsp.AdjMatD.Cost(q,t) + tsp.AdjMatD.Cost(t,pp[t]) <= delta -MY_EPS ){
                A.push_back(t);
            }
        }
        if(A.size() > 0){   // se eh possivel ir a algum terminal que tenha posto proximo
            //cout << "A.size=" << A.size() << "\t";
            unsigned long randomIndex = rand() % A.size();  // vai para algum dos terminais
            DNode v = A[randomIndex];
            rota.push_back(v);
            usoAtualCombustivel += tsp.AdjMatD.Cost(q,v);
            q = v;
            for(auto it = terminaisNaoCobertos.begin(); it != terminaisNaoCobertos.end(); it++){
                if( *it == v){ terminaisNaoCobertos.erase(it); break; }
            }
        }
        else{   // nao eh possivel ir a nenhum terminal
            vector<DNode> B; // lista de postos alcancaveis do ultimo vertice
            for(auto p : postos){
                if( usoAtualCombustivel + tsp.AdjMatD.Cost(q,p) <= delta -MY_EPS && p != q ){
                    B.push_back(p);
                }
            }

            //cout << "B.size=" << B.size() << " gasto atual: " << usoAtualCombustivel << " node corrente: " << tsp.g.id(q) << " nao cobertos=" << terminaisNaoCobertos.size() << endl;

            if(B.size() > 0){
                unsigned long randomIndex = rand() % B.size(); // vai para algum posto
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
 
/// <sumary>
/// Escolhe duas arestas partindo de um posto e tenta permutar 
/// </sumary>
vector<DNode> neighborTwoOpt (vector<DNode> rota, TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, int delta){
    DNode u, v; // nodes postos
    DNode x, y; // nodes terminais
    vector<DNode> postosNoCircuito;
    
    // Lista de postos no circuito encontrado
    for (auto p : postos) {
        // Procura na rota
        if (std::find(rota.begin(), rota.end(), p) != rota.end()) {
            postosNoCircuito.push_back(p);
        }
    }

    // Enquanto ainda tiver postos para escolher no circuito
    while( postosNoCircuito.size() != 0){
        unsigned long randomIndex = rand() % postosNoCircuito.size();
        u = postosNoCircuito[randomIndex];
        
        // Pega a posicao do posto na rota        
        unsigned long index_u = (std::find(rota.begin(), rota.end(), u) - rota.begin());
        unsigned long index_v = index_u+1;
        
        // Seleciona a aresta que saiu do posto
        double radius = tsp.AdjMatD.Cost(u,v);
        
        vector<DNode> A; // lista de terminais dentro do raio da aresta (u,v)
        
        // Escolhe os vertices terminais que estao no raio
        for (auto t : terminais) {
            // Se esta no caminho
            if (std::find(rota.begin(), rota.end(), t) != rota.end()) {
                // Se esta dentro do raio
                if (tsp.AdjMatD.Cost(u,t) <= radius) {
                    A.push_back(t);
                }
            }    
        }
            
        if(A.size() > 0){   // se existem terminais no raio
            //cout << "A.size=" << A.size() << "\t";
            randomIndex = rand() % A.size();  // vai para algum dos terminais
            x = A[randomIndex];
            unsigned long index_x = (std::find(rota.begin(), rota.end(), x) - rota.begin());
            unsigned long index_y = index_x+1;
            
            // pega nova rota
            vector<DNode> novaRota = twoOptSwap(rota, index_v, index_x);
            
            // verifica se nova rota eh valida
            if (validaRota(novaRota, tsp, terminais, postos, delta)) {
                return novaRota;
            } else {
                // pega nova rota
                novaRota = twoOptSwap(rota, index_v, index_y);
                
                // verifica se nova rota eh valida
                if (validaRota(novaRota, tsp, terminais, postos, delta)) {
                    return novaRota;
                } else {
                    // falhou!
                    // tenta mexer em outro posto
                    for (auto it=postosNoCircuito.begin(); it!= postosNoCircuito.end(); ++it) {
                        if (*it == u) {
                            postosNoCircuito.erase(it);
                            continue;
                        }
                    }
                }
            }
        }
    }
    
    // falhou a heuristica 2OPT
    return rota.clear();
} 

vector<DNode> twoOptSwap (vector<DNode> rota, unsigned long i, unsigned long k) {
    // rota usando arestas
    vector<DNode> novaRota;

    // 1. Copia rota[0] ate rota[i-1] na nova rota
    for ( int c = 0; c <= i - 1; ++c )
    {
        novaRota.push_back( rota[c] );
    }
     
    // 2. Copia rota[i] ate rota[k] ao contrario na nova rota
    int dec = 0;
    for ( int c = i; c <= k; ++c )
    {
        novaRota.push_back( rota[k - dec] );
        dec++;
    }
 
    // 3. Copia de rota[k+1] ate o fim da rota na nova rota 
    for ( int c = k + 1; c < size; ++c )
    {
        novaRota.push_back( rota[c] );
    }
    
    return novaRota;
    
}

bool validaRota(vector<DNode> rota, TSP_Data_R tsp, const vector<DNode> terminais, const vector<DNode> postos, int delta) {
    vector<DNode> terminaisNaoCobertos = terminais;
    double usoAtualCombustivel = 0.0;
    DNode q = rota[0];

    // Verifica node por node
    for (auto it=rota.begin(); it!=INVALID; ++it) {
        DNode u = *it;
        
        // checa distancia
        if ( usoAtualCombustivel + tsp.AdjMatD.Cost(q,u) >= delta ) {
            // Rota Invalida - estourou o limite de combustivel
            return false;
        } else if (u == q) {
            // Rota Invalida - chegou no source sem passar por todos os nodes
            return false;
        } else {
            // Atualiza combustivel gasto e verifica se eh terminal
            usoAtualCombustivel += tsp.AdjMatD.Cost(q,u);
            
            // se for terminal, cobre ele
            for (auto v=terminaisNaoCobertos.begin(); v!=INVALID; ++v) {
                if (u == *v) { terminaisNaoCobertos.erase(v); break;}
            }
            
            // se for um posto, abastece
            for (auto v=postos.begin(); v!=INVALID; ++v) {
                if (u == *v) {usoAtualCombustivel = 0.0; break;}
            }
        }
    }
    
    // terminou de analisar a rota, ve se cobriu todos os terminais
    if (terminaisNaoCobertos.size() > 0) {
        // Rota Invalida - nao cobriu todos os terminais
        return false;
    }
    
    return true;
}

// ********** ALTERE DAQUI PARA CIMA (HEURISTICA DE VIZINHANÇA) *****************
