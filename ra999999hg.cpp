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

// Rodando com multiplas threads
// #ifndef _OPENMP
// #define _OPENMP
// #include <omp>
// #endif

// Inclui a API do BRKGA
#include <algorithm>
#include "brkgaAPI/BRKGA.h"
#include "brkgaAPI/MTRand.h"

// Inclui o decoder
#include "TSPRDecoder.h"
#include "TSPRSolver.h"

bool check_feasible(vector<DNode> sol, TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, 
                    const DNode source, int delta); 

// ATENÇÃO: Não modifique a assinatura deste método.
bool heuristica_hg999999(TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos,
                         const DNode source,
                         int delta, int maxTime, vector<DNode> &sol, double &lbound) {

    // Inicializa as variaveis que serao utilizadas
    const clock_t begin = clock();
    
    const DNode s = source;

    // Decoder
    TSPRDecoder decoder(tsp, terminais, postos, source, delta);

    const long unsigned rngSeed = time(0);  // seed para gerar rand nums
    MTRand rng(rngSeed);    // gerador de numeros aleatoreos

    // caracteristicas do BRKGA
    const unsigned n = terminais.size();    // tamanho dos chromossomos
    const unsigned p = 256;     // tamanho da populacao
    const double pe = 0.10;     // fracao da populao que sera elite
    const double pm = 0.10;     // fracao da populacao que sera mutante
    const double rhoe = 0.70;   // probabilidade de herdar do pai elite
    const unsigned K = 5;       // numero de populacoes independentes
    const unsigned MAXT = 4;    // numero de threads para decode paralelo

    // inicializa a heuristica BRKGA
    BRKGA< TSPRDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);
    
    // BRKGA inner loop (evolution) configuration: Exchange top individuals
	const unsigned X_INTVL = 100;	    // troca os melhores individuos a cada 100 geracoes
	const unsigned X_NUMBER = 2;	    // troca os dois melhores
	const unsigned MAX_GENS = 5000;	    // roda para 1000 geracoes 

    // configuracao de evolucao do BRKGA: estrategia de restart
    unsigned relevantGeneration = 1;    // ultima geracao relevante
    const unsigned RESET_AFTER = 2000;
    vector< double > bestChromosome;
    double bestFitness = std::numeric_limits< double >::max();

    // Imprime informacoes sobre o multi-threading
    #ifdef  _OPENMP
        cout << "Running for " << MAX_GENS << " generation using " << MAXT
            << " out of " << omp_get_max_threads()
            << " available thread units..." << endl;
    #endif
    #ifndef _OPENMP
   	    cout << "Running for " << MAX_GENS
				<< " generations without multi-threading..." << endl;
    #endif

    // Roda o loop de evolucoes
    unsigned generation = 1;        // geracao corrente
    do {
        algorithm.evolve();         // evolui a populacao em uma geracao

        // Bookeeping: melhorou a melhor solucao
        if (algorithm.getBestFitness() < bestFitness) {
            // Salva a melhgor solucao para usar depois na cadeia evolutiva
            relevantGeneration = generation;
            bestFitness = algorithm.getBestFitness();
            bestChromosome = algorithm.getBestChromosome();

            cout << "\t" << generation
                << ") Improved best solution so far: "
                << bestFitness << endl;
        }

        // Estrategia de evolucao: restart
        if (generation - relevantGeneration > RESET_AFTER) {
            algorithm.reset();      // reinicia o algoritmo com chaves randomicas
            relevantGeneration = generation;

            cout << "\t" << generation << ") Reset at generation "
                << generation << endl;
        }

        // Estrategia de evolucao: troca dos melhores individuos entre populacoes
        if (generation % X_INTVL == 0 && relevantGeneration != generation) {
            algorithm.exchangeElite(X_NUMBER);

            cout << "\t" << generation 
                << ") Exchange top individuals." << endl;
        }

        // Proxima geracao?
        ++generation;
    } while (generation < MAX_GENS);

    // Imprime o fitness dos top 10 individuos de cada populacao:
	cout << "Fitness of the top 10 individuals of each population:" << endl;
	const unsigned bound = std::min(p, unsigned(10));	// makes sure we have 10 individuals
	for(unsigned i = 0; i < K; ++i) {
		cout << "Population #" << i << ":" << endl;
		for(unsigned j = 0; j < bound; ++j) {
			cout << "\t" << j << ") "
					<< algorithm.getPopulation(i).getFitness(j) << endl;
		}
	}

    // Reconstroi a melhor solucao
    TSPRSolver bestSolution (tsp, terminais, postos, s, delta, bestChromosome);
    bool retVal;

    // Se encontrou solucao reconstroi
    if (bestSolution.getTourDistance() != DBL_MAX) {
        // Imprime sua distancia
        cout << "Best solution found has objective value = "
                << bestSolution.getTourDistance() << endl;

        // Imprime o melhor tour
        cout << "Best tour:";
        for(auto pt : bestSolution.getTour()) {
            cout << " " << tsp.g.id(pt);
        }
        cout << endl;
        retVal = true;
    } else {
        // caso contrario, chora
        cout << "Unable to find a feasible solution!" << endl;
        retVal = false;
    }
    
    // Checa para ver se de fato é uma solucao viavel
    // Checando
    cout << "Double checking the solution" << endl;
    retVal = check_feasible(bestSolution.getTour(), tsp, terminais, postos, source, delta);
        
    if (retVal) {
        cout << "Feasible solution" << endl;
    }
    
    // Atribui para a solucao
    sol = bestSolution.getTour();

	const clock_t end = clock();
	cout << "BRKGA run finished in " << (end - begin) / double(CLOCKS_PER_SEC) << " s." << endl;

    return retVal;
}

bool check_feasible(vector<DNode> sol, TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos, 
                    const DNode source, int delta) {
    bool retVal = true;
    double fuelUsed = 0.0;
    
    vector<DNode> terminaisFaltantes = terminais;
    
    // Percorre o circuito
    DNode u = *(sol.begin());
    DNode v;
    DNode lastStation = source;
    vector<DNode>::iterator lts = sol.begin();
    for (vector<DNode>::iterator it=sol.begin()+1; it!=sol.end(); ++it) {
        bool isPosto = false;
        v = *it;
        
        // soma com o combustivel ja gasto
        fuelUsed += tsp.AdjMatD.Cost(u,v);
        
        // ve se nao ultrapassa o limite de combustivel
        if (fuelUsed < delta) {
            // ve se e um posto
            for (auto p : postos) {
                if (p==v) {
                    lastStation = v;
                    lts = it;
                    fuelUsed = 0.0;
                    u=v;
                    isPosto = true;
                    break;
                }
            }
            
            // se nao e posto e um terminal
            if (!isPosto) {
                for (auto t = terminaisFaltantes.begin(); t!=terminaisFaltantes.end(); ++t) {
                    if (*t==v) { terminaisFaltantes.erase(t); u=v; break; }
                }
            }
        } else {
            cout << "Infeasible solution by fuel limit!" << endl;
            
            // imprime o trecho que estourou o limite de combustivel
            for (;lts!=it; ++lts) {
                cout << " " << tsp.g.id(*lts);
            }
            cout << "\tFuel used: " << fuelUsed << "\t Fuel limit: "
                << delta << endl;
            
            retVal = false;
            break;
        }
    }
    
    if (retVal && terminaisFaltantes.size() != 0) {
        cout << "Infeasible solution for not covering all the terminals" << endl;
        // imprime os terminais faltantes
        for (auto t : terminaisFaltantes) {
            cout << " " << tsp.g.id(t);
        }
        cout << endl;
        retVal = false;
    }
    
    return retVal;
}

// ********** ALTERE DAQUI PARA CIMA (HEURISTICA GENETICA) *****************