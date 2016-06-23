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

// Inclui a API do BRKGA
#include <algorithm>
#include "brkgaAPI/BRKGA.h"
#include "brkgaAPI/MTRand.h"

// Inclui o decoder
#include "TSPRDecoder.h"
#include "TSPRSolver.h"

// ATENÇÃO: Não modifique a assinatura deste método.
bool heuristica_hg999999(TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos,
                         const DNode source,
                         int delta, int maxTime, vector<DNode> &sol, double &lbound) {

    // Inicializa as variaveis que serao utilizadas
    const clock_t begin = clock();

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
    const unsigned K = 3;       // numero de populacoes independentes
    const unsigned MAXT = 2;    // numero de threads para decode paralelo

    // inicializa a heuristica BRKGA
    BRKGA< TSPRDecoder, MTRand > algorithm(n, p, pe, pm, rhoe, decoder, rng, K, MAXT);

    // configuracao de evolucao do BRKGA: estrategia de restart
    unsigned relevantGeneration = 1;    // ultima geracao relevante
    const unsigned RESET_AFTER = 2000;
    vector< double > bestChromossome;
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
            bestChromossome = algorithm.getBestChromossome();

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
    TSPRSolver bestSolution (tsp, bestChromossome);
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
    

	const clock_t end = clock();
	cout << "BRKGA run finished in " << (end - begin) / double(CLOCKS_PER_SEC) << " s." << endl;

    return retVal;
}
// ********** ALTERE DAQUI PARA CIMA (HEURISTICA GENETICA) *****************