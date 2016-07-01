#include "group_tsp.h"
#include <chrono>
#include <limits>
#include <cmath>
#include <cstdio>
#include <cfloat>
#include <numeric>
#include <gurobi_c++.h>
#include <lemon/list_graph.h>
#include <lemon/gomory_hu.h>
#include <lemon/adaptors.h>
#include <lemon/dfs.h>
#include "mygraphlib.h"

#define BC_EPS 0.001
#define BC_INF 1000000000000.0

//=============================================================================
/**
 * This is the type used to obtain the pointer to the problem data.
 * This pointer is stored in the branch and cut tree.
 * And when we define separation routines,
 * we can recover the pointer and access the problem data again.
 */
class GTSP_Data {
public:
	const ListGraph &g;
	GTSP_Data(const ListGraph &graph) :
			g(graph) {
	}
};

//=============================================================================

class subtourelim: public GRBCallback {
	GTSP_Data &gtsp;
	ListGraph::EdgeMap<GRBVar>& x;
	ListGraph::NodeMap<GRBVar>& y;

	double (GRBCallback::*solution_value)(GRBVar);
public:
	subtourelim(GTSP_Data &tsp, ListGraph::EdgeMap<GRBVar>& x,
			ListGraph::NodeMap<GRBVar>& y) :
			gtsp(tsp), x(x), y(y) {
	}
protected:
	void callback() {
		if (where == GRB_CB_MIPSOL) // if integer solution
		{
			solution_value = &subtourelim::getSolution;
		} else
		// try to comment the next two lines and test the execution time
		if ((where == GRB_CB_MIPNODE)
				&& (getIntInfo(GRB_CB_MIPNODE_STATUS) == GRB_OPTIMAL)) // if node with optimal fractional solution
				{
			solution_value = &subtourelim::getNodeRel;
		} else
			return;
		try {

			typedef ListGraph::EdgeMap<double> capacityType;
			capacityType capacity(gtsp.g);

			for (ListGraph::EdgeIt e(gtsp.g); e != INVALID; ++e) {
				capacity[e] = getSolution(x[e]);
			}

			typedef ListGraph::NodeMap<double> SolutionY;
			SolutionY sol(gtsp.g);
			for (NodeIt u(gtsp.g); u != INVALID; ++u) {
				sol[u] = getSolution(y[u]);
			}

			GomoryHu<ListGraph, capacityType> ght(gtsp.g, capacity);
			ght.run();

			// The Gomory-Hu tree is given as a rooted directed tree. Each node has
			// an arc that points to its father. The root node has father -1.
			// Remember that each arc in this tree represents a cut and the value of
			// the arc is the weight of the corresponding cut. So, if an arc has weight
			// less than 2 than we found a violated cut and in this case, we insert the
			// corresponding constraint.
			double vcut;
			for (NodeIt u(gtsp.g); u != INVALID; ++u) {

				for (NodeIt v(gtsp.g, u); v != INVALID; ++v) {

					if (u != v) {
						// Check if nodes are in solution
						if (sol[u] >= (1 - BC_EPS) && sol[v] >= (1 - BC_EPS)) {
							vcut = ght.minCutValue(u, v);
							// If there is a edge between nodes, go to next.
							if (vcut > 2.0 * (sol[u] + sol[v] - 1) - BC_EPS)
								continue;
							// Otherwise add a restriction, stating that node
							// should at least degree equals 2.
							GRBLinExpr expr = 0;

							for (GomoryHu<ListGraph, EdgeWeight>::MinCutEdgeIt a(
									ght, u, v); a != INVALID; ++a)
								expr += x[a];
							addLazy(expr >= 2 * (y[u] + y[v] - 1));
						}

					}

				}
			}

		} catch (...) {
			cout << "Error during callback..." << endl;
		}
	}
};

//=============================================================================

int ExactGroupTSP(const ListGraph &g, const ListGraph::EdgeMap<double>& weights,
		const vector<set<ListGraph::Node> > &S, vector<ListGraph::Node> &sol,
		long max_time, double &best_time, double &LB, string &alg_info) {
	/**
	 * Realiza branch-and-cut para o Group TSP.
	 *
	 * Entrada:
	 * @param   g           grafo simples utilizado
	 * @param   weights     pesos das arestas
	 * @param   S           vetor de grupos de vertices (ver def. do problema)
	 * @param   max_time    tempo maximo (em seg) que o procedimento deve ocorrer
	 *
	 * Saida:
	 * @param   sol         sequencia de vertices que representa ciclo
	 * @param   best_time   momento em que solucao atual foi encontrada
	 * @param   LB          limite inferior encontrado para custo otimo
	 * @param   alg_info    informacoes de execucao do algoritmo, ex: cadeia de
	 *                      heuristicas utilizadas
	 *
	 * @return          0 = nao foi possivel encontrar solucao
	 *                  1 = solucao encontrada, mas nao necessariamente otima
	 *                  2 = solucao otima encontrada
	 */

	int seed = 1;
	srand48(seed);

	ListGraph::EdgeMap<GRBVar> x(g); // Variables that will represent if edge is IN (x=1) of OUT (x=0) the solution.
	ListGraph::NodeMap<GRBVar> y(g); // Variables that will represent if node is IN (y=1) of OUT (y=0) the solution.
	GRBEnv env = GRBEnv(); // Create Gurobi environment.
	GRBModel model = GRBModel(env); // Creater Gurobi Model.

	model.getEnv().set(GRB_IntParam_LazyConstraints, 1); //using lazy constraints
	model.getEnv().set(GRB_IntParam_Seed, seed);

	model.set(GRB_StringAttr_ModelName, "Group TSP"); // set problem's name
	model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // set a minimization problem

	// Add one binary variable for each edge and also sets its cost in the objective function
	char name[1000];
	NodeName vname(g);
	for (ListGraph::EdgeIt e(g); e != INVALID; ++e) {
		sprintf(name, "x_%s_%s", vname[g.u(e)].c_str(), vname[g.v(e)].c_str());
		x[e] = model.addVar(0.0, 1.0, weights[e], GRB_BINARY, name);
	}

	// Add one binary variable for each node and also sets its cost in the objective function
	for (ListGraph::NodeIt v(g); v != INVALID; ++v) {
		y[v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY);
	}

	model.update(); // run update to use model inserted variables

	// Add degree constraint for each node in the solution (sum of solution edges incident to a node is 2)
	for (ListGraph::NodeIt v(g); v != INVALID; ++v) {
		GRBLinExpr expr;
		for (ListGraph::IncEdgeIt e(g, v); e != INVALID; ++e)
			expr += x[e];
		model.addConstr(expr == 2 * y[v]);
	}

	// Add constraint for each Set provided (at least one node from each set should be in solution)
	set<ListGraph::Node>::iterator v;
	for (size_t i = 0; i < S.size(); i++) {
		GRBLinExpr expr;
		for (v = S[i].begin(); v != S[i].end(); v++)
			expr += y[*v];
		model.addConstr(expr >= 1);

	}

	GTSP_Data gtsp(g);

	try {
		// Process any pending model modifications.
		model.update();
		if (max_time >= 0)
			model.getEnv().set(GRB_DoubleParam_TimeLimit, max_time);

		// Set callback to be used.
		subtourelim cb = subtourelim(gtsp, x, y);
		model.setCallback(&cb);

		// Process any pending model modifications.
		model.update();
		model.optimize();

		LB = model.get(GRB_DoubleAttr_ObjBound);
		best_time = model.get(GRB_DoubleAttr_Runtime);

		int origin = 0;
		Node left = g.nodeFromId(origin);
		Node center = g.nodeFromId(origin);
		Node right;
		double solE;
		double sol_center;
		double sol_right;
		sol.clear();
		// Calculates Path.
		do {
			sol_center = y[center].get(GRB_DoubleAttr_X);

			if (sol_center > 0.99) {

				for (ListGraph::IncEdgeIt e(g, center); e != INVALID; ++e) {
					right = g.oppositeNode(center, e);

					sol_right = y[right].get(GRB_DoubleAttr_X);

					if (sol_right > 0.99) {

						solE = x[e].get(GRB_DoubleAttr_X);

						if (g.id(right) != g.id(left) && (solE > 0.99)) {
							sol.push_back(right);
							left = center;
							center = right;
							break;
						}

					}

				}

			} else {
				origin++;
				left = g.nodeFromId(origin);
				center = g.nodeFromId(origin);
			}

		} while (g.id(right) != origin);

	} catch (GRBException e) {
		cout << "Error code: " << e.getErrorCode() << "\n";
		cout << "Error code: " << e.getMessage() << "\n";
		return 0;
	}

	if (sol.size() > 0) {
		return 2;
	}

	return 2;

}

/**
 * Calculates weights for random process to each node in graph.
 * @param g[in]				Instance graph.
 * @param S[in]	Times 		Set of nodes from instance problem.
 * @param nodeWeight[out]	Weights for random choosing process.
 * @param weightsSum[in]	Sum of weights for each set.
 */
void calcNodeWeights(const ListGraph &g, const vector<set<ListGraph::Node> > &S,
		ListGraph::NodeMap<int> &nodeWeight, vector<int> &weightsSum) {

	set<ListGraph::Node>::iterator v;
	int cnt;
	double setWeightSum = 0;

	// It will create a weightSum for each set.
	if (weightsSum.size() != S.size())
		weightsSum.resize(S.size());

	for (size_t i = 0; i < S.size(); i++) {
		setWeightSum = 0;
		for (v = S[i].begin(); v != S[i].end(); v++) {

			cnt = 0;
			// Calculates weight of each node.
			for (ListGraph::IncEdgeIt e(g, *v); e != INVALID; ++e) {
				cnt++;
			}
			nodeWeight[*v] = cnt;
			// Calculates sum of weights for nodes in each set.
			setWeightSum += cnt;
		}
		weightsSum[i] = setWeightSum;
	}

}

/**
 * Choose a random node based on number of incident edges of each node.
 *
 * @param S				Set that represents a node group.
 * @param nodeWeight	Map containing weight of each node based on edges.
 * @param weightsSum	Sum of weight for the current group.
 */
ListGraph::Node chooseRandomNode(const set<ListGraph::Node> &S,
		ListGraph::NodeMap<int> &nodeWeight, int &weightsSum) {

	int rnd = rand() % weightsSum + 1;
	set<ListGraph::Node>::iterator v;

	for (v = S.begin(); v != S.end(); ++v) {

		if (rnd <= nodeWeight[*v]) {
			return *v;
		}
		rnd -= nodeWeight[*v];

	}

	cerr << "Error choosing random node.\n";
	exit(1);

}

/**
 * Change current node with a node from the same set choose randomly according
 * to its number of edges.
 *
 * @param[out] current		Current node reference. It will contains the new node choose randomly.
 * @param[in]  S			Vector of Sets representinf groups of nodes in problem.
 * @param[in]  nodeWeight	Map containing weight of each node based on edges.
 * @param[in]  nodeWeight	Vector containing sum of weights of each node group.
 */
void changeNodeRondom(ListGraph::Node &current,
		const vector<set<ListGraph::Node> > &S,
		ListGraph::NodeMap<int> &nodeWeight, vector<int> &weightsSum) {

	// Loop to the node group.
	for (size_t i = 0; i < S.size(); ++i) {

		if (S[i].find(current) != S[i].end()) {
			current = chooseRandomNode(S[i], nodeWeight, weightsSum[i]);
			return;
		}
	}

}

/**
 * Perform a local search in solution, trying change nodes.
 *
 * @param sol[out]			Current solution that could be altered in local search.
 * @param cur_min_cost[out]	Current minimun cost value. Can be altered.
 * @param S[in]	Times 		Set of nodes from instance problem.
 * @param nodeWeight[in]	Weights for random choosing process.
 * @param weightsSum[in]	Sum of weights for each set.
 * @param g[in]				Instance graph.
 * @param weights[in]		Edges weight.
 * @param maxTimes[in]		Maximum times to perform local search.
 * @param end[out]			Time when possible solution is found.
 */
void localSearch(vector<ListGraph::Node> &sol, double &cur_min_cost,
		const vector<set<ListGraph::Node> > &S,
		ListGraph::NodeMap<int> &nodeWeight, vector<int> &weightsSum,
		const ListGraph &g, const ListGraph::EdgeMap<double>& weights,
		const int maxTimes, chrono::time_point<chrono::system_clock> &end) {

	int searchs = 0;
	int i;
	int sz = sol.size();
	double cost = 0.0;
	vector<ListGraph::Node> cur(sol); // candidate to solution.

	while (searchs < maxTimes) {
		i = rand() % sz;
		changeNodeRondom(cur[i], S, nodeWeight, weightsSum);

		/* compute cost */
		cost = 0.0;
		for (int i = 0; i < sz; i++) {
			/* find the cheapest edge available */
			double mn = numeric_limits<double>::infinity();
			for (ListGraph::Edge e = findEdge(g, cur[i], cur[(i + 1) % sz]);
					e != INVALID;
					e = findEdge(g, cur[i], cur[(i + 1) % sz], e)) {
				mn = min(mn, weights[e]);
			}
			cost += mn;
		}
		// Change solution if current cost is better
		if (cost < cur_min_cost) {
			cur_min_cost = cost;
			sol = cur;
			end = chrono::system_clock::now();
		}
		++searchs;
	}

}

/**
 *
 * Add a random node to solution according to its degree. Nodes with
 * bigger degrees have bigger chances to be chosen.
 *
 * @param g[in]				Instance graph.
 * @param nodeWeight[in]	Weights for random choosing process.
 * @param seen[out]			Set that contains nodes in already in solution.
 * @param curr[out]			Current solution that could be altered in local search.
 *
 */
void addRandomNodes(const ListGraph &g, ListGraph::NodeMap<int> &nodeWeight,
		set<ListGraph::Node> &seen, vector<ListGraph::Node> &cur,
		const int &overallSum) {
	int rnd = rand() % overallSum + 1;

	for (NodeIt v(g); v != INVALID; ++v) {

		if (rnd <= nodeWeight[v]) {
			if (seen.find(v) == seen.end()) {
				seen.insert(v);
				cur.push_back(v);
				return;
			} else {
				addRandomNodes(g, nodeWeight, seen, cur, overallSum);
				return;
			}
		}
		rnd -= nodeWeight[v];

	}

	cerr << "Error adding random node.\n";
	exit(1);

}

int HeuristicGroupTSP(const ListGraph &g,
		const ListGraph::EdgeMap<double>& weights,
		const vector<set<ListGraph::Node> > &S, vector<ListGraph::Node> &sol,
		long max_time, double &best_time, double &LB, string &alg_info) {
	/**
	 * Computa solucao heuristica para o Group TSP.
	 *
	 * Entrada:
	 * @param   g           grafo simples utilizado
	 * @param   weights     pesos das arestas
	 * @param   S           vetor de grupos de vertices (ver def. do problema)
	 * @param   max_time    tempo maximo (em seg) que o procedimento deve ocorrer
	 *
	 * Saida:
	 * @param   sol         sequencia de vertices que representa ciclo
	 * @param   best_time   momento em que solucao atual foi encontrada
	 * @param   LB          limite inferior encontrado para custo otimo
	 * @param   alg_info    informacoes de execucao do algoritmo, ex: cadeia de
	 *                      heuristicas utilizadas
	 *
	 * @return          0 = nao foi possivel encontrar solucao
	 *                  1 = solucao encontrada, mas nao necessariamente otima
	 *                  2 = solucao otima encontrada
	 */

	/*****************************************************************************
	 ALGORITMO DE EXEMPLO
	 Enquanto ha tempo:
	 1) Seleciona aleatoriamente um vertice (unico) por grupo.
	 2) Permuta aleatoriamente os vertices selecionados.
	 3) Compara a solucao encontrada com a otima, substituindo-a se for melhor.
	 *****************************************************************************/
	const double threshold = 0.01;
	chrono::time_point<chrono::system_clock> start, end;
	start = chrono::system_clock::now();
	srand((unsigned int) time(NULL));

	best_time = 0;
	LB = 0;

	double cur_min_cost = numeric_limits<double>::infinity();

	ListGraph::NodeMap<int> nodesWeight(g);
	vector<int> weightsSum(S.size());
	calcNodeWeights(g, S, nodesWeight, weightsSum);
	int overAllSum = accumulate(weightsSum.begin(), weightsSum.end(), 0);

	int maxTimes = 30;

	while (end = chrono::system_clock::now(), (double) (chrono::nanoseconds(
			end - start)).count() / 1.0e9 + threshold < (double) max_time) {
		set<ListGraph::Node> seen;
		vector<ListGraph::Node> cur;
		for (const auto s : S) {
			if (s.size() == 0)
				continue;
			auto it = s.begin();
			advance(it, rand() % s.size());
			ListGraph::Node choice = *it;
			if (seen.find(choice) == seen.end()) {
				seen.insert(choice);
				cur.push_back(choice);
			}
		}

		// guarantee solution size at least 3.
		while (cur.size() < 3) {
			addRandomNodes(g, nodesWeight, seen, cur, overAllSum);
		}

		random_shuffle(cur.begin(), cur.end());
		/* compute cost */
		int sz = (int) cur.size();
		double cost = 0.0;
		for (int i = 0; i < sz; i++) {
			/* find the cheapest edge available */
			double mn = numeric_limits<double>::infinity();
			for (ListGraph::Edge e = findEdge(g, cur[i], cur[(i + 1) % sz]);
					e != INVALID;
					e = findEdge(g, cur[i], cur[(i + 1) % sz], e)) {
				mn = min(mn, weights[e]);
			}
			cost += mn;
		}
		if (cost < cur_min_cost) {
			cur_min_cost = cost;
			sol = cur;


			end = chrono::system_clock::now();

			// perform a local search changing nodes randomly.
			localSearch(sol, cur_min_cost, S, nodesWeight, weightsSum, g,
					weights, maxTimes, end);

			best_time = (double) (chrono::nanoseconds(end - start)).count()
					/ 1.0e9;
			cout << cur_min_cost << ";" << best_time << "\n";
		}
	}

	if (isfinite(cur_min_cost)) {
		return 1;
	} else {
		return 0;
	}
}
