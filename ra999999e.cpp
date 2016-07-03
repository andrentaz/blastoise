// MC658 - Projeto e Análise de Algoritmos III
// ---- LABORATÓRIO TSPR ----
// Prof: Flavio Keidi Miyazawa
// PED: Rafael Arakaki

#include <iostream>
#include <float.h>
#include <gurobi_c++.h>
#include <lemon/list_graph.h>
#include <lemon/unionfind.h>
#include <lemon/hao_orlin.h>
#include <lemon/adaptors.h>
#include <lemon/connectivity.h>
#include "mygraphlib.h"
#include "algoritmos.h"

using namespace lemon;
using namespace std;

// -----------------------------------------------------------------------------
// Classe auxiliar para usar o Gumory Hu
// is stored in the branch and cut tree. And when we define separation routines,
// we can recover the pointer and access the problem data again.
class TSP_Data {
public:
    TSP_Data(ListGraph &graph,
        NodeStringMap &nodename,
        NodePosMap &posicaox,
        NodePosMap &posy,
        EdgeValueMap &eweight);
    ListGraph &g;
    int NNodes,NEdges;
    int max_perturb2opt_it; // maximum number of iterations for heuristic Perturb2OPT
    NodeStringMap &vname;
    EdgeStringMap ename;
    NodeColorMap vcolor;
    EdgeColorMap ecolor;
    EdgeValueMap &weight;
    NodePosMap &posx;
    NodePosMap &posy;
    AdjacencyMatrix AdjMat; // adjacency matrix
    vector<Node> BestCircuit; // vector containing the best circuit found
    double BestCircuitValue;
};

TSP_Data::TSP_Data(ListGraph &graph,
		   NodeStringMap &nodename,
		   NodePosMap &posicaox,
		   NodePosMap &posicaoy,
		   EdgeValueMap &eweight):
    g(graph),
    vname(nodename),
    ename(graph),
    vcolor(graph),
    ecolor(graph),
    weight(eweight),
    posx(posicaox),
    posy(posicaoy),
    AdjMat(graph,eweight,MY_INF), //
    BestCircuit(countEdges(graph)) 
{
    NNodes=countNodes(this->g);
    NEdges=countEdges(this->g);
    BestCircuitValue = DBL_MAX;
    max_perturb2opt_it = 3000; // default value
}

class subtourelim: public GRBCallback
{ 
    TSP_Data &tsp;
    ListGraph::EdgeMap<GRBVar>& x;
    const Node source;
    ListGraph::NodeMap<bool> &postos;
    int delta;
    double (GRBCallback::*solution_value)(GRBVar);
public:
    subtourelim(TSP_Data &tsp, ListGraph::EdgeMap<GRBVar>& x, const Node source, ListGraph::NodeMap<bool> &postos, int delta)
        : tsp(tsp),x(x),source(source),postos(postos),delta(delta)  {    }
protected:
    void callback()
    { // --------------------------------------------------------------------------------
    // get the correct function to obtain the values of the lp variables
    if  (where==GRB_CB_MIPSOL) // if this condition is true, all variables are integer
        {solution_value = &subtourelim::getSolution;}
    else if ((where==GRB_CB_MIPNODE) &&  
        (getIntInfo(GRB_CB_MIPNODE_STATUS)==GRB_OPTIMAL))// node with optimal fractional solution
        {solution_value = &subtourelim::getNodeRel;}
    else return; // return, as this code do not take advantage of the other options
    // --------------------------------------------------------------------------------
    // Stores the edges with fractional values and integer values
    vector<Edge> FracEdges,OneEdges;
    // produces a subgraph h of g, with edges e with x[e]==1 
    // contracted, so we can apply Gomory-Hu tree in a small graph
    ListGraph::EdgeMap<bool> one_filter(tsp.g, false); // start without any edge
    ListGraph::EdgeMap<bool> non_zero_filter(tsp.g, false); // start without any edge
    for (EdgeIt e(tsp.g); e != INVALID; ++e) {
        if ((this->*solution_value)(x[e]) > 1-MY_EPS)
            OneEdges.push_back(e); // stores the edges with x[e]==1
        else if ((this->*solution_value)(x[e]) > MY_EPS)
            FracEdges.push_back(e); // includes edges with 0 < x[e] < 1
    }// define the subgraph with edges that have x[e]==1

    try {
        // --------------------------------------------------------------------------------
        // Use union-find to contract nodes (to obtain graph where each componente of g is contracted)
        //for (int i=0;i<tsp.NNodes;i++) UFIndexToNode[i]=INVALID;
        ListGraph::NodeMap<int> aux_map(tsp.g);
        UnionFind<ListGraph::NodeMap<int> > UFNodes(aux_map);
        for (NodeIt v(tsp.g); v!=INVALID; ++v) UFNodes.insert(v);
        for (vector<Edge>::iterator e_it=OneEdges.begin(); e_it != OneEdges.end(); ++e_it)
            UFNodes.join(tsp.g.u(*e_it),tsp.g.v(*e_it));// No problem if they are in a same component
        // --------------------------------------------------------------------------------
        // Put in a separate set all edges that are not inside a component 
        vector<Edge> CrossingEdges;
        for (EdgeIt e(tsp.g); e != INVALID; ++e) 
        if (UFNodes.find(tsp.g.u(e)) != UFNodes.find(tsp.g.v(e)))
            CrossingEdges.push_back(e);
        // --------------------------------------------------------------------------------
        // Generate an inverted list UFIndexToNode to find the node that represents a component
        vector<bool> ComponentIndex(tsp.NNodes);
        vector<Node> Index2h(tsp.NNodes);
        for(int i=0;i<tsp.NNodes;i++) ComponentIndex[i]=false;
        for (NodeIt v(tsp.g); v!=INVALID; ++v) ComponentIndex[UFNodes.find(v)]=true;
        // --------------------------------------------------------------------------------
        // Generate graph of components, add one node for each component and edges
        ListGraph h;
        EdgeValueMap h_capacity(h); 
        for(int i=0;i<tsp.NNodes;i++)  // add nodes to the graph h
            if (ComponentIndex[i]) Index2h[i]=h.addNode(); 
        for (vector<Edge>::iterator e_it=FracEdges.begin(); e_it != FracEdges.end(); ++e_it){
            Node  u = tsp.g.u(*e_it),              v = tsp.g.v(*e_it),
            hu = Index2h[UFNodes.find(u)],   hv = Index2h[UFNodes.find(v)];
            Edge a = h.addEdge(hu , hv );         // add edges to the graph h
            h_capacity[a] = (this->*solution_value)(x[*e_it]);
        }
        // create a map to terminals
        NodeBoolMap HasTerminal(h);
        for (NodeIt uit(tsp.g); uit!=INVALID; ++uit)
        {
            // if it is not a station
            if(!postos[uit]) 
            {
                HasTerminal[Index2h[UFNodes.find(uit)]] = true;
            }
        }
        // add the source to the 'terminal'
        HasTerminal[Index2h[UFNodes.find(source)]] = true;
        // --------------------------------------------------------------------------------
        GomoryHu<ListGraph, EdgeValueMap> ght(h, h_capacity);   ght.run();
        // The Gomory-Hu tree is given as a rooted directed tree. Each node has
        // an arc that points to its father. The root node has father -1.
        // Remember that each arc in this tree represents a cut and the value of
        // the arc is the weight of the corresponding cut. So, if an arc has weight
        // less than 2, then we found a violated cut and in this case, we insert the
        // corresponding constraint.

        NodeBoolMap cutmap(h);
        for (NodeIt u(h); u != INVALID; ++u) {
            GRBLinExpr expr = 0;
            if (ght.predNode(u)==INVALID) continue; // skip the root node
            if (ght.predValue(u) > 2.0 - MY_EPS) continue; // value of the cut is good
            
            // this cut violate this conditions
            // check if there is a terminal in both components
            if (HasTerminal[u] && HasTerminal[ght.predNode(u)])
            {
                ght.minCutMap(u, ght.predNode(u), cutmap);  // now, we have a violated cut
            
                // Percorre as arestas que cruzam alguma componente e insere as que pertencem ao corte
                for (vector<Edge>::iterator e_it=CrossingEdges.begin();e_it!=CrossingEdges.end();++e_it){
                    Node u=tsp.g.u(*e_it), v=tsp.g.v(*e_it),
                    hu = Index2h[UFNodes.find(u)], hv=Index2h[UFNodes.find(v)];
                    if (cutmap[hu] != cutmap[hv])
                    expr += x[*e_it];
                }
                addLazy( expr >= 2 );
            }
        }

        // Insere a restricao de sub caminho -------------------------------
        // -----------------------------------------------------------------
        // Mapeia as arestas da solucao
        int uncovered = 0;
        EdgeBoolMap covered(tsp.g, false);
        EdgeBoolMap active(tsp.g, false);
        for (EdgeIt a(tsp.g); a!=INVALID; ++a)
        {
            if (getSolution(x[a]) > 1 - MY_EPS)
            {
                active[a] = true;
                uncovered++;
            }
            else
            {
                active[a] = false;
            }
        }
        
        // -----------------------------------------------------------------
        // Percorre os caminhos e insere as restricoes de consumo
        GRBLinExpr expr = 0;
        Node u = source;
        bool changed = false;
        while (uncovered > 0)
        {
            // u := vertice atual
            // v := vertice destino
            // a := arco percorrido
            for (ListGraph::IncEdgeIt ait(tsp.g,u); ait!=INVALID; ++ait)
            {
                // verifica se o arco esta na solucao
                Edge a(ait);
                if (active[a] && !covered[a])
                {
                    // soma a distancia percorrida
                    changed = true;
                    expr += x[a]*tsp.weight[a];
                    
                    // verifica se o node destino e um posto
                    Node v = tsp.g.v(a);
                    if (u == v) 
                    {
                        v = tsp.g.u(a);
                    }
                    if (postos[v])
                    {
                        // se for um posto, termina aqui uma restricao de consumo
                        addLazy( expr <= delta);
                        expr = 0;
                    }
                    
                    // atualiza o valor do vertice atual e diminui as arestas
                    // que ainda nao foram cobertas
                    u = v;
                    covered[a] = true;
                    uncovered--;
                    break;                        
                }
            }

            // check the changed
            if (changed) 
            {
                changed = false;
            }
            else // it is a dead end 
            {
                break;
            }
        }
        // cout << "TERMINANDO CALLBACK" << endl;
    } catch (...) {
      cout << "Error during callback..." << endl;
    }
  }
};

// ********** ALTERE DAQUI PARA BAIXO (MÉTODO BRANCH & BOUND/CUT) *****************
// -----------------------------------------------------------------------------
// Converte a solucao do PLI para um vetor de arestas
bool convertSol(ListGraph::EdgeMap<GRBVar> &x, vector<DNode> &sol, TSP_Data_R &tsp, TSP_Data &utsp, Node u, EdgeBoolMap &cover, int uncovered)
{
    
    if (uncovered==0)
        return true;

    bool changed = false;
    for (IncEdgeIt e(utsp.g, u); e!=INVALID; ++e)
    {
        if (BinaryIsOne(x[e].get((GRB_DoubleAttr_X))))
        {
            if (cover[e]) 
            {
                // get the nodes id
                Node v = utsp.g.v(e);
                if (u == v)
                {
                    v = utsp.g.u(e);
                }

                // push the solution
                sol.push_back(tsp.g.nodeFromId(utsp.g.id(v)));
                cover[e] = false;
                uncovered--;
                changed = convertSol(x, sol, tsp, utsp, v, cover, uncovered);
                if (!changed)
                {
                    cover[e] = true;
                    sol.pop_back();
                    uncovered++;
                }
            }
        }
    }

    return changed;
}


// -----------------------------------------------------------------------------
// Obs: As variáveis e restrições abaixo são apenas exemplo e não necessariamente deverão ser usadas no laboratório.

// Otimiza o problema TSP-R através de branch & bound/cut.
// Retorna true se encontrar alguma solução viável, false do contrário.

// ATENÇÃO: Não modifique a assinatura deste método.
bool brach_and_bound999999(TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos,
                           const DNode source,
                           int delta, int maxTime, vector<DNode> &sol, double &lbound){
    // Converte o TSP direcionado para um nao direcionado com duas arestas
    ListGraph graph;
    EdgeValueMap weights(graph);

    // Adiciona os nos
    for (ListDigraph::NodeIt u(tsp.g); u!=INVALID; ++u)
    {
        Node v = graph.addNode();
    }

    // Adiciona as arestas
    for (ListDigraph::ArcIt ait(tsp.g); ait!=INVALID; ++ait)
    {
        // pega os dois nos incidentes
        Arc a(ait);
        DNode u = tsp.g.source(a);
        DNode v = tsp.g.target(a);

        // cria a mesma aresta no grafo não direcionado
        Node gu = graph.nodeFromId(tsp.g.id(u));
        Node gv = graph.nodeFromId(tsp.g.id(v));

        // insere a aresta no grafo nao direcionado
        Edge e = graph.addEdge(gu, gv);
        
        // Atribui pesos as arestas
        weights[e] = tsp.weight[a];
    }

    NodeStringMap nodename(graph);
    NodePosMap posicaox(graph);
    NodePosMap posicaoy(graph);

    TSP_Data utsp(graph, nodename, posicaox, posicaoy, weights);

    // utiliza o convertido
    ListGraph::EdgeMap<GRBVar> x(graph);
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    // TODO: [Opcional] Comente a linha abaixo caso não queira inserir cortes durante a execução do B&B
    model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

    model.getEnv().set(GRB_IntParam_Seed, 0);
    model.set(GRB_StringAttr_ModelName, "TSPR - TSP with Refueling"); // name to the problem
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a minimization problem

    // Add one binary variable for each arc and also sets its cost in the objective function
    for (EdgeIt e(utsp.g); e!=INVALID; ++e) {
        char name[100];
        
        Edge edge(e);
        unsigned uid = utsp.g.id(utsp.g.u(edge));
        unsigned vid = utsp.g.id(utsp.g.v(edge));

        sprintf(name,"x_%s_%s",tsp.vname[tsp.g.nodeFromId(uid)].c_str(),tsp.vname[tsp.g.nodeFromId(vid)].c_str());
        x[e] = model.addVar(0.0, 1.0, utsp.weight[e],GRB_BINARY,name);
    }
    model.update(); // run update to use model inserted variables

    // converte os terminais e os postos
    vector<Node> uterminais;
    for (auto t : terminais)
    {
        unsigned tid = tsp.g.id(t);
        uterminais.push_back(utsp.g.nodeFromId(tid));
    }

    NodeBoolMap upostos(utsp.g, false);
    for (auto p: postos)
    {
        unsigned pid = tsp.g.id(p);
        // upostos.push_back(utsp.g.nodeFromId(pid));
        upostos[utsp.g.nodeFromId(pid)] = true;
    }

    // Adicione restrições abaixo

    // (1) Nós terminais devem ser visitados exatamente uma vez
    for (auto v : uterminais) {
        GRBLinExpr expr = 0;
        for (IncEdgeIt e(utsp.g,v); e!=INVALID; ++e){
            expr += x[e];
        }
        model.addConstr(expr == 2 );
    }

    // (3) Nó source sempre presente no início do caminho
    Node usource = utsp.g.nodeFromId(tsp.g.id(source));
    GRBLinExpr expr = 0;
    for (IncEdgeIt e(utsp.g,usource); e!=INVALID; ++e){
        expr += x[e];
    }
    model.addConstr(expr >= 1 );

    try {
        model.update(); // Process any pending model modifications.
        //if (maxTime >= 0) model.getEnv().set(GRB_DoubleParam_TimeLimit,maxTime);

        subtourelim cb = subtourelim(utsp , x, usource, upostos, delta);
        model.setCallback(&cb);

        // TODO: [Opcional] Pode-se utilizar o valor de uma solução heurística p/ acelerar o algoritmo B&B (cutoff value).
        //cutoff = tsp.BestCircuitValue-MY_EPS;
        double cutoff = 0.0;
        if (cutoff > MY_EPS) model.getEnv().set(GRB_DoubleParam_Cutoff, cutoff );
        model.update(); // Process any pending model modifications.
        model.optimize();

        // Obtém o status da otimização
        int status = model.get(GRB_IntAttr_Status);
        if(status == GRB_INFEASIBLE || status == GRB_INF_OR_UNBD){
            cout << "Modelo inviavel ou unbounded." << endl;
            return false;
        }

        // Limitante inferior e superior do modelo
        //lbound = model.get(GRB_DoubleAttr_ObjBoundC);

        if( model.get(GRB_IntAttr_SolCount) <= 0 ){
            cout << "Modelo nao encontrou nenhuma solucao viavel no tempo. LowerBound = " << lbound << endl;
            return false;
        }
        else if (status == GRB_OPTIMAL){
            if(verbose) cout << "O modelo foi resolvido ate a otimalidade." << endl;
        }
        else {
            if(verbose) cout << "O modelo encontrou uma solucao sub-otima (i.e. nao ha garantia de otimalidade)." << endl;
        }

        double custo_solucao = model.get(GRB_DoubleAttr_ObjVal);
        
        int uncovered=0;
        EdgeBoolMap cover(utsp.g, false);
        for (EdgeIt e(utsp.g); e!=INVALID; ++e)
        {
            if (BinaryIsOne(x[e].get(GRB_DoubleAttr_X)))
            {
                cover[e] = true;
                uncovered++;
            }
        }
        sol.push_back(tsp.g.nodeFromId(utsp.g.id(usource)));        
        convertSol(x, sol, tsp, utsp, usource, cover, uncovered);

        // Calculo manual do custo da solução (deve ser igual ao ObjVal do Gurobi).
        double soma=0.0;
        ArcName aname(tsp.g);
        vector<Arc> edgesSol;
        ArcColorMap acolor(tsp.g);
        // if( verbose ) cout << "####### " << endl << "Edges of Solution (B&B):" << endl;
        // for (EdgeIt e(utsp.g); e!=INVALID; ++e){
        //     if (BinaryIsOne(x[e].get(GRB_DoubleAttr_X))){ // Note que se este método serve para variáveis binárias, p/ inteiras terá de usar outro método.
        //         soma += utsp.weight[e];
        //         edgesSol.push_back(tsp.g.arcFromId(utsp.g.id(e)));
        //         if( verbose) cout << "(" << tsp.vname[tsp.g.nodeFromId(utsp.g.id(utsp.g.u(e)))] << "," << tsp.vname[tsp.g.nodeFromId(utsp.g.id(utsp.g.v(e)))] << ")" << endl;
        //         acolor[tsp.g.arcFromId(utsp.g.id(e))] = BLUE;
        //     }
        // }
        // if( verbose ) cout << "####### " << endl;
        if( verbose ) cout << "####### " << endl << "Edges of Solution (B&B):" << endl;
        DNode u = sol[0];
        for (int i=1; i<sol.size(); i++) 
        {
            DNode v = sol[i];
            soma += tsp.AdjMatD.Cost(u,v);
            if ( verbose ) cout << "(" << tsp.vname[u] << "," << tsp.vname[v] << ")" << endl;
            u = v;
        }
        if( verbose ) cout << "####### " << endl;

        if( verbose ) cout << "Custo calculado pelo B&B = "<< soma << " / " << custo_solucao << endl;
        if( verbose ){
            cout << "Caminho encontrado a partir do vértice de origem (" << tsp.vname[source] << "): ";
            for(auto node : sol){
                cout << tsp.vname[node] << " ";
            } // Obs: O caminho é gerado a partir do nó source, se o conjunto de arestas retornado pelo B&B for desconexo, o caminho retornado por 'path_search' será incompleto.
            cout << endl << "Custo calculado da solucao (caminho a partir do no origem) = " << solutionCost(tsp, sol) << endl;
            ostringstream out;
            out << "TSP with Refueling B&B, cost= " << custo_solucao;
            ViewListDigraph(tsp.g, tsp.vname, tsp.posx, tsp.posy, tsp.vcolor, acolor, out.str());
        }
        return true;
    }
    catch(GRBException e) {
        cerr << "Gurobi exception has been thrown." << endl;
        cerr << "Error code = " << e.getErrorCode() << endl;
        cerr << e.getMessage();
    }
    catch (...) {
        cout << "Model is infeasible"  << endl;
        return false;
    }
    return false;
}
// ********** ALTERE DAQUI PARA CIMA (MÉTODO BRANCH & BOUND/CUT) *****************
