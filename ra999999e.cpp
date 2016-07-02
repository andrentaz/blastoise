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

// Rotina de callback para inserir cortes durante a execução do B&B
class subtourelim: public GRBCallback
{
    TSP_Data_R &tsp;
    ListDigraph::ArcMap<GRBVar>& x;
    
    // Atributos para inserir restricoes de limite de caminho
    const DNode source;
    const vector<DNode> postos;
    int delta;    

    double (GRBCallback::*solution_value)(GRBVar);
public:
    subtourelim(TSP_Data_R &tsp, ListDigraph::ArcMap<GRBVar>& x, 
        const DNode source, const vector<DNode> postos, int delta) 
            : tsp(tsp),x(x),source(source),postos(postos),delta(delta)
    {    }
protected:
    void callback()
    {
        // -----------------------------------------------------------------
        // O metodo de callback tem duas funcoes:
        // 1. ele tem que garantir que o circuito e conexo
        // 2. ele tem que garantir que o caminho respeita a restricao de 
        // distancia entre os postos
        // -----------------------------------------------------------------
        if  (where==GRB_CB_MIPSOL) {
            // if mip solution (in this case, x is integer, ***but can be unconnected***)
            // getSolution is the function that obtain the values of x in MIPSOL
            solution_value = &subtourelim::getSolution;

        } else // it is not a candidate solution
            // if node with optimal fractional solution
            // getNodeRel obtain the optimum relaxation for that node in MIPNODE when OPTIMAL
            // Exercise: try to comment the next two lines and test the execution time
        if ((where==GRB_CB_MIPNODE) && (getIntInfo(GRB_CB_MIPNODE_STATUS)==GRB_OPTIMAL)){
            solution_value = &subtourelim::getNodeRel; }
        else return;
        // ******* ALTERE DAQUI PARA BAIXO (BRANCH AND CUT) ********** //
        try {
            cout << "INICIANDO CALLBACK" << endl;            
            // -----------------------------------------------------------------
            // Cria um grafo representando cada um dos circuitos fechados
            vector<Arc> FracArcs, OneArcs;
            for (ArcIt a(tsp.g); a!=INVALID; ++a)
            {
                if (getSolution(x[a]) > 1-MY_EPS) 
                {
                    OneArcs.push_back(a);   // Arestas com valor x[a]==1
                }
                else if (getSolution(x[a]) > MY_EPS)
                {
                    FracArcs.push_back(a);  // Arestas com valores fracionarios
                }
            }   // definindo os subgrafos

            // -----------------------------------------------------------------
            // Utiliza o algoritmo UnionFind para contrair as componentes conexas
            Digraph::NodeMap<int> aux_map(tsp.g);
            UnionFind< Digraph::NodeMap<int> > UFNodes(aux_map);
            for (Digraph::NodeIt v(tsp.g); v!=INVALID; ++v)
            {
                UFNodes.insert(v);
            }

            // Cria aqui as componentes contraidas
            for (vector<Arc>::iterator ait=OneArcs.begin(); ait!=OneArcs.end(); ++ait) 
            {
                // Mesmo se estiverem na mesma componente contraida
                UFNodes.join(tsp.g.source(*ait),tsp.g.target(*ait));
            }
            
            // Coloca num outro conjunto de arcos que cruzam as componentes
            vector<Arc> CrossingArcs;
            for (ArcIt a(tsp.g); a!=INVALID; ++a)
            {
                if (UFNodes.find(tsp.g.source(a)) != UFNodes.find(tsp.g.target(a)))
                {
                    CrossingArcs.push_back(a);
                }
            }
            // -----------------------------------------------------------------
            // Gera uma lista para obter de maneira rapida o no que representa
            // a componente
            vector<bool> ComponentIndex(tsp.NNodes, false);
            for (Digraph::NodeIt v(tsp.g); v!=INVALID; ++v)
            {
                ComponentIndex[UFNodes.find(v)] = true;
            }
            // -----------------------------------------------------------------
            // Gera o grafo de componentes conexas, adiciona um super no e uma aresta
            vector<DNode> Index2h(tsp.NNodes);
            Digraph h;
            ArcValueMap h_capacity(h);
            for (int i=0; i<tsp.NNodes; i++)
            {
                if (ComponentIndex[i])
                {
                    Index2h[i] = h.addNode();
                }
            }
            for (vector<Arc>::iterator ait=FracArcs.begin(); ait!=FracArcs.end(); ++ait)
            {
                DNode u = tsp.g.source(*ait);
                DNode v = tsp.g.target(*ait);
                DNode hu = Index2h[UFNodes.find(u)];
                DNode hv = Index2h[UFNodes.find(v)];

                Arc a = h.addArc(hu, hv);   // adiciona o arco ao grafo h
                h_capacity[a] = getSolution(x[*ait]);
            }
            // -----------------------------------------------------------------
            // Insere as restricoes de conexidade
            Tolerance<double> tolerance;
            HaoOrlin<Digraph, ArcValueMap, Tolerance<double> > hocut(h, h_capacity, tolerance);

            // Inicia o corte e coloca 

            // calcula o corte minimo de cada vertice de h (componente de g)
            for (Digraph::NodeIt uit(h); uit!=INVALID; ++uit)
            {
                DNode u(uit);
                hocut.init(u);
                hocut.calculateOut();
                Digraph::NodeMap<bool> outCutMap(h);

                // pega o corte minimo de saida e o seu valor
                double outcut = hocut.minCutMap(outCutMap);

                hocut.init(u);
                hocut.calculateIn();
                Digraph::NodeMap<bool> inCutMap(h);

                // pega o corte minimo de entrada e o seu valor
                double incut = hocut.minCutMap(inCutMap);

                // Pula 
                if (outcut > 1 - MY_EPS) 
                    if (incut > 1 - MY_EPS)
                        if (outcut - incut < MY_EPS)
                            continue;   // o corte e valido
                
                // Trata o caso de um corte que viola a condicao
                // Percorre os arcos que cruzam componentes e insere as que pertence ao corte de saida
                GRBLinExpr exprout = 0;
                for (vector<Arc>::iterator ait=CrossingArcs.begin(); ait!=CrossingArcs.end(); ++ait)
                {
                    DNode u = tsp.g.source(*ait);
                    DNode v = tsp.g.target(*ait);
                    DNode hu = Index2h[UFNodes.find(u)];
                    DNode hv = Index2h[UFNodes.find(v)];

                    // verifica quais arcos estao no corte de saida
                    if (outCutMap[hu]!=outCutMap[hv])
                    {
                        exprout += x[*ait];
                    }
                }
                addLazy( exprout >= 1);
                
                // Percorre os arcos que cruzam componentes e insere as que pertence ao corte de entrada
                GRBLinExpr exprin = 0;
                for (vector<Arc>::iterator ait=CrossingArcs.begin(); ait!=CrossingArcs.end(); ++ait)
                {
                    DNode u = tsp.g.source(*ait);
                    DNode v = tsp.g.target(*ait);
                    DNode hu = Index2h[UFNodes.find(u)];
                    DNode hv = Index2h[UFNodes.find(v)];

                    // verifica quais arcos estao no corte de saida
                    if (inCutMap[hu]!=inCutMap[hv])
                    {
                        exprin += x[*ait];
                    }
                }
                addLazy( exprin >= 1);

                // Por fim, faz com que o numero dos arcos que saem e os que entram sejam iguais
                addLazy( exprout - exprin == 0);

            }
            // -----------------------------------------------------------------
            // Insere a restricao de sub caminho -------------------------------
            // -----------------------------------------------------------------
            // Mapeia as arestas da solucao
            int uncovered = 0;
            Digraph::ArcMap<bool> active(tsp.g);
            for (ArcIt a(tsp.g); a!=INVALID; ++a)
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
            
            // Mapeia os postos
            Digraph::NodeMap<bool> stations(tsp.g, false);
            for (auto p : postos)
            {
                stations[p] = true;
            }
            // -----------------------------------------------------------------
            // Percorre os caminhos e insere as restricoes de consumo
            GRBLinExpr expr = 0;
            DNode u = source;
            while (uncovered > 0)
            {
                // u := vertice atual
                // v := vertice destino
                // a := arco percorrido
                for (Digraph::OutArcIt ait(tsp.g,u); ait!=INVALID; ++ait)
                {
                    // verifica se o arco esta na solucao
                    Arc a(ait);
                    if (active[a])
                    {
                        // soma a distancia percorrida
                        expr += x[a]*tsp.weight[a];
                        
                        // verifica se o node destino e um posto
                        DNode v = tsp.g.target(a);
                        if (stations[v])
                        {
                            // se for um posto, termina aqui uma restricao de consumo
                            addLazy( expr <= delta);
                            expr = 0;
                        }
                        
                        // atualiza o valor do vertice atual e diminui as arestas
                        // que ainda nao foram cobertas
                        u = v;
                        uncovered--;
                        break;                        
                    }
                }
            }
            cout << "TERMINANDO CALLBACK" << endl;
        }
        catch (...) {
            cout << "Error during callback..." << endl;
        }
        // ******* ALTERE DAQUI PARA CIMA (BRANCH AND CUT) ********** //
    }
};

// ********** ALTERE DAQUI PARA BAIXO (MÉTODO BRANCH & BOUND/CUT) *****************
// Obs: As variáveis e restrições abaixo são apenas exemplo e não necessariamente deverão ser usadas no laboratório.

// Otimiza o problema TSP-R através de branch & bound/cut.
// Retorna true se encontrar alguma solução viável, false do contrário.

// ATENÇÃO: Não modifique a assinatura deste método.
bool brach_and_bound999999(TSP_Data_R &tsp, const vector<DNode> &terminais, const vector<DNode> &postos,
                           const DNode source,
                           int delta, int maxTime, vector<DNode> &sol, double &lbound){
    ListDigraph::ArcMap<GRBVar> x(tsp.g);
    GRBEnv env = GRBEnv();
    GRBModel model = GRBModel(env);

    // TODO: [Opcional] Comente a linha abaixo caso não queira inserir cortes durante a execução do B&B
    model.getEnv().set(GRB_IntParam_LazyConstraints, 1);

    model.getEnv().set(GRB_IntParam_Seed, 0);
    model.set(GRB_StringAttr_ModelName, "TSPR - TSP with Refueling"); // name to the problem
    model.set(GRB_IntAttr_ModelSense, GRB_MINIMIZE); // is a minimization problem

    // Add one binary variable for each arc and also sets its cost in the objective function
    for (ArcIt e(tsp.g); e!=INVALID; ++e) {
        char name[100];
        sprintf(name,"x_%s_%s",tsp.vname[tsp.g.source(e)].c_str(),tsp.vname[tsp.g.target(e)].c_str());
        x[e] = model.addVar(0.0, 1.0, tsp.weight[e],GRB_BINARY,name);
    }
    model.update(); // run update to use model inserted variables

    // Adicione restrições abaixo

    // (1) Nós terminais devem ser visitados exatamente uma vez
    for (auto v : terminais) {
        GRBLinExpr expr;
        for (InArcIt e(tsp.g,v); e!=INVALID; ++e){
            expr += x[e];
        }
        model.addConstr(expr == 1 );
    }

    // (2) No máximo sai um arco de cada terminal
    for (auto v : terminais) {
        GRBLinExpr expr;
        for (OutArcIt e(tsp.g,v); e!=INVALID; ++e){
            expr += x[e];
        }
        model.addConstr(expr <= 1 );
    }

    // (3) Nó source sempre presente no início do caminho
    GRBLinExpr expr;
    for (OutArcIt e(tsp.g,source); e!=INVALID; ++e){
        expr += x[e];
    }
    model.addConstr(expr >= 1 );

    try {
        model.update(); // Process any pending model modifications.
        if (maxTime >= 0) model.getEnv().set(GRB_DoubleParam_TimeLimit,maxTime);

        subtourelim cb = subtourelim(tsp , x, source, postos, delta);
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

        // Calculo manual do custo da solução (deve ser igual ao ObjVal do Gurobi).
        double soma=0.0;
        vector<Arc> edgesSol;
        ArcName aname(tsp.g);
        ArcColorMap acolor(tsp.g);
        if( verbose) cout << "####### " << endl << "Edges of Solution (B&B):" << endl;
        for (ArcIt e(tsp.g); e!=INVALID; ++e){
            if (BinaryIsOne(x[e].get(GRB_DoubleAttr_X))){ // Note que se este método serve para variáveis binárias, p/ inteiras terá de usar outro método.
                soma += tsp.weight[e];
                edgesSol.push_back(e);
                if( verbose) cout << "(" << tsp.vname[tsp.g.source(e)] << "," << tsp.vname[tsp.g.target(e)] << ")" << endl;
                acolor[e] = BLUE;
            }
        }
        if( verbose ) cout << "####### " << endl;

        if( verbose ) cout << "Custo calculado pelo B&B = "<< soma << " / " << custo_solucao << endl;
        sol = path_search(tsp, source, edgesSol);
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
