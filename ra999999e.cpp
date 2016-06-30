// MC658 - Projeto e Análise de Algoritmos III
// ---- LABORATÓRIO TSPR ----
// Prof: Flavio Keidi Miyazawa
// PED: Rafael Arakaki

#include <iostream>
#include <float.h>
#include <gurobi_c++.h>
#include <lemon/list_graph.h>
#include "mygraphlib.h"
#include "algoritmos.h"

using namespace lemon;
using namespace std;

// Rotina de callback para inserir cortes durante a execução do B&B
class subtourelim: public GRBCallback
{
    TSP_Data_R &tsp;
    ListDigraph::ArcMap<GRBVar>& x;

    double (GRBCallback::*solution_value)(GRBVar);
public:
    subtourelim(TSP_Data_R &tsp, ListDigraph::ArcMap<GRBVar>& x) : tsp(tsp),x(x)
    {    }
protected:
    void callback()
    {
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
        // -----------------------------------------------------------------
        // Guarda os arcos com valor fracionario e inteiro
        vector<Edges> FracArcs, OneArcs;

        // cria um subgrafo h com arestas com x[e]==1
        // assim pode-se aplicar Gomory-Hu num grafo pequeno
        ListDigraph::ArcMap<bool> one_filter(tsp.g, false);    // inicia sem nenhuma aresta
        ListDigraph::ArcMap<bool> non_zero_filter(tsp.g, false); // inicia sem nenhuma aresta

        for (ArcIt a(tsp.g); a!=INVALID; ++a) {
            if ((this->*solution_value)(x[a]) > 1-MY_EPS)
                OneArcs.push_back(a);    // guarda arcos com x[a]==1
            else if((this->*solution)(x[a]) > MY_EPS) 
                FracArcs.push_back(a);  // inclui arcos com valor 0 < x[a] < 1
        } // define o subgrafo com arcos tais que x[a]==1

        try {
            // Para obter os valores da solução fracionária:
            // for (ArcIt e(tsp.g); e != INVALID; ++e) {
            //     (this->*solution_value)(x[e]); // valor fracionário da variável correspondente à aresta 'e'.
            // }

            // Para inserir os cortes utilize addLazy().

            // -----------------------------------------------------------------
            // Utiliza union-find para contrair nos (obtem um grafo onde cada componente de g e contraida)
            ListDigraph::NodeMap<int> aux_map(tsp.g);
            UnionFind<ListDigraph::NodeMap<int> > UFNodes(auxmap);
            for (NodeIt v(tsp.g); v!=INVALID; ++v) { 
                UFNodes.insert(v);
            }
            for (vector<Arc>::iterator a_it=OneArcs.begin(); a_it!=OneArcs.end(); ++a_it) {
                UFNodes.join(tsp.g.u(*a_it), tsp.g.v(*a_it));   // Sem problema se sao a mesma componente
            }
            // -----------------------------------------------------------------
            // Coloca num conjunto separado todas as arestas que nao estao numa componente
            vector<Arc> CrossingArcs;
            for (ArcIt a(tsp.g); a!=INVALID; ++a) {
                if (UFNodes.find(tsp.g.u(a)) != UFNodes.find(tsp.g.v(e)))
                    CrossingArcs.push_back(a);
            }

            // -----------------------------------------------------------------
            // Gera uma lista invertida UFIndexToNode para encontrar o no que representa a componente
            vector<bool> ComponentIndex(tsp.NNodes);
            vector<DNode> Index2h(tsp.NNodes);
            for (int=0;i<tsp.NNodes;i++) {
                ComponentIndex[i]=false;
            }
            for (NodeIt v(tsp.g); v!=INVALID; ++v) {
                ComponentIndex[UFNodes.find(v)]=true;
            }
            // -----------------------------------------------------------------
            // Gera um grafo de componentes, adiciona um no de cada componente e arcos
            ListDigraph h;
            ArcValueMap h_capacity(h);

            for (int i=0;i<tsp.NNodes;++i) {
                if (ComponentIndex[i]) {
                    Index2h[i]=h.addNode();
                }
            }
            for (vector<Arc>::iterator a_it=FracArcs.begin(); a_it != FracArcs.end(); ++a_it) {
                    DNode u = tsp.g.u(*a_it);
                    DNode v = tsp.g.v(*a_it);
                    DNode hu = Index2h[UFNodes.find(u)];
                    DNode hv = Index2h[UFNodes.find(v)];
                    Arc a = h.addArc(hu, hv);   // Adiciona arcos para o grafo h
                    h_capacity[a] = (this->*solution_value)(x[*a_it]);
            }
            // -----------------------------------------------------------------
            GomoryHu<ListDigraph, ArcValueMap> ght(h, h_capacity);  ght.run();
            // A arvore de Gomory-Hu e dada como um grafo direcionado com raiz, Cada no tem
            // um arco que aponta para o seu pai. A raiz tem pai -1.
            // Lembre que cada arco nesta arvore representa um corte e o valor do
            // arco eh o peso do corte correspondente. Entao, se um arco tem peso
            // menor do que 2, entao achamos um corte violado e nesse caso, inserimos
            // a restricao correspondente
            
            NodeBoolMap cutmap(h);
            for (NodeIt u(h); u!=INVALID; ++u) {
                GRBLinExpr expr = 0;
                if (ght.predNode(u)==INVALID) continue; // pula a raiz
                if (ght.predValue(u) > 2.0 - MY_EPS) continue; // valor de corte permitido
                ght.minCutMap(u, ght.predNode(u), cutmap);  // agora temos um corte que viola

                // Percorre as arestas que cruzam alguma componente e insera as que pertencem ao corte
                for (vector<Arc>::iterator a_it=CrossingArcs.begin(); a_it!=CrossingArcs.end(); ++a_it) {
                    DNode u=tsp.g.u(*a_it), v=tsp.g.v(*a_it);
                    DNode hu=Index2h[UFNodes.find(u)], hv=Index2h[UFNodes.find(v)];

                    if (cutmap[hu]!=cutmap[hv]) {
                        expr += x[a_it];
                    }
                }
                addLazy( expr >= 2);
            }
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

        subtourelim cb = subtourelim(tsp , x);
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
