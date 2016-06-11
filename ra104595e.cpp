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
        try {
            // ******* ALTERE DAQUI PARA BAIXO (BRANCH AND CUT) ********** //
            // Para obter os valores da solução fracionária:
            for (ArcIt e(tsp.g); e != INVALID; ++e) {
                (this->*solution_value)(x[e]); // valor fracionário da variável correspondente à aresta 'e'.
            }

            // Para inserir os cortes utilize addLazy().

            // ******* ALTERE DAQUI PARA CIMA (BRANCH AND CUT) ********** //
        }
        catch (...) {
            cout << "Error during callback..." << endl;
        }
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
        lbound = model.get(GRB_DoubleAttr_ObjBoundC);

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