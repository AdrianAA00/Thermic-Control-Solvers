#include "Solvers.h"

int StableStationarySolver(VectorXd& T, double error_T, int maxIter_T)
{
    int size;                         //Number on nodes
    int boundaryNodes;                //Number of nodes with boundary condition

    SparseMatrix<double> kl_S;        //Vectors and matrices incluiding boundary conditions
    SparseMatrix<double> kr_S;
    VectorXd T0;
    VectorXd QL;
    VectorXd c;
    double time = 0;

    //*****************************************************************************************************************************************************
    //                                                           LU Factorizing kl 
    //*****************************************************************************************************************************************************
    VectorXd b;                                              //Independent term in the fix point iterating method
    SparseMatrix<double> kt;

    //*****************************************************************************************************************************************************
    //                                  Iterating for solving KL*incre_T = F*incre_T: fix point iteration method
    //*****************************************************************************************************************************************************
    MatrixXd T_diago;
    MatrixXd T_diago_3;
    SparseMatrix<double> T_diago_3_S;
    VectorXd incre_T;
    VectorXd incre_T0;
    VectorXd T_4;

    int count_T = 0;

    do
    {
        if (count_T == 0)
        {
            ObjectsDefinition(size, boundaryNodes, kl_S, kr_S, T0, QL, c, T, time);

            SparseMatrix<double> Aux;
            Aux = kl_S.selfadjointView<Upper>();
            kl_S = Aux;
            Aux = kr_S.selfadjointView<Upper>();
            kr_S = Aux;
        }

        if (count_T == 0)
        {
            // T0.setConstant(1);               //Do not start at 0 Kelvin 
            T = T0;                             //Initial iteration point
        }

        T_diago = T.asDiagonal();
        T_diago_3 = T_diago.array().pow(3);
        T_diago_3_S = T_diago_3.sparseView();
        T_diago_3.resize(0, 0);

        T_4 = T.array().pow(4);

        b = -(kl_S * T + QL + kr_S * T_4);
        kt = kl_S + 4 * kr_S * T_diago_3_S;            // solve lu in each iteration

        //Eigen solver
        
        //BiCGSTAB <SparseMatrix<double>> solver;
        //SparseLU
        //SimplicialLDLT  -> Much faster
        //BiCGSTAB
        //ConjugateGradient 
        
        //MKL solver

        Eigen::PardisoLU< SparseMatrix<double>> solver;
        //PardisoLU
        //PardisoLDLT 
        //PardisoLLT

        solver.compute(kt);

        incre_T = solver.solve(b);                                                               //solving system. kl was previously lu factorized

        count_T++;
        T0 = T;                                                                                  //save previous iteration value for stopping criterion
        T = T + incre_T;

        //std::cout << "\n Error \n" << (T - T0).norm() << ", number of iterations T = " << count_T << "\n";

    } while ((T - T0).norm() > error_T && count_T < maxIter_T);                                  //stop criterion for temperature

    return 0;
}
