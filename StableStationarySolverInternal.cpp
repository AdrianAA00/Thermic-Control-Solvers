#include "Solvers.h"

int StableStationarySolverInternal(SparseMatrix<double> kl_S, SparseMatrix<double> kr_S, VectorXd QL, VectorXd& T, double error_T, int maxIter_T)
{
    //*****************************************************************************************************************************************************
    //                                                           LU Factorizing kl 
    //*****************************************************************************************************************************************************
    VectorXd b;                                              //Independent term in the fix point iterating method
    SparseMatrix<double> kt;

    //*****************************************************************************************************************************************************
    //                                  Iterating for solving KL*incre_T = F*incre_T: fix point iteration method
    //*****************************************************************************************************************************************************
    VectorXd T0;
    MatrixXd T_diago;
    MatrixXd T_diago_3;
    SparseMatrix<double> T_diago_3_S;
    VectorXd incre_T;
    VectorXd incre_T0;
    VectorXd T_4;

    int count_T = 0;

    do
    {
        T_diago = T.asDiagonal();                      // Start iterating with previous temperature for faster convergence
        T_diago_3 = T_diago.array().pow(3);
        T_diago_3_S = T_diago_3.sparseView();
        T_diago_3.resize(0, 0);

        T_4 = T.array().pow(4);

        b = -(kl_S * T + QL + kr_S * T_4);
        kt = kl_S + 4 * kr_S * T_diago_3_S;            // solve lu in each iteration

        //-> Eigen solver
        //BiCGSTAB <SparseMatrix<double>> solver;
        //SparseLU
        //SimplicialLDLT  -> Much faster
        //BiCGSTAB
        //ConjugateGradient 

        //-> MKL solver
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
