#include "Solvers.h"

int StationarySolver(VectorXd& T, double error_T, int maxIter_T, double error_fix_point, int maxIter_fix_point)
{
    int size;                         //Number on nodes
    int boundaryNodes;                //Number of nodes with boundary condition

    SparseMatrix<double> kl_S;        //Vectors and matrices incluiding boundary conditions
    SparseMatrix<double> kr_S;
    VectorXd T0;
    VectorXd QL;
    VectorXd c;
    double time = 0;

    //**********************************************************************************************************************************************************
    //                             LU Factorizing kl and direct solving (no fix point method)  Faster vs Non-Stable
    //**********************************************************************************************************************************************************
    MatrixXd T_diago;
    MatrixXd T_diago_3;
    SparseMatrix<double> T_diago_3_S;
    VectorXd incre_T;
    VectorXd incre_T0;
    VectorXd T_4;

    int count_fix_point;
    int count_T = 0;

    do
    {
        if (count_T == 0)
        { 
            ObjectsDefinition(size, boundaryNodes, kl_S, kr_S, T0, QL, c, T, time);
            T_diago_3_S = kl_S.selfadjointView<Upper>();
            kl_S = T_diago_3_S;

            T_diago_3_S = kr_S.selfadjointView<Upper>();
            kr_S = T_diago_3_S;
        }
        

        T0.setConstant(0);;                    //Do not start at 0 Kelvin 

        if (count_T == 0)
        {
            T = T0;                             //Initial iteration point
        }

        // LU Factorizing
        VectorXd b;                                                  //Independent term in the fix point iterating method

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

        solver.compute(kl_S);
    
        T_diago = T.asDiagonal();
        T_diago_3 = T_diago.array().pow(3);
        T_diago_3_S = T_diago_3.sparseView();
        T_diago_3.resize(0, 0);

        T_4 = T.array().pow(4);
        incre_T = VectorXd::Zero(size);                  // Start point of unknown variable we are calculating
        count_fix_point = 0;                             // Record number of iterations  

        do
        {
            b = -(kl_S * T + QL + kr_S * T_4 + 4 * kr_S * (T_diago_3_S * incre_T));
            incre_T0 = incre_T;                                                                                 //Save previous iteration value for stopping criterion
            incre_T = solver.solve(b);                                                                          //Solving system. kl was previously LU factorized
            count_fix_point++;
        } while ((incre_T - incre_T0).norm() > error_fix_point && count_fix_point < maxIter_fix_point);         //Stop criterion for fix point iteration method

        std::cout << "\n incre_T = \n" << ", number of iterations incre_T = " << count_fix_point << "\n";
        count_T++;
        T0 = T;                                                                                  //Save previous iteration value for stopping criterion
        T = T + incre_T;
        std::cout << "\n Error = \n" << incre_T.norm() << ", number of iterations T = " << count_T << "\n";

    } while ((T - T0).norm() > error_T && count_T < maxIter_T);                                  //Stop criterion for temperature

    return 0;
}
