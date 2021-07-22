// C++ 2017
// Author: Adrián Antón Álvarez & IDR

//**********************************************Description********************************************************
// Solver for non linear systems of equation. Used for thermic control problems.
// This code solves the problems with fix point iteration method and with sparse matrices. This allows to manage 
// bigger size of matrices.
// Sample matrices are used to demonstrate how it works.
// Eigen library needed.
//*****************************************************************************************************************
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <chrono>

using namespace std::chrono;
using namespace Eigen;

int main()
{
    int size = 4;                    //Number on nodes

    //error_fix_point should be higher than error_T
    double error_fix_point = 10E-10;                    //Error during solving with fix point method (incre_T)
    int maxIter_fix_point = 10E1;                       //Max iterations fix point method (incre_T)

    double error_T = 10E-8;                             //Error during solving T
    int maxIter_T = 10E1;                               //Max iterations T

    MatrixXd kl(size, size);
    MatrixXd kr(size, size);
    VectorXd T(size);
    VectorXd QL(size);


    //**************************************************************************************************************
    // Sample matrices
    kl << -0.0025, 0, 0.0025, 0,     //Conductive terms
        0., -0.09, 0, 0.09,
        0.0025, 0, -0.0165, 0.01,
        0, 0.09, 0.01, -0.1;

    kr << -0.0055, 0, 0.0055, 0,     //Radiative terms
        0, -0.004, 0.004, 0,
        0.0055, 0.004, -0.0365, 0.027,
        0, 0, 0.027, -0.027;

    kr = kr * 10E-8;

    T << 10, 10, 10, 10;    // 10 Kelvin as start point
    VectorXd T0(size);

    QL << 50, 100, 75, 75;  // External thermal loads
    QL = QL * 10E-4;

    std::cout << "kl" << "\n" << kl << "\n";
    std::cout << "\n" << "kr" << "\n" << kr << "\n";

    // Converting dense matrices into sparse
    typedef Eigen::SparseMatrix<double> SparseMat;
    typedef Eigen::SparseVector<double> VectorMat;

    SparseMat kl_(size, size), kr_(size, size), kl_LU_(size, size);

    kl_ = kl.sparseView();
    kr_ = kr.sparseView();

    //--------------------------------------------------Starts timing------------------------------------------------
    // Get starting timepoint
    auto start = high_resolution_clock::now();

    //****************************************************************************************************************
    // Factorizing kl_ 
    VectorXd b(size);        //Independent term in the fix point iterating method
    kl_LU_ = kl_;     //Save kl for future operations
    SparseLU<SparseMat> LUsolve(kl_LU_);


    //std::cout << "\n Matrix A after decomposition:\n" << kl_LU << "\n";


    //****************************************************************************************************************
    //Iterating for solving KL*incre_T = F*incre_T: fix point iteration method
    MatrixXd T_diago(size, size);
    MatrixXd T_diago_3(size, size);
    VectorXd T_4 = T.array().pow(4);
    VectorXd incre_T(size);
    VectorXd incre_T0(size);

    int count_fix_point;
    int count_T = 0;

    do
    {
        T_diago = T.asDiagonal();
        T_diago_3 = T_diago.array().pow(3);
        T_4 = T.array().pow(4);
        incre_T << 0, 0, 0, 0;                  // Start point of unknown variable we are calculating
        count_fix_point = 0;                    // Record number of iterations  


        do
        {
            b = -(kl * T + QL + kr * T_4 + 4 * kr * T_diago_3 * incre_T);
            incre_T0 = incre_T;                                                                                  //Save previous iteration value for stopping criterion
            incre_T = LUsolve.solve(b);                                                                          //Solving system. kl was previously LU factorized
            count_fix_point++;
        } while ((incre_T - incre_T0).norm() > error_fix_point && count_fix_point < maxIter_fix_point);          //Stop criterion for fix point iteration method

        //std::cout << "\n incre_T = \n" << incre_T << ", number of iterations incre_T = "<< count_fix_point << "\n";

        count_T++;
        T0 = T;                                                                                                  //Save previous iteration value for stopping criterion
        T = T + incre_T;
        //std::cout << "\n T = \n" << T << ", number of iterations T = " << count_T << "\n";

    } while ((T - T0).norm() > error_T && count_T < maxIter_T);                                                  //Stop criterion for temperature

    std::cout << "\n T:\n" << T << "\n";

    // Get ending timepoint
    auto stop = high_resolution_clock::now();
    //--------------------------------------------------Finishes timing------------------------------------------------

    // Get duration of execution.
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "\n Time taken by function: \n " << duration.count() << " microseconds";

    return 0;
}
