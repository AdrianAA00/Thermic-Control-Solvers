// C++ 2017
// Author: Adrián Antón Álvarez & IDR

//**********************************************Description********************************************************
// Solver for non linear systems of equation. Used for thermic control problems.
// This code solves the problems of stability for certain values of the matrices.
// Sample matrices are used to demonstrate how it works.
// Eigen library needed.
//*****************************************************************************************************************

#include <iostream>
#include <Eigen/Dense>
#include <chrono>

using namespace std::chrono;
using namespace Eigen;

int main()
{
    int size = 4;                    //Number on nodes

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
    VectorXd T0;

    QL << 50, 100, 75, 75;  // External thermal loads
    QL = QL * 10E-5;

    std::cout << "kl" << "\n" << kl << "\n";
    std::cout << "\n" << "kr" << "\n" << kr << "\n";

    // Get starting timepoint
    auto start = high_resolution_clock::now();


    //****************************************************************************************************************
    // LU Factorizing kl 
    VectorXd b(size);                                               //Independent term in the fix point iterating method
    MatrixXd kt(size, size);
    MatrixXd kt_LU(size, size);                                      //Save kl for future operations



    //****************************************************************************************************************
    //Iterating for solving KL*incre_T = F*incre_T: fix point iteration method

    MatrixXd T_diago(size, size);
    MatrixXd T_diago_3(size, size);
    VectorXd T_4 = T.array().pow(4);
    VectorXd incre_T(size);
    VectorXd incre_T0(size);

    int count_T = 0;

    do
    {
        T_diago = T.asDiagonal();
        T_diago_3 = T_diago.array().pow(3);
        T_4 = T.array().pow(4);
        incre_T = VectorXd::Zero(size);         // Start point of unknown variable we are calculating

        b = -(kl * T + QL + kr * T_4);

        kt = kl + 4 * kr * T_diago_3;            //solve lu in each iteration
        kt_LU = kt;
        PartialPivLU<Ref<MatrixXd> > lu(kt_LU);

        incre_T = lu.solve(b);                                                                   //solving system. kl was previously lu factorized

        //std::cout << "\n incre_T = \n" << incre_T;

        count_T++;
        T0 = T;                                                                                  //save previous iteration value for stopping criterion
        T = T + incre_T;
        //std::cout << "\n T = \n" << T << ", number of iterations T = " << count_T << "\n";

    } while ((T - T0).norm() > error_T && count_T < maxIter_T);                                  //stop criterion for temperature

    std::cout << "\n T:\n" << T;


    // Get ending timepoint
    auto stop = high_resolution_clock::now();

    // Get duration of execution.
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "\n Time taken by function: \n " << duration.count() << " microseconds";


    return 0;
}
