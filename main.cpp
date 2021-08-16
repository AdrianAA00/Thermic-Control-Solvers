/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                 //
//                                                                          C++ 2017                                                                               //
//                                                                                                                                                                 //
//                                                                Author: Adrián Antón Álvarez & IDR                                                               //                    
//                                                                                                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//************************************************************************ DESCRIPTION ****************************************************************************//
//  1) Eigen library needed -> 3.4_rc1 version or more recent                                                                                                      //
//  2) Definition and inicialization of objects needed during the program                                                                                          //
//      2.1.) BoundaryNodesRemover -> Converts hiperstatic system with boundary conditions to a non-singular system removing nodes with a fixed temperature        //
//      2.2.) BoundaryNodesRemoverSparse -> Does the same as BoundaryNodesRemover but with sparse matrices                                                         //
//      2.3.) ObjectsDefinition -> Sample function which provides objects needed during execution                                                                  //
//      2.4.)  Size -> Retrieves size of matrices without boundary conditions                                                                                      //
//                                                                                                                                                                 //
//  3) Solvers:                                                                                                                                                    //
//      3.1.) StationarySolver -> Static solver: faster but non stable. Uses fix point iteration                                                                   //
//      3.2.) StableStationarySolver -> Static solver: slower but more robust. Uses direct LU factorization in each iteration                                      //  
//      3.3.) StableStationarySolverInternal -> Static solver used for the implicit temporal schemes                                                               //                                                                                                       
//      3.3.) Temporalsolver -> Dynamic solver. Adjustable time step and range of time. Three solvers:                                                             //
//              3.3.1.) Euler                                                                                                                                      //
//              3.3.2.) AB2 -> Adams Bashforth 2                                                                                                                   //
//              3.3.3.) RK4 -> Runge Kutta 4                                                                                                                       //
//                                                                                                                                                                 //
// 4) Thermal Model:                                                                                                                                               //
//      4.1.) ThermalModel -> Reads Gls, Grs, T0, QL and C and creates the matrices for solving the problem                                                        //
//                                                                                                                                                                 //
// 5) Header Files:                                                                                                                                                //
//      5.1.) Solvers                                                                                                                                              //
//      5.2.) ThermalModel                                                                                                                                         //
//                                                                                                                                                                 //
//*****************************************************************************************************************************************************************//
#include "Solvers.h"

int main()
{
    Eigen::initParallel();

    // Get starting timepoint
    auto start = high_resolution_clock::now();

    // ****************************************************************************************************************************************************
    //                                                                SOLVING EXAMPLE
    // ****************************************************************************************************************************************************
    std::cout << "\n" << "---------------------------------------------------------------------------------------------" << "\n";
    std::cout << "\n" << "|                            Author: Adrian Anton Alvarez & IDR                             |" << "\n";
    std::cout << "\n" << "---------------------------------------------------------------------------------------------" << "\n";


    std::cout << "\n" << "---------------------------------------------------------------------------------------------" << "\n";
    std::cout << "\n" << "|                                      SOLVING EXAMPLE                                      |" << "\n";
    std::cout << "\n" << "---------------------------------------------------------------------------------------------" << "\n";

    // ----------------Temporal Solver-----------------------
    // Parameters of the solver
    string typeSolver;

    double t_Incre = 10E-2;
    int maxIter = 10;

    // Set the size of the matrix without boundary conditions
    int size = Size();

    //Matrix used for storing the evolution of temperature
    MatrixXd T(size, maxIter);

    // ------------------Temporal Solvers------------------------
    //      typeSolver:
    //      CN -> Crank Nicolson Implicit   ->  More stable / Slower
    //      RK4 -> Runge Kutta 4 Explicit 4th orden
    //      AB2 -> Adam Bashford 2nd order Explicit
    //      Euler -> Euler 1st order Explicit

    typeSolver = "Euler"; 
    Temporalsolver(typeSolver, t_Incre, maxIter, T);
    std::cout << "\n" << "RK4:" << "\n" << T.col(maxIter - 1) << "\n";


    //// ----------------Stationary Solver-----------------------
    //VectorXd T_StatS;
    //StableStationarySolver(T_StatS);
    //std::cout << "\n" << " StableStationarySolver:" << "\n" << T_StatS << "\n";

    //VectorXd T_Stat;
    //StationarySolver(T_Stat);  //Highly unstable
    //std::cout << "\n" << "StationarySolver:" << "\n" << T_Stat << "\n";

    // Get ending timepoint
    auto stop = high_resolution_clock::now();

    // Get duration of execution.
    auto duration = duration_cast<microseconds>(stop - start);
    std::cout << "\n Time taken by function: \n " << duration.count() << " microseconds";

    return 0;
}
