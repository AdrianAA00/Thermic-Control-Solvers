# Thermic-Control-Solvers

Brief Description: high efficiency solvers developed to deal with thermic control problems from satellites. Below are explained all the files needed for running correctly the solvers. The user should only modify the txt files which define the thermic problem for each specific case: GLs, GRxBB and Nodes. In this example the solvers are solving the model of a flat plate with boundary conditions. The data should be provided in txt files with the exact same format as the given txt files. User should also provide the total number of nodes where it´s indicated in the main file.

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                                                                                                                                                 //
//                                                                          C++ 2017                                                                               //
//                                                                                                                                                                 //
//                           :rocket:                         Author: Adrián Antón Álvarez & IDR                      :rocket:                                     //           
//                                                                                                                                                                 //
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//************************************************************************ DESCRIPTION ****************************************************************************//
 1) Eigen library needed -> 3.4_rc1 version or more recent                                                                                     
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

Hope you find it useful :blush:

