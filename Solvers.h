#pragma once


#define EIGEN_USE_MKL_ALL


#include <omp.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <chrono>
#include <iostream>
#include <Eigen/PardisoSupport>


using namespace std::chrono;
using namespace Eigen;
using namespace std;

//Definition and initialization of objects needed during the program
int ObjectsDefinition(int& size, int& boundaryNodes, SparseMatrix<double>& kl_S, SparseMatrix<double>& kr_S, VectorXd& T0, VectorXd& QL, VectorXd& c, VectorXd T, double time);
int BoundaryNodesRemover(int boundaryNodes, int& size, MatrixXd& kl, MatrixXd& kr, VectorXd& T0, VectorXd& QL, VectorXd& c);
int BoundaryNodesRemoverSparse(int boundaryNodes, int& size, SparseMatrix<double>& kl, SparseMatrix<double>& kr, VectorXd& T0, VectorXd& QL, VectorXd& c);
int Size();

//Solvers
int Temporalsolver(string typeSolver, double t_Incre, int maxIter, MatrixXd& T);
int StableStationarySolver(VectorXd& T, double error_T = 10E-3, int maxIter_T = 10E2);                                                             // Set Error and max iterations of the execution
int StationarySolver(VectorXd& T, double error_T = 10E-3, int maxIter_T = 10E4, double error_fix_point = 10E-8, int maxIter_fix_point = 10E1);     // Set Error and max iterations of the execution of fix iteration method and in temperature
int StableStationarySolverInternal(SparseMatrix<double> kl_S, SparseMatrix<double> kr_S, VectorXd QL, VectorXd& T, double error_T = 10E-6, int maxIter_T = 10E3);    // Set Error and max iterations of the execution


