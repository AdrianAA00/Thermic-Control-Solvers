#include "Solvers.h"

int Size()
{
    int size;                         //Number on nodes
    int boundaryNodes;                //Number of nodes with boundary condition

    SparseMatrix<double> kl_S;        //Vectors and matrices incluiding boundary conditions
    SparseMatrix<double> kr_S;
    VectorXd T0;
    VectorXd QL;
    VectorXd c;
    VectorXd T;
    double time = 0;

    ObjectsDefinition(size, boundaryNodes, kl_S, kr_S, T0, QL, c, T, time);

    return size;
}