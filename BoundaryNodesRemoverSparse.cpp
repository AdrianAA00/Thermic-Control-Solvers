#include "Solvers.h"

int BoundaryNodesRemoverSparse(int boundaryNodes, int& size, SparseMatrix<double>& kl, SparseMatrix<double>& kr, VectorXd& T0, VectorXd& QL, VectorXd& c)
{
    // ****************************************************************************************************************************************************
    //           Converts hiperstatic system with boundary conditions to a non-singular system removing nodes with a fixed temperature
    // ****************************************************************************************************************************************************
    SparseMatrix<double> k_short;
    SparseMatrix<double> k_;
    VectorXd T0_short;
    VectorXd T_4;
    size = kl.rows();

    k_ = kl.selfadjointView<Upper>();
    k_short = k_.topRightCorner(size, boundaryNodes);
    T0_short = T0.bottomRows(boundaryNodes);
    QL += k_short * T0_short;

    T_4 = T0_short.array().pow(4);
    k_ = kr.selfadjointView<Upper>();
    k_short = k_.topRightCorner(size, boundaryNodes);
    QL += k_short * T_4;

    kl.conservativeResize(size - boundaryNodes, size - boundaryNodes);
    kr.conservativeResize(size - boundaryNodes, size - boundaryNodes);
    T0.conservativeResize(size - boundaryNodes);
    QL.conservativeResize(size - boundaryNodes);
    c.conservativeResize(size - boundaryNodes);
    size -= boundaryNodes;

    return 0;
}