#include "Solvers.h"

int BoundaryNodesRemover(int boundaryNodes, int& size, MatrixXd& kl, MatrixXd& kr, VectorXd& T0, VectorXd& QL, VectorXd& c)
{
    // ****************************************************************************************************************************************************
    //           Converts hiperstatic system with boundary conditions to a non-singular system removing nodes with a fixed temperature
    // ****************************************************************************************************************************************************
    MatrixXd k_short;
    VectorXd T0_short;
    VectorXd T_4;
    size = kl.rows();

    k_short = kl.topRightCorner(size, boundaryNodes).selfadjointView<Upper>();
    T0_short = T0.bottomRows(boundaryNodes);
    QL += k_short * T0_short;

    T_4 = T0_short.array().pow(4);
    k_short = kr.topRightCorner(size, boundaryNodes).selfadjointView<Upper>();
    QL += k_short * T_4;

    kl.conservativeResize(size - boundaryNodes, size - boundaryNodes);
    kr.conservativeResize(size - boundaryNodes, size - boundaryNodes);
    T0.conservativeResize(size - boundaryNodes);
    QL.conservativeResize(size - boundaryNodes);
    c.conservativeResize(size - boundaryNodes);
    size -= boundaryNodes;

    return 0;
}