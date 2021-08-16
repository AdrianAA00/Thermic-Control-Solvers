#include "Eigen/Sparse"
#include <iostream>
#include "chrono"

#include "Solvers.h"
#include "ThermalModel.h"


int ObjectsDefinition(int& size, int& boundaryNodes, SparseMatrix<double>& kl_S, SparseMatrix<double>& kr_S, VectorXd& T0, VectorXd& QL, VectorXd& c, VectorXd T, double time)
{
    //New model
    Thermal_Model model;

    //******************************************* INDICATE HERE THE TOTAL NUMBER OF NODES OF THE MODEL *****************************************
    model.n = 1574;

    //Folder with the files: 'GLs.dat', 'GRs.dat' and 'Nodes.dat'
    std::string model_folder = "D:\\Uni\\3ºGIA\\Beca colaboración\\Codigo Beca\\DataReading";

    //Read Nodes, GLs y GRs
    read_thermal_model(model_folder, model);

    size = model.n - model.nb;
    boundaryNodes = model.nb;

    kl_S.resize(size, size);            //Vectors and matrices without boundary conditions
    kr_S.resize(size, size);
    T0.resize(size);
    QL.resize(size);
    c.resize(size);

    //*****************************************************************************************************************************************************
    //                                                          Removing boundaryNodes
    // ****************************************************************************************************************************************************

    BoundaryNodesRemoverSparse(model.nb, size, model.kl, model.kr, model.T, model.Q, model.C);      // Function which retrieves vectors and matrices without boundary nodes

    kl_S = model.kl;
    model.kl.resize(0, 0);    //Releasing memory

    kr_S = model.kr;
    model.kr.resize(0, 0);

    T0 = model.T;
    model.T.resize(0);

    QL = model.Q;
    model.Q.resize(0);

    c = model.C;
    model.C.resize(0);

    return 0;
}