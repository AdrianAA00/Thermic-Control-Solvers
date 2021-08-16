#include "Solvers.h"
#include <Eigen/PardisoSupport>

int Temporalsolver(string typeSolver, double t_Incre, int maxIter, MatrixXd& T)
{
    // ****************************************************************************************************************************************************
    //                          EULER               //                     Solved by Euler integration. It is the fastest but less stable method.     
    // ****************************************************************************************************************************************************
    if (typeSolver == "Euler")
    {
        int size;                         //Number on nodes
        int boundaryNodes;                //Number of nodes with boundary condition

        SparseMatrix<double> kl_S;        //Vectors and matrices incluiding boundary conditions
        SparseMatrix<double> kr_S;
        VectorXd T0;
        VectorXd QL;
        VectorXd c;
        double time = 0;

        int count;
        VectorXd T_4;
        VectorXd T_col;
        MatrixXd c_diago;

        count = 0;

        do
        {
            ObjectsDefinition(size, boundaryNodes, kl_S, kr_S, T0, QL, c, T.col(count), time);

            if (count == 0)
            {
                T.col(0) = T0;                   // Initial conditions
            }

            c = c.cwiseInverse();
            c_diago = c.asDiagonal();
            T_col = T.col(count);
            T_4 = T_col.array().pow(4);

            T.col(count + 1) = T.col(count) + t_Incre * c_diago * (kl_S.selfadjointView<Upper>() * T.col(count) + QL + kr_S.selfadjointView<Upper>() * T_4);          //Finite differential equation

            time += t_Incre;
            count++;
            std::cout << "count" << count << "\n";
        } while (count + 1 < maxIter);
    }

    // ****************************************************************************************************************************************************
    //      ADAMS BASHFORTH 2      //      Second order explicit method --- U(n+1) <- U(n) + Dt/2 ( 3 F(n)-F(U(n-1) ) --- x2 times slower than Euler
    // ****************************************************************************************************************************************************
    else if (typeSolver == "AB2")
    {
        int size;                         //Number on nodes
        int boundaryNodes;                //Number of nodes with boundary condition

        SparseMatrix<double> kl_S;        //Vectors and matrices incluiding boundary conditions
        SparseMatrix<double> kr_S;
        VectorXd T0;
        VectorXd QL;
        VectorXd c;
        double time = 0;

        int count;
        VectorXd T_4_0;
        VectorXd T_col_0;
        VectorXd F_0;
        MatrixXd c_diago;

        VectorXd T_4_1;
        VectorXd T_col_1;
        VectorXd F_1;

        count = 0;

        do
        {
            ObjectsDefinition(size, boundaryNodes, kl_S, kr_S, T0, QL, c, T.col(count), time);

            c = c.cwiseInverse(); 
            c_diago = c.asDiagonal();

            if (count == 0)
            {
                T.col(0) = T0;                             // Initial conditions
            }

            //count
            T_col_0 = T.col(count);
            T_4_0 = T_col_0.array().pow(4);
            F_0 = c_diago * (kl_S.selfadjointView<Upper>() * T.col(count) + QL + kr_S.selfadjointView<Upper>() * T_4_0);

            if (count == 0)
            {
                T.col(1) = T0 + t_Incre * F_0;             //Aproximation of second temperature term with Euler
            }

            //count + 1
            T_col_1 = T.col(count + 1);
            T_4_1 = T_col_1.array().pow(4);
            F_1 = c_diago * (kl_S.selfadjointView<Upper>() * T.col(count + 1) + QL + kr_S.selfadjointView<Upper>() * T_4_1);

            T.col(count + 2) = T.col(count + 1) + t_Incre * (3. / 2. * F_1 - 1. / 2. * F_0);          //Finite differences equation

            time += t_Incre;
            count++;
            std::cout << "count" << count << "\n";
        } while (count + 2 < maxIter);

    }

    // ****************************************************************************************************************************************************
    //                                            RUNGE KUTTA 4         // Fourth order explicit method. x4 times slower than Euler
    // ****************************************************************************************************************************************************
    else if (typeSolver == "RK4")
    {

        int size;                         //Number on nodes
        int boundaryNodes;                //Number of nodes with boundary condition

        SparseMatrix<double> kl_S;        //Vectors and matrices incluiding boundary conditions
        SparseMatrix<double> kr_S;
        VectorXd T0;
        VectorXd QL;
        VectorXd c;
        double time = 0;

        int count;
        VectorXd T_4;
        VectorXd T_col;
        MatrixXd c_diago;

        VectorXd k1, k2, k3, k4;

        count = 0;

        do
        {
            ObjectsDefinition(size, boundaryNodes, kl_S, kr_S, T0, QL, c, T.col(count), time);

            if (count == 0)
            {
                T.col(0) = T0;                   // Initial conditions
            }

            c = c.cwiseInverse();
            c_diago = c.asDiagonal();

            // k1
            T_col = T.col(count);
            T_4 = T_col.array().pow(4);
            k1 = c_diago * (kl_S.selfadjointView<Upper>() * T_col + QL + kr_S.selfadjointView<Upper>() * T_4);

            // k2
            T_col = T.col(count) + t_Incre / 2 * k1;
            T_4 = T_col.array().pow(4);
            k2 = c_diago * (kl_S.selfadjointView<Upper>() * T_col + QL + kr_S.selfadjointView<Upper>() * T_4);

            // k3
            T_col = T.col(count) + t_Incre / 2 * k2;
            T_4 = T_col.array().pow(4);
            k3 = c_diago * (kl_S.selfadjointView<Upper>() * T_col + QL + kr_S.selfadjointView<Upper>() * T_4);

            // k4
            T_col = T.col(count) + t_Incre * k3;
            T_4 = T_col.array().pow(4);
            k4 = c_diago * (kl_S.selfadjointView<Upper>() * T_col + QL + kr_S.selfadjointView<Upper>() * T_4);

            T.col(count + 1) = T.col(count) + 1. / 6. * t_Incre * (k1 + 2 * k2 + 2 * k3 + k4);          //Finite differential equation

            time += t_Incre;
            count++;
            std::cout << "count" << count << "\n";
        } while (count + 1 < maxIter);
    }


    else if (typeSolver == "CN")
    {
        int size;                         //Number on nodes
        int boundaryNodes;                //Number of nodes with boundary condition

        SparseMatrix<double> kl_S;          //Vectors and matrices incluiding boundary conditions
        SparseMatrix<double> kr_S;
        VectorXd T0;
        VectorXd QL;
        VectorXd c;
        double time = 0;

        SparseMatrix<double> kle_S;
        SparseMatrix<double> kre_S;
        VectorXd QLe;

        int count;
        VectorXd T_4;
        VectorXd T_col;
        SparseMatrix<double> c_diago;

        count = 0;

        do
        {
            ObjectsDefinition(size, boundaryNodes, kl_S, kr_S, T0, QL, c, T.col(count), time);

            if (count == 0)
            {
                T.col(0) = T0;                    //Initial conditions
            }

            size = T0.size();

            
            SparseMatrix<double> AuxS;
            AuxS = kl_S.selfadjointView<Upper>();
            kl_S = AuxS;
            AuxS = kr_S.selfadjointView<Upper>();
            kr_S = AuxS;
     
            
            MatrixXd Aux;
            c = c.cwiseInverse();
            Aux = c.asDiagonal();
            c_diago = Aux.sparseView();
            

            // QLe
            T_col = T.col(count); 
            T_4 = T_col.array().pow(4);
            QLe = t_Incre / 2. * c_diago * (kl_S * T_col + 2 * QL + kr_S * T_4) + T_col;

           
            //kle
            SparseMatrix<double> I(size, size);
            I.setIdentity();
            kle_S = t_Incre / 2. * c_diago * kl_S - I;
            

            // kre
            kre_S = t_Incre / 2. * c_diago * kr_S;

            StableStationarySolverInternal(kle_S, kre_S, QLe, T_col);
            time += t_Incre;
            count++;

            T.col(count) = T_col;

        } while (count + 1 < maxIter);

    }

    return 0;
}