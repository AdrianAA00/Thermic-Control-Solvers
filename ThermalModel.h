#pragma once

#define EIGEN_USE_MKL_ALL
#include "Eigen/Sparse"


struct Thermal_Model
{
	int n;
	int nb;
	double time;
	Eigen::VectorXd T;
	Eigen::VectorXd Q;
	Eigen::VectorXd C;
	Eigen::SparseMatrix<double> kl;
	Eigen::SparseMatrix<double> kr;

};

void read_nodes(std::string fpath, int& n, int& nb, Eigen::VectorXd& T, Eigen::VectorXd& Q, Eigen::VectorXd& C);
void read_couplings(std::string fpath, Eigen::SparseMatrix<double>& k);
void read_thermal_model(std::string folder, Thermal_Model& model); 