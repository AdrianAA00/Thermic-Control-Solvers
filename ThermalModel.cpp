#include "Eigen/Sparse"
#include <fstream>
#include <iostream>

#include "ThermalModel.h" 


void read_nodes(std::string fpath, int& n, int& nb, Eigen::VectorXd& T, Eigen::VectorXd& Q, Eigen::VectorXd& C)
{
	std::ifstream read_file(fpath);
	assert(read_file.is_open());
	char node_type;

	nb = 0;

	T.resize(n);
	Q.resize(n);
	C.resize(n);

	for (int i = 0; i < n; i++)
	{
		//Read the data
		read_file >> node_type >> C[i] >> Q[i] >> T[i];

		// Number of boundary nodes
		if (node_type == 'B') {
			nb++;
		}
	}
	read_file.close();

	return;
}

void read_couplings(std::string fpath, Eigen::SparseMatrix<double>& k)
{
	std::ifstream read_file(fpath);
	assert(read_file.is_open());

	std::string line;

	int i = 0;
	int j = 0;
	double coupling_val = 0;

	typedef Eigen::Triplet<double> Trip;
	std::vector<Trip> trp;

	while (!read_file.eof())//(!read_file.eof())
	{
		// Create the triplets
		// (index, index, value)
		read_file >> i >> j >> coupling_val;
		
		trp.push_back(Trip(i, j, coupling_val));

		//Values of the diagonal;
		trp.push_back(Trip(i, i, -coupling_val));
		trp.push_back(Trip(j, j, -coupling_val));
	}

	//Assign them to the sparse Eigen matrix
	k.setFromTriplets(trp.begin(), trp.end());

	read_file.close();

	return;
}

void read_thermal_model(std::string folder, Thermal_Model& model) {

	//Read nodes
	read_nodes("Nodes.dat", model.n, model.nb, model.T, model.Q, model.C);

	//Resize coupling matrices
	model.kl.resize(model.n, model.n);
	model.kr.resize(model.n, model.n);

	//Read linear couplings
	read_couplings("GLs.dat", model.kl);

	//Read radiative couplings
	read_couplings("GRxBB.txt", model.kr);
	model.kr = 5.670E-8 * model.kr;
}