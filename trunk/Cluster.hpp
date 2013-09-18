#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <vector>
#include <ANN/ANN.h>
#include <eigen3/Eigen/Dense>
using namespace Eigen;
class Cluster
{
	public:
		Cluster(class Matrix<float, Dynamic, Dynamic> & cords);
		void GetCluster(unsigned int);
		~Cluster();
	private:
		ANNpointArray datapts;
		ANNkd_tree* kdTree;
		unsigned int nPts;
		unsigned int dim;
	//	std::vector<std::vector<double> > & datapts;
};
#endif
