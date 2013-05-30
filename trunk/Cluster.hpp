#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <vector>
#include <ANN/ANN.h>
class Cluster
{
	public:
		Cluster(std::vector<std::vector<double> > & surfcords);
		void GetClusters();
		~Cluster();
	private:
		ANNpointArray datapts;
		ANNkd_tree* kdTree;
		unsigned int nPts;
		unsigned int dim;
	//	std::vector<std::vector<double> > & datapts;
};
#endif
