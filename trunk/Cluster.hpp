#ifndef CLUSTER_HPP
#define CLUSTER_HPP

#include <vector>
#include <ANN/ANN.h>
#include <eigen3/Eigen/Dense>
using namespace Eigen;
struct ClusInfo
{
	ClusInfo(unsigned int cid = 0) : id(cid) {};
	unsigned int id;
	std::vector<unsigned int> mem;
};
class Cluster
{
	public:
		Cluster(class Matrix<float, Dynamic, Dynamic> & cords);
		Cluster(class CompMgr* cmg, float d);
		void GetMembers(unsigned int id, std::vector<unsigned int> & mem, float d);
		std::vector<unsigned int> & GetClusters(float d);
		~Cluster();
		std::vector<unsigned int> clidarr;
	private:
		ANNpointArray datapts;
		ANNkd_tree* kdTree;
		unsigned int nPts;
		unsigned int dim;
		std::vector<ClusInfo> clusters;
};
#endif
