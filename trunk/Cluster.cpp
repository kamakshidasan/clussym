#include "Cluster.hpp"
#include <stdio.h>

Cluster::Cluster(Matrix<float, Dynamic, Dynamic> & cords) : clidarr(cords.rows(),0)
{
	clusters.push_back(ClusInfo());
	nPts = cords.rows();
	if(nPts)
	{
		dim = cords.cols();
		datapts = annAllocPts(nPts, dim);

		for(unsigned int i = 0; i < nPts; i++)
		{
			for(unsigned int d = 0; d < dim; d++)
				datapts[i][d] = cords(i,d);
			std::cout<<"Cord "<<i<<" : "<<std::endl<<cords.row(i)<<std::endl;
		}

		kdTree = new ANNkd_tree(datapts, nPts, dim);	
	}
}
std::vector<unsigned int> & Cluster::GetClusters()
{
	for(unsigned int i = 0; i < nPts; i++)
	{
		if(clidarr[i] == 0)
		{
			ClusInfo clinfo;
			clinfo.id = clusters.size();
			GetMembers(i, clinfo.mem);
			clusters.push_back(clinfo);
		}
	}
	return clidarr;
}
void Cluster::GetMembers(unsigned int id, std::vector<unsigned int> & mem)
{
	unsigned int sz = kdTree->annkFRSearch(datapts[id], 0.1, 0);
	printf("For pt %d sphere contains %d\n", id, sz);
	ANNidxArray nnIdx = new ANNidx[sz];						// allocate near neigh indices
	ANNdistArray dists = new ANNdist[sz];						// allocate near neighbor dists
	kdTree->annkFRSearch(datapts[id], 0.1, sz, nnIdx, dists);

	for (unsigned int j = 0; j < sz; j++) 
	{
		std::cout << "\t" << j << "\t" << nnIdx[j] << "\t" << dists[j] << "\n";
		clidarr[nnIdx[j]] = clusters.size();
	}
	delete [] nnIdx;
	delete [] dists;

}

Cluster::~Cluster()
{
	if(nPts)
	{
		delete kdTree;
		annDeallocPts(datapts);
		annClose();
	}	
}
