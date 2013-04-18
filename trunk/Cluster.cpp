#include "Cluster.hpp"
#include <stdio.h>
Cluster::Cluster(std::vector<std::vector<double> > & surfcords) 
{
	nPts = surfcords.size();
	if(nPts)
	{
		dim = surfcords[0].size();
		datapts = annAllocPts(nPts, dim);

		for(unsigned int i = 0; i < nPts; i++)
		{
			for(unsigned int d = 0; d < dim; d++)
				datapts[i][d] = surfcords[i][d]/surfcords[i][dim-1];
		}

		kdTree = new ANNkd_tree(datapts, nPts, dim);	
	}
}

void Cluster::GetClusters()
{
	for(unsigned int i = 0; i < nPts; i++)
	{
		unsigned int sz = kdTree->annkFRSearch(datapts[i], 1.2, 0);
		printf("For pt %d sphere contains %d\n", i, sz);
		ANNidxArray nnIdx = new ANNidx[sz];						// allocate near neigh indices
		ANNdistArray dists = new ANNdist[sz];						// allocate near neighbor dists
		kdTree->annkFRSearch(datapts[i], 1.2, sz, nnIdx, dists);

		for (unsigned int j = 0; j < sz; j++) 
		{
			std::cout << "\t" << j << "\t" << nnIdx[j] << "\t" << dists[j] << "\n";
		}
		delete [] nnIdx;
		delete [] dists;
	}

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
