#include "Cluster.hpp"
#include <stdio.h>

Cluster::Cluster(Matrix<float, Dynamic, Dynamic> & cords)
{
	nPts = cords.rows();
	if(nPts)
	{
		dim = cords.cols();
		datapts = annAllocPts(nPts, dim);

		for(unsigned int i = 0; i < nPts; i++)
		{
			for(unsigned int d = 0; d < dim; d++)
				datapts[i][d] = cords(i,d);
		}

		kdTree = new ANNkd_tree(datapts, nPts, dim);	
	}
}

void Cluster::GetCluster(unsigned int id)
{
	unsigned int sz = kdTree->annkFRSearch(datapts[id], 1, 0);
	printf("For pt %d sphere contains %d\n", id, sz);
	ANNidxArray nnIdx = new ANNidx[sz];						// allocate near neigh indices
	ANNdistArray dists = new ANNdist[sz];						// allocate near neighbor dists
	kdTree->annkFRSearch(datapts[id], 1, sz, nnIdx, dists);

	for (unsigned int j = 0; j < sz; j++) 
	{
		std::cout << "\t" << j << "\t" << nnIdx[j] << "\t" << dists[j] << "\n";
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
