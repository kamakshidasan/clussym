#include "Cluster.hpp"
#include "CompMgr.hpp"
#include <stdio.h>

Cluster::Cluster(CompMgr* cmp, float dist) : clidarr(cmp->comps.size(), 0)
{
	clusters.push_back(ClusInfo());
	nPts = cmp->comps.size();
	if(nPts)
	{
		dim = cmp->comps[0]->cords.size()+1;
		datapts = annAllocPts(nPts, dim);

		for(unsigned int i = 0; i < nPts; i++)
		{
			float cordid;
//			std::cout<<"Cord "<<i<<" : ";
			for(unsigned int d = 0; d < dim-1; d++)
			{
				cordid = datapts[i][d] = cmp->comps[i]->cords[d];
//				std::cout<<cordid<<" ";
			}
			cordid = datapts[i][dim-1] = cmp->comps[i]->fnid*dist;
	//		std::cout<<cordid<<std::endl;
		}

		kdTree = new ANNkd_tree(datapts, nPts, dim);	
	}
}
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
//			std::cout<<"Cord "<<i<<" : "<<std::endl<<cords.row(i)<<std::endl;
		}

		kdTree = new ANNkd_tree(datapts, nPts, dim);	
	}
}
std::vector<unsigned int> & Cluster::GetClusters(float d)
{
	for(unsigned int i = 0; i < nPts; i++)
	{
		if(clidarr[i] == 0)
		{
			ClusInfo clinfo;
			clinfo.id = clusters.size();
			GetMembers(i, clinfo.mem, d);
			clusters.push_back(clinfo);
		}
	}
	return clidarr;
}
void Cluster::GetMembers(unsigned int id, std::vector<unsigned int> & mem, float d)
{
	unsigned int sz = kdTree->annkFRSearch(datapts[id], d, 0);
//	printf("For pt %d sphere of radius %f contains %d\n", id, d, sz);
	ANNidxArray nnIdx = new ANNidx[sz];						// allocate near neigh indices
	ANNdistArray dists = new ANNdist[sz];						// allocate near neighbor dists
	kdTree->annkFRSearch(datapts[id], d, sz, nnIdx, dists);

	for (unsigned int j = 0; j < sz; j++) 
	{
//		std::cout << "\t" << j << "\t" << nnIdx[j] << "\t" << dists[j] << "\n";
		clidarr[nnIdx[j]] = clusters.size();
		mem.push_back(nnIdx[j]);
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
