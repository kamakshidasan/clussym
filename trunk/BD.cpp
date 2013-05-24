#include "BD.hpp"
#include "Utils.hpp"
#include <algorithm>
extern "C" {
#include "tourtre.h"
}


void BD::GetNeighbours(unsigned int k, std::vector<unsigned int> & nbrs, unsigned int ftype)
{
	unsigned int x,y,z;
//	if(ftype == 2)
//		GetCofaces(k+1,nbrs);
//	else
	{
		DeIndex(k,x,y,z);

		assert(x >= 0 && x < SIZEX);
		assert(y >= 0 && y < SIZEY);
		assert(z >= 0 && z < SIZEZ);

		if(y > 0)
		{
			nbrs.push_back(Index(x,y-1,z));
			if(z > 0)
				nbrs.push_back(Index(x,y-1,z-1));
			if(z < SIZEZ-1)
				nbrs.push_back(Index(x,y-1,z+1));
		}
		if(y < SIZEY-1)
		{
			nbrs.push_back(Index(x,y+1,z));
			if(z > 0)
				nbrs.push_back(Index(x,y+1,z-1));
			if(z < SIZEZ-1)
				nbrs.push_back(Index(x,y+1,z+1));
		}
		if(z > 0)
			nbrs.push_back(Index(x,y,z-1));
		if(z < SIZEZ-1)
			nbrs.push_back(Index(x,y,z+1));

		if(x > 0)
		{
			nbrs.push_back(Index(x-1,y,z));
			if(z > 0)
				nbrs.push_back(Index(x-1,y,z-1));
			if(z < SIZEZ-1)
				nbrs.push_back(Index(x-1,y,z+1));
			if(y < SIZEY-1)
			{
				nbrs.push_back(Index(x-1,y+1,z));
				if(z > 0)
					nbrs.push_back(Index(x-1,y+1,z-1));
				if(z < SIZEZ-1)
					nbrs.push_back(Index(x-1,y+1,z+1));
			}
			if(y > 0)
			{
				nbrs.push_back(Index(x-1,y-1,z));
				if(z > 0)
					nbrs.push_back(Index(x-1,y-1,z-1));
				if(z < SIZEZ-1)
					nbrs.push_back(Index(x-1,y-1,z+1));
			}
		}


		if(x < SIZEX -1)
		{
			nbrs.push_back(Index(x+1,y,z));
			if(z > 0)
				nbrs.push_back(Index(x+1,y,z-1));
			if(z < SIZEZ - 1)
				nbrs.push_back(Index(x+1,y,z+1));
			if(y < SIZEY - 1)
			{
				nbrs.push_back(Index(x+1,y+1,z));
				if(z > 0)
					nbrs.push_back(Index(x+1,y+1,z-1));
				if(z < SIZEZ - 1)
					nbrs.push_back(Index(x+1,y+1,z+1));
			}
			if(y > 0)
			{
				nbrs.push_back(Index(x+1,y-1,z));
				if(z > 0)
					nbrs.push_back(Index(x+1,y-1,z-1));
				if(z < SIZEZ - 1)
					nbrs.push_back(Index(x+1,y-1,z+1));
			}

		}
	}

}


size_t neighbours(size_t v, size_t * ar, void* data)
{
	BD* bd = (BD*) data;
	std::vector<unsigned int> nbrs;
	bd->GetNeighbours(v, nbrs, 0);
	memcpy(ar, &nbrs[0], sizeof(unsigned int)*nbrs.size());
	return nbrs.size();
}

my_float value(size_t v, void *data)
{
	BD* bd = (BD*) data;
	return bd->m_vlist[v].w;
}
bool FnCmp::operator()(const unsigned int & l, const unsigned int & r)
{
	return ((*pverts)[l].w < (*pverts)[r].w);
}



BD::BD(std::vector<Vertex> & verts) : m_vlist(verts)
{
}
void BD::BuildBD()
{
	size_t nv = m_vlist.size();
	size_t *ar = new size_t [nv];

	std::vector<unsigned int> vidx(nv);
	for(unsigned int i = 0; i < vidx.size(); i++)
	{
		vidx[i] = i;
	}

	std::sort(vidx.begin(), vidx.end(), FnCmp(&m_vlist));

	for ( int i = 0; i < nv; ++i )
	{
		ar[i] = vidx[i];
	}	
	ctContext* ctx = ct_init(nv, ar, value, neighbours, this);
	ct_sweepAndMerge( ctx );
	ctArc ** arcs = ct_arcMap(ctx );
	ctArc** arcsOut;
	ctNode** nodesOut;
	size_t numarcs, numNodes;
	ct_arcsAndNodes(*arcs, &arcsOut, &numarcs, &nodesOut, &numNodes);
	ctBranch* root = ct_decompose( ctx );

	//create branch decomposition
	ctBranch ** brvertmap;
	brvertmap = ct_branchMap(ctx);

}
