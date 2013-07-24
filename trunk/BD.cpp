#include "BD.hpp"
#include "Utils.hpp"
#include <algorithm>
extern "C" {
#include "tourtre.h"
}

int SIZEX = 64;
int SIZEY = 64;
int SIZEZ = 64;

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
	for(unsigned int i = 0; i < nbrs.size(); i++)
		ar[i] = nbrs[i];
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

SymBranch* BD::BuildSymTree(ctBranch* b, SymBranch* node, unsigned int & brid, unsigned int & ch, unsigned int & h)
{
	node->ht = 1;
	node->ext = b->extremum;
	node->sad = b->saddle;
	node->sadw = m_vlist[node->sad].w;
	node->extw = m_vlist[node->ext].w;
//	node->cl.push_back(SplitBr (0,node->sad, 0));
	node->orgbr = node;
	bridsarr[brid] = node;
	b->data = node;
	node->bid = brid;
	unsigned int totch = 0;
	unsigned int ht = 1;
	for ( ctBranch * c = b->children.head; c != NULL; c = c->nextChild )
	{
		unsigned int chch = 0;
		SymBranch* chnode = new SymBranch;
		chnode->par = node;
		chnode->lev = node->lev+1;
		node->ch.push_back(chnode);
		brid++;
		BuildSymTree(c, chnode, brid, chch, ht);
		if(node->ht < chnode->ht+1) node->ht = chnode->ht+1;
		totch += chch;
	}
	ch = totch + 1;
	node->totch = ch;
	return node;
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

	for(unsigned int i = 0; i < nv; ++i )
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

	unsigned int brid = 1;
	bridsarr = std::vector<SymBranch*>(numNodes,0);
	symroot = new SymBranch;
	symroot->par = NULL;
	symroot->lev = 1;
	SymBranch* node = BuildSymTree(root, symroot,brid, symroot->totch, symroot->ht);
	numbr = brid;

	vtobrmap = std::vector<int>(m_vlist.size(), 0);
	for(unsigned int i = 0; i < m_vlist.size(); i++)
	{
		vtobrmap[i] = ((SymBranch*)brvertmap[i]->data)->bid;
	}

       	ct_cleanup( ctx );
}

std::vector<int> & BD::GetVertMap()
{
	return vtobrmap;	
}
