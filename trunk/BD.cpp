#include "BD.hpp"
#include "Utils.hpp"
#include <algorithm>

int SIZEX, SIZEY, SIZEZ;
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

BDBranch* BD::BuildSymTree(ctBranch* b, BDBranch* node, unsigned int & brid, unsigned int & ch, unsigned int & h)
{
	node->ht = 1;
	node->ext = b->extremum;
	node->sad = b->saddle;
	node->sadw = m_vlist[node->sad].w;
	node->extw = m_vlist[node->ext].w;
//	node->cl.push_back(SplitBr (0,node->sad, 0));
	node->orgbr = node;
	bridsarr[brid] = brid;
	b->data = &bridsarr[brid];
	node->brid = brid;
	unsigned int totch = 0;
	unsigned int ht = 1;
	for ( ctBranch * c = b->children.head; c != NULL; c = c->nextChild )
	{
		unsigned int chch = 0;
		BDBranch* chnode = new BDBranch;
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
	bridsarr = std::vector<unsigned int>(numNodes,0);
	ctBranch* root = ct_decompose( ctx );

	//create branch decomposition
	ctBranch ** brvertmap;
	brvertmap = ct_branchMap(ctx);
	
	unsigned int brid = 1;
	BDBranch* bdroot = new BDBranch;
	bdroot->par = NULL;
	bdroot->lev = 1;
	BDBranch* node = BuildSymTree(root, bdroot, brid, bdroot->totch, bdroot->ht);
	vtobrmap = std::vector<unsigned int>(m_vlist.size(), 0);
	for(unsigned int i = 0; i < m_vlist.size(); i++)
	{
		vtobrmap[i] = *(int*)brvertmap[i]->data;
	}

	assert(node == bdroot);
	ct_cleanup( ctx );



}
