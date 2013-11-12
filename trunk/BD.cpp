#include "BD.hpp"
#include "Utils.hpp"
#include <algorithm>
extern "C" {
#include "tourtre.h"
}
#include <string>
#include <stdio.h>
#include <vtkExtractEdges.h>
#include <vtkTriangleFilter.h>
#include <vtkIdList.h>

void BD::GetTNeighbours(unsigned int id, std::vector<unsigned int> & nbrs, unsigned int ftype)
{
	vtkSmartPointer<vtkIdList> cellIdList =	vtkSmartPointer<vtkIdList>::New();
	mesh->GetPointCells(id, cellIdList);
	for(vtkIdType i = 0; i < cellIdList->GetNumberOfIds(); i++)
	{
		vtkSmartPointer<vtkIdList> pointIdList = vtkSmartPointer<vtkIdList>::New();
		mesh->GetCellPoints(cellIdList->GetId(i), pointIdList);
		if(pointIdList->GetId(0) != id)
		      nbrs.push_back(pointIdList->GetId(0));
		else
		      nbrs.push_back(pointIdList->GetId(1));
	}
}
void BD::GetNeighbours(unsigned int k, std::vector<unsigned int> & nbrs, unsigned int ftype)
{
	unsigned int x,y,z;
//	if(ftype == 2)
//		GetCofaces(k+1,nbrs);
//	else
	{
		DeIndex(k,x,y,z, SIZEX, SIZEY, SIZEZ);

		assert(x >= 0 && x < SIZEX);
		assert(y >= 0 && y < SIZEY);
		assert(z >= 0 && z < SIZEZ);

		if(y > 0)
		{
			nbrs.push_back(Index(x,y-1,z, SIZEX, SIZEY, SIZEZ));
			if(z > 0)
				nbrs.push_back(Index(x,y-1,z-1, SIZEX, SIZEY, SIZEZ));
			if(z < SIZEZ-1)
				nbrs.push_back(Index(x,y-1,z+1, SIZEX, SIZEY, SIZEZ));
		}
		if(y < SIZEY-1)
		{
			nbrs.push_back(Index(x,y+1,z, SIZEX, SIZEY, SIZEZ));
			if(z > 0)
				nbrs.push_back(Index(x,y+1,z-1, SIZEX, SIZEY, SIZEZ));
			if(z < SIZEZ-1)
				nbrs.push_back(Index(x,y+1,z+1, SIZEX, SIZEY, SIZEZ));
		}
		if(z > 0)
			nbrs.push_back(Index(x,y,z-1, SIZEX, SIZEY, SIZEZ));
		if(z < SIZEZ-1)
			nbrs.push_back(Index(x,y,z+1, SIZEX, SIZEY, SIZEZ));

		if(x > 0)
		{
			nbrs.push_back(Index(x-1,y,z, SIZEX, SIZEY, SIZEZ));
			if(z > 0)
				nbrs.push_back(Index(x-1,y,z-1, SIZEX, SIZEY, SIZEZ));
			if(z < SIZEZ-1)
				nbrs.push_back(Index(x-1,y,z+1, SIZEX, SIZEY, SIZEZ));
			if(y < SIZEY-1)
			{
				nbrs.push_back(Index(x-1,y+1,z, SIZEX, SIZEY, SIZEZ));
				if(z > 0)
					nbrs.push_back(Index(x-1,y+1,z-1, SIZEX, SIZEY, SIZEZ));
				if(z < SIZEZ-1)
					nbrs.push_back(Index(x-1,y+1,z+1, SIZEX, SIZEY, SIZEZ));
			}
			if(y > 0)
			{
				nbrs.push_back(Index(x-1,y-1,z, SIZEX, SIZEY, SIZEZ));
				if(z > 0)
					nbrs.push_back(Index(x-1,y-1,z-1, SIZEX, SIZEY, SIZEZ));
				if(z < SIZEZ-1)
					nbrs.push_back(Index(x-1,y-1,z+1, SIZEX, SIZEY, SIZEZ));
			}
		}


		if(x < SIZEX -1)
		{
			nbrs.push_back(Index(x+1,y,z, SIZEX, SIZEY, SIZEZ));
			if(z > 0)
				nbrs.push_back(Index(x+1,y,z-1, SIZEX, SIZEY, SIZEZ));
			if(z < SIZEZ - 1)
				nbrs.push_back(Index(x+1,y,z+1, SIZEX, SIZEY, SIZEZ));
			if(y < SIZEY - 1)
			{
				nbrs.push_back(Index(x+1,y+1,z, SIZEX, SIZEY, SIZEZ));
				if(z > 0)
					nbrs.push_back(Index(x+1,y+1,z-1, SIZEX, SIZEY, SIZEZ));
				if(z < SIZEZ - 1)
					nbrs.push_back(Index(x+1,y+1,z+1, SIZEX, SIZEY, SIZEZ));
			}
			if(y > 0)
			{
				nbrs.push_back(Index(x+1,y-1,z, SIZEX, SIZEY, SIZEZ));
				if(z > 0)
					nbrs.push_back(Index(x+1,y-1,z-1, SIZEX, SIZEY, SIZEZ));
				if(z < SIZEZ - 1)
					nbrs.push_back(Index(x+1,y-1,z+1, SIZEX, SIZEY, SIZEZ));
			}

		}
	}

}


size_t neighbours(size_t v, size_t * ar, void* data)
{
	BD* bd = (BD*) data;
	std::vector<unsigned int> nbrs;
	bd->GetTNeighbours(v, nbrs, 0);
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

SymBranch* BD::BuildSymTree(ctBranch* b, SymBranch* node, unsigned int & brid, 
		unsigned int & ch, unsigned int & h)
{
	node->ht = 1;
	node->ext = b->extremum;
	node->sad = b->saddle;
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



BD::BD(std::vector<Vertex> & verts, vtkSmartPointer<vtkUnstructuredGrid> tgrid) : m_vlist(verts)
{
	vtkSmartPointer<vtkExtractEdges> extractEdges =
		vtkSmartPointer<vtkExtractEdges>::New();
	extractEdges->SetInput(tgrid);
	extractEdges->Update();

	mesh = extractEdges->GetOutput();
}
BD::BD(std::vector<Vertex> & verts, int dimx, int dimy, int dimz) : m_vlist(verts), 
	SIZEX(dimx), SIZEY(dimy), SIZEZ(dimz)
{
}

void BD::UpdateSymTree(SymBranch* b, std::vector<unsigned int> & sadidx)
{

	std::list<SymBranch*>::iterator bit = b->ch.begin();
	b->csz = b->sz;
	for(; bit != b->ch.end(); bit++)
	{
		UpdateSymTree(*bit, sadidx);
		b->csz += (*bit)->csz;
	}
	if(b->csz > 120) sadidx.push_back(b->sad);
}
void BD::BuildBD(std::vector<unsigned int> & sadidx)
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
	sadidx.push_back(root->extremum);
	numbr = brid;

	vtobrmap = std::vector<int>(m_vlist.size(), 0);
	for(unsigned int i = 0; i < m_vlist.size(); i++)
	{
		SymBranch* br = (SymBranch*)brvertmap[i]->data;
		vtobrmap[i] = br->bid;
		br->sz++;
	}

	UpdateSymTree(symroot, sadidx);

       	ct_cleanup( ctx );
}

std::vector<int> & BD::GetVertMap()
{
	return vtobrmap;	
}
bool BD::BrType(unsigned int bid, int type)
{
	SymBranch* br = bridsarr[bid];
	unsigned int ext = br->ext;
	unsigned int sad = br->sad;

	if(type == -1 && (bid == 1 || m_vlist[ext].w <= m_vlist[sad].w))
		return true;
	else if(type == 1 && (bid == 1 || m_vlist[ext].w > m_vlist[sad].w))
		return true;
	else
		return false;
}
void BD::AppendExportMask(unsigned int, std::vector<unsigned int> & mask)
{
}
void BD::MaskBranches(SymBranch* br, std::vector<unsigned int> & brmask)
{

	brmask[br->bid] = 1;

	std::list<SymBranch*>::iterator bit = br->ch.begin();
	for(; bit != br->ch.end(); bit++)
	{
		MaskBranches(*bit, brmask);
	}

}
void BD::SetBrMask(unsigned int bid, std::vector<unsigned int> & brmask)
{
	brmask = std::vector<unsigned int> (bridsarr.size(), 0);
	SymBranch* br = bridsarr[bid];
	MaskBranches(br, brmask);
}
void BD::SetVertMask(unsigned int clid, unsigned int cid, std::vector<unsigned int> & vmask, std::vector<unsigned int> & brmask, float fval)
{
	const float minval = m_vlist[symroot->ext].w;
	const float maxval = m_vlist[symroot->sad].w;
	char clnm[50];
	sprintf(clnm, "%d-%d.raw",clid, cid);
	FILE* fpminbin = fopen(clnm, "w+b");

	vmask = std::vector<unsigned int> (m_vlist.size(), 0);
	for(unsigned int i = 0; i < m_vlist.size(); i++)
	{
		if(brmask[vtobrmap[i]] && m_vlist[i].w <= fval)
			vmask[i] = 1;

		if(vmask[i])
		{
			unsigned char val = (m_vlist[i].w - minval)/(maxval - minval)*(255);
			fwrite(&val, sizeof(unsigned char), 1, fpminbin);
		}
		else
		{
			unsigned char val = 0;
			val = 255;
			fwrite(&val, sizeof(unsigned char), 1, fpminbin);
		}
	}
	fclose(fpminbin);

}
