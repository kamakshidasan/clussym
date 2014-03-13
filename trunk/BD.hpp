#ifndef BD_HPP
#define BD_HPP

#include "Elements.hpp"
#include <boost/unordered_map.hpp>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
extern "C" {
#include "tourtre.h"
}

extern unsigned int fsz;
class FnCmp
{
	public:
	FnCmp(std::vector<Vertex>*  pv) : pverts(pv) {};
	bool operator ()(const unsigned int & l, const unsigned int & r);
	std::vector<Vertex>* pverts;

};
class ArcCmp
{
	public:
	ArcCmp(std::vector<Vertex>*  pv) : pverts(pv) {};
	bool operator ()(ctArc* const &  l, ctArc* const & r);
	std::vector<Vertex>* pverts;

};
struct NodeData
{
	NodeData(unsigned int s = 0, unsigned int n = 0) : sz(s), num(n) {};
	unsigned int sz;
	unsigned int num;
};
struct SymBranch
{
	SymBranch() : sz(0), csz(0), orgbr(0), totch(0), ht(1) {};
	SymBranch(unsigned int extn, unsigned int sadn, unsigned int fsz, unsigned int levn,
			unsigned int bridn, SymBranch* pbr) :
		ext(extn), sad(sadn), lev(levn), bid(bridn), sz(0), csz(0), par(pbr), orgbr(0), totch(0), ht(1)
		{
	//		cl.push_back(SplitBr(0,sad,cln));
		};

	unsigned int ext;
	unsigned int sad;
	unsigned int lev;
	unsigned int bid;
	unsigned int sz;
	unsigned int csz;
	SymBranch* par;
	SymBranch* orgbr;
	std::list<SymBranch*> ch;
	unsigned int totch;
	unsigned int ht;
	boost::unordered_map<unsigned int, unsigned int> comps;
	unsigned int topcomp;
	
};

class BD
{
public:
	BD(std::vector<Vertex> & verts, int dimx, int dimy, int dimz);
	BD(std::vector<Vertex> & verts, vtkSmartPointer<vtkUnstructuredGrid> tgrid);
	void BuildBD(std::vector<unsigned int> & sadidx, std::vector<float> & isovals, float alpha);
	unsigned int NumBr() { return numbr; };
	SymBranch* BuildSymTree(class ctBranch* b, SymBranch* node, unsigned int & brid, 
			unsigned int & ch, unsigned int & ht);
	void GetNeighbours(unsigned int k, std::vector<unsigned int> & nbrs, unsigned int ftype);
	void GetTNeighbours(unsigned int k, std::vector<unsigned int> & nbrs, unsigned int ftype);
	SymBranch* GetBranch(unsigned int bid) { return bridsarr[bid]; };
	void AppendExportMask(unsigned int, std::vector<unsigned int> & mask);
	bool BrType(unsigned int bid, int type);	
	void SetBrMask(unsigned int bid, std::vector<unsigned int> & brmask, float fval);
	void SetVertMask(unsigned int cid, unsigned int bid, std::vector<unsigned int> & vmask, std::vector<unsigned int> & brmask, float fval);
	void MaskBranches(SymBranch* br, std::vector<unsigned int> & brmask, float fval);
	void UpdateSymTree(SymBranch* b, std::vector<unsigned int> & sadidx);
	void Sample(std::vector<float> & isovals, float alpha);
	std::vector<int> & GetVertMap();
	std::vector<Vertex> & m_vlist;
	std::vector<int> vtobrmap;
	std::vector<int> arcdata;
	std::vector<NodeData> nodedata;
	std::vector<SymBranch*> bridsarr;
	SymBranch* symroot;
	unsigned int numbr;
	int SIZEX;
	int SIZEY;
	int SIZEZ;
	unsigned int first;
	vtkSmartPointer<vtkPolyData> mesh;
	ctArc** arcsOut;
	size_t numarcs;
	ctNode** nodesOut;
	size_t numNodes;
};
#endif

