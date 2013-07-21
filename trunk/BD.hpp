#ifndef BD_HPP
#define BD_HPP

#include "Elements.hpp"

class FnCmp
{
	public:
	FnCmp(std::vector<Vertex>*  pv) : pverts(pv) {};
	bool operator ()(const unsigned int & l, const unsigned int & r);
	std::vector<Vertex>* pverts;

};

struct SymBranch
{
	SymBranch() {};
	SymBranch(unsigned int extn, unsigned int sadn, unsigned int cln, unsigned int levn,
			unsigned int bridn, SymBranch* pbr) :
		ext(extn), sad(sadn), lev(levn), bid(bridn), par(pbr), orgbr(0), totch(0), ht(1)
		{
	//		cl.push_back(SplitBr(0,sad,cln));
		};

	unsigned int ext;
	unsigned int sad;
	my_float extw;
	my_float sadw;
//	std::vector<SplitBr> cl;
	unsigned int lev;
	unsigned int bid;
	SymBranch* par;
	SymBranch* orgbr;
	std::list<SymBranch*> ch;
	unsigned int totch;
	unsigned int ht;
	
};

class BD
{
public:
	BD(std::vector<Vertex> & verts);
	SymBranch* BuildSymTree(class ctBranch* b, SymBranch* node, unsigned int & brid, unsigned int & ch, unsigned int & ht);
	void BuildBD();
	void GetNeighbours(unsigned int k, std::vector<unsigned int> & nbrs, unsigned int ftype);
	std::vector<int> & GetVertMap();
	std::vector<Vertex> & m_vlist;
	std::vector<int> vtobrmap;
	std::vector<SymBranch*> bridsarr;
	SymBranch* symroot;
};
#endif

