#ifndef BD_HPP
#define BD_HPP

#include "Elements.hpp"
extern "C" {
#include "tourtre.h"
}

class FnCmp
{
	public:
	FnCmp(std::vector<Vertex>*  pv) : pverts(pv) {};
	bool operator ()(const unsigned int & l, const unsigned int & r);
	std::vector<Vertex>* pverts;

};

struct BDBranch
{
	BDBranch() {};
	BDBranch(unsigned int extn, unsigned int sadn, unsigned int cln, unsigned int levn,
			unsigned int bridn, BDBranch* pbr) :
		ext(extn), sad(sadn), lev(levn), brid(bridn), par(pbr), orgbr(0), totch(0), ht(1)
		{
			//cl.push_back(SplitBr(0,sad,cln));
		};

	unsigned int ext;
	unsigned int sad;
	my_float extw;
	my_float sadw;
//	std::vector<SplitBr> cl;
	unsigned int lev;
	unsigned int brid;
	BDBranch* par;
	BDBranch* orgbr;
	std::list<BDBranch*> ch;
	unsigned int totch;
	unsigned int ht;
	
};

class BD
{
public:
	BD(std::vector<Vertex> & verts);
	void BuildBD();
	void GetNeighbours(unsigned int k, std::vector<unsigned int> & nbrs, unsigned int ftype);
	BDBranch* BuildSymTree(ctBranch* b, BDBranch* node, unsigned int & brid, unsigned int & ch, unsigned int & h);
	std::vector<Vertex> & m_vlist;
	std::vector<unsigned int> bridsarr;
	std::vector<unsigned int> vtobrmap;
};
#endif

