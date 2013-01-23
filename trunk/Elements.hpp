#ifndef ELEMENTS_HPP
#define ELEMENTS_HPP

#include <string.h>
#include <assert.h>
#include <list>
#define my_float double
#ifdef FIVE
const unsigned int dim = 5;
#else
const unsigned int dim = 4;
#endif
struct EdgeCount
{
	EdgeCount(const unsigned int & cnt = 0, my_float s = 0, my_float d = 0.0) : count(cnt), sum(s), dist(d){};
	unsigned int count;
	my_float sum;
	my_float dist;
	my_float pers;
};

struct Vertex
{
	my_float xyz[3];
	my_float w;
	unsigned int cof;
//	unsigned int marked;
	std::list<unsigned int> reg;
//	unsigned int tm;
//	bool bdry;
	Vertex(my_float v[3], const my_float & field = 0)
	{
//		bdry = false;
//		marked = 0; //tm = 0; 
		w = field;
		cof = 0;
		memcpy(xyz, v, sizeof(xyz));
	}
	const Vertex & operator = (const Vertex & vtx)
	{
		this->xyz[0] = vtx.xyz[0];
		this->xyz[1] = vtx.xyz[1];
		this->xyz[2] = vtx.xyz[2];
		this->w = vtx.w;
		return *this;
	}
//	unsigned int GetTm()
//	{ return tm; }
//	void SetTm(unsigned int newtm)
//	{ assert(newtm > tm); tm = newtm;}
//	void SetEmbed() { embed = true; };
	void SetEmbed() { ; };
//	bool Embed() { return embed; };
	bool Embed() { return true; };
	void SetCoface(unsigned int e) { cof = e; }
	unsigned int GetCoface(){return cof;}
	
//	Matrix<my_float> quadric;
};
struct DirEdge
{
	unsigned int node[2];
	DirEdge(unsigned int n[2]) 
	{ memcpy(node, n, sizeof(node));}
	bool operator < (const DirEdge & rhs) const
	{
		return node[0] < rhs.node[0] || 
			((node[0] == rhs.node[0]) && (node[1] < rhs.node[1]));
	}
	bool operator == (const DirEdge & rhs) const
	{
		return node[0] == rhs.node[0] && node[1] == rhs.node[1];
	}
};

struct Edge
{
	unsigned int node[2];
	Edge(unsigned int n[2]) 
	{ memcpy(node, n, sizeof(node));}
	bool operator < (const Edge & rhs) const
	{
		return node[0] < rhs.node[0] || 
			((node[0] == rhs.node[0]) && (node[1] < rhs.node[1]));
	}
	bool operator == (const Edge & rhs) const
	{
		return node[0] == rhs.node[0] && node[1] == rhs.node[1];
	}
};

struct Tri
{
	unsigned int node[3];
	Tri(unsigned int n[3]) 
	{ memcpy(node, n, sizeof(node));}
	bool operator < (const Tri & rhs) const
	{
		return node[0] < rhs.node[0] || 
			((node[0] == rhs.node[0]) && (node[1] < rhs.node[1])) 
			|| (node[0] == rhs.node[0] && node[1] == rhs.node[1]
					&& node[2] < rhs.node[2]);
	}
};

struct Tet
{
	unsigned int node[4];
	Tet() { memset(node, 0, sizeof(node));}
	Tet(unsigned int n[4]) 
	{ memcpy(node, n, sizeof(node));}
	bool operator < (const Tet & r) const
	{
		return (node[0] < r.node[0]) || (node[0] == r.node[0] 
				&& node[1] < r.node[1]) ||
			(node[0] == r.node[0] && node[1] == r.node[1]
			&& node[2] < r.node[2]) ||
		(node[0] == r.node[0] && node[1] == r.node[1] &&
		node[2] == r.node[2] && node[3] < r.node[3]);
	}
};

#endif
