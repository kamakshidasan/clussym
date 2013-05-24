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


class BD
{
public:
	BD(std::vector<Vertex> & verts);
	void BuildBD();
	void GetNeighbours(unsigned int k, std::vector<unsigned int> & nbrs, unsigned int ftype);
	std::vector<Vertex> & m_vlist;
};
#endif

