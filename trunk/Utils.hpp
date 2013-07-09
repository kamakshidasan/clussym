#ifndef UTILS_HPP
#define UTILS_HPP

#include "Elements.hpp"
#include <set>
#include <math.h>

extern int SIZEX;
extern int SIZEY;
extern int SIZEZ;

template <class T> inline long Round(T x) 
{
	return ((x) >=0 ? (long)((x)+0.5) : (long)((x)-0.5));
}
template <class T> inline T Min(T a, T b)
{
	return a < b ? a : b;
}
template <class T> inline T Max(T a, T b)
{
	return a > b ? a : b;
}
inline void Sub(my_float u[3], my_float v[3], my_float res[3])
{
	res[0] = u[0] - v[0];
	res[1] = u[1] - v[1];
	res[2] = u[2] - v[2];
}
inline my_float NormDot(my_float u[3], my_float v[3])
{
	return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2])/(sqrt(u[0]*u[0]+u[1]*u[1]+u[2]*u[2])*sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]));
}

inline my_float Dot(my_float u[3], my_float v[3])
{
	return (u[0]*v[0] + u[1]*v[1] + u[2]*v[2]);
}

inline void Cross(my_float u[3], my_float v[3], my_float n[3])
{
	n[0] = u[1]*v[2] - u[2]*v[1];
	n[1] = u[2]*v[0] - u[0]*v[2];
	n[2] = u[0]*v[1] - u[1]*v[0];
}

inline void MatAdd(my_float* A, my_float* B, my_float* C, unsigned int n, int sign = 1)
{
	for(unsigned int i = 0; i < n; i++)
	{
		if(sign > 0)
			C[i] = A[i] + B[i];
		else if(sign <= 0)
			C[i] = A[i] - B[i];
	}
}
inline void MultMv(my_float *NewV, const my_float* M, const my_float *V, unsigned int n)
{
	for(unsigned int i = 0; i < n; i++)
	{
		NewV[i] = 0;
		for(unsigned int j = 0; j < n; j++)
		{
			NewV[i] += M[i+j*n]*V[j];
		}
	}
}

/*
 * If (a > b) Swap(a,b)
 */
inline void Order(unsigned int & a, unsigned int & b)
{
	if(a > b)
	{
		unsigned int t = a;
		a = b;
		b = t;
	}
}

inline void SplitEdges(const Edge & ed, std::set<unsigned > & nodes)
{
	
	nodes.insert(ed.node[0]);
	nodes.insert(ed.node[1]);
}
inline void SplitTri(const Tri & tri, std::set<Edge> & edges)
{
	unsigned int node[2];


	node[0] = tri.node[0];
	node[1] = tri.node[1];
	Order(node[0], node[1]);
	Edge ed(node);
	edges.insert(ed);


	node[0] = tri.node[1];
	node[1] = tri.node[2];
	Order(node[0], node[1]);
	ed = Edge(node);
	edges.insert(ed);

	node[0] = tri.node[2];
	node[1] = tri.node[0];
	Order(node[0], node[1]);
	ed = Edge(node);
	edges.insert(ed);
}
inline void SplitTri(unsigned int tri[], unsigned int edges[][2])
{


	edges[0][0] = tri[0];
	edges[0][1] = tri[1];
	Order(edges[0][0], edges[0][1]);
	edges[1][0] = tri[0];
	edges[1][1] = tri[2];
	Order(edges[1][0], edges[1][1]);
	edges[2][0] = tri[1];
	edges[2][1] = tri[2];
	Order(edges[2][0], edges[2][1]);

}

unsigned int Index(unsigned int x, unsigned int y, unsigned int z);

void DeIndex(unsigned int i, unsigned int & x, unsigned int & y, unsigned int & z);
#endif
