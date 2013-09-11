#include "CompMgr.hpp"
#include <assert.h>
#include <math.h>

CompMgr::CompMgr(unsigned int fsz) : fnmap(fsz)
{
}
void CompMgr::AddComp(CompNode* c)
{
	assert(c->id == comps.size());
	comps.push_back(c);
	unsigned int sz = fnmap[c->fnid].size();
	for(unsigned int i = 0; i < sz; i++)
	{
		unsigned int compidx = fnmap[c->fnid][i];
		CompNode* other = comps[compidx];
		float val = c->Vote(other);
		c->votes[other->id] = val;

	}
	fnmap[c->fnid].push_back(c->id);
}

void CompMgr::BuildSimMatrix(Matrix<float, Dynamic, Dynamic> & A)
{
	std::vector<CompNode*>::iterator it = comps.begin();
	for(; it != comps.end(); it++)
	{
		CompNode* c = *it;
		boost::unordered_map<unsigned int, float>::iterator itm = c->votes.begin();
		for(; itm != c->votes.end(); ++itm)
		{
			A(c->id, itm->first) = itm->second;
			A(itm->first, c->id) = itm->second;
		}
		A(c->id, c->id) = 1.0;
	}

}

void CompMgr::Cluster()
{
	unsigned int csz = comps.size();
	Matrix<float, Dynamic, Dynamic> A = Matrix<float, Dynamic, Dynamic>::Zero(csz, csz);
	BuildSimMatrix(A);
	SelfAdjointEigenSolver<Eigen::Matrix<float, Dynamic, Dynamic> > eigs(A);
}
float CompNode::Vote(CompNode* other)
{
	float diff = 0;
	std::vector<float>::iterator it1 = cords.begin();
	std::vector<float>::iterator it2 = other->cords.begin();
	for(; it1 != cords.end(); it1++, it2++)
	{
		diff += (*it1 - *it2)*(*it1 - *it2);
	}
	return exp(-diff);
}

