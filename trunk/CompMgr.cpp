#include "CompMgr.hpp"
#include "BD.hpp"
#include <assert.h>
#include <math.h>
#include <lemon/matching.h>

CompMgr::CompMgr(unsigned int fsz, BD* pbd) : fnmap(fsz), bd(pbd)
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
void CompMgr::SetParent(CompNode* c, unsigned int pfid)
{
	unsigned int bid = c->bid;
	bool done = false;
	SymBranch* b = bd->GetBranch(bid);
	while(!done)
	{
		boost::unordered_map<unsigned int, unsigned int>::iterator cit;
		cit = b->comps.find(pfid);
		if(cit != b->comps.end())
		{
			done = true;
			comps[cit->second]->ch.push_back(c);
			c->par = comps[cit->second];
		}
		else
		{
			b = b->par;
			assert(b->par);
		}
	}
}
float CompMgr::Match(CompNode* c1, CompNode* c2, Matrix<float, Dynamic, Dynamic> & A)
{
	typedef lemon::ListGraph LGraph;
	typedef LGraph::Edge LEdge;
	typedef LGraph::Node LNode;
	typedef float EdgeAttr;
	typedef unsigned int NodeAttr;
	typedef LGraph::EdgeMap<EdgeAttr> LEdgeMap;
	typedef LGraph::NodeMap<NodeAttr> LNodeMap;
	LGraph g;
	LNodeMap ndmap(g);
	LEdgeMap edmap(g);

	std::list<CompNode*>::iterator it1 = c1->ch.begin();
	std::list<CompNode*>::iterator it2 = c2->ch.begin();

	boost::unordered_map<unsigned int, LNode> nodes;
	bool create = true;

	for(; it1 != c1->ch.end(); it1++)
	{
		LNode n1 = g.addNode();
		ndmap[n1] = (*it1)->id;
		for(; it2 != c2->ch.end(); it2++)
		{
			unsigned int id2 = (*it2)->id;
			if(create)
			{
				LNode n2 = g.addNode();
				ndmap[n2] = id2;
				nodes[id2] = n2;
			}
			LNode n2 = nodes[id2];
			LEdge e = g.addEdge(n1,n2);
			edmap[e] = A(ndmap[n1],ndmap[n2]);
		}
	}
	
	lemon::MaxWeightedMatching<LGraph,LEdgeMap> mwm(g, edmap);
      	mwm.run();
	float tot = mwm.matchingWeight();
	return tot;
}

void CompMgr::UpSweep(Matrix<float, Dynamic, Dynamic> & U)
{
	for(unsigned int fidx = 0; fidx < fnmap.size() - 1; fidx++)
	{
		unsigned int fidxsz = fnmap[fidx].size();
		for(unsigned int i = 0; i < fnmap[fidx].size(); i++)
		{
			unsigned int cid = fnmap[fidx][i];
			SetParent(comps[cid], fidx+1);

		}

		for(unsigned int i = 0; i < fnmap[fidx].size() - 1; i++)
		{	
			unsigned int cid1 = fnmap[fidx][i];
			for(unsigned int j = i+1; j < fnmap[fidx].size(); j++)
			{
				unsigned int cid2 = fnmap[fidx][j];
				float fval = Match(comps[cid1], comps[cid2], U);
				U(cid1, cid2) += fval;	

			}
			//ClearChildren(fnmap[fidx][i]);
		}
		//ClearChildren(fnmap[fidx][fidxsz]);
	}
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
	std::cout<<A<<std::endl;
	SelfAdjointEigenSolver<Eigen::Matrix<float, Dynamic, Dynamic> > eigs(A);
	for(unsigned int i = 0; i < csz; i++)
	{

		std::cout<<"Eigen Vector "<<i<<std::endl<<eigs.eigenvalues()[i]*eigs.eigenvectors().col(i).transpose()<<std::endl;
	}
	Matrix<float, Dynamic, Dynamic> U = A;
	UpSweep(U);
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

