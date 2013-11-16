#include "CompMgr.hpp"
#include "Cluster.hpp"
#include "BD.hpp"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <lemon/matching.h>

CompMgr::CompMgr(std::vector<float> & fnvals, BD* pbd) : fnmap(fnvals.size()), fvals(fnvals), bd(pbd)
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
		if(bd->BrType(c->bid, -1) && bd->BrType(other->bid, -1))
		{
			float orgval = c->Vote(other);
			float val = 0.0;
			//if(orgval > 0.98) val = 1.0;
			//else val = 0.0;
			c->votes[other->id] = orgval;
			printf("Vote(%d %d) = %f %f\n", c->id, other->id, val, orgval);
		}
	}
	fnmap[c->fnid].push_back(c->id);
}
void CompMgr::SetParent(CompNode* c, unsigned int ppfid)
{
	unsigned int pfid = ppfid;
	unsigned int bid = c->bid;
	bool done = false;
	SymBranch* b = bd->GetBranch(bid);
	while(!done)
	{
		if(pfid <= b->topcomp)
		{
			boost::unordered_map<unsigned int, unsigned int>::iterator cit;
			cit = b->comps.find(pfid);
			if(cit != b->comps.end() && cit->second != -1)
			{
				done = true;
				comps[cit->second]->ch.push_back(c);
				c->par = comps[cit->second];
			}
			else
			{
				pfid++;
			}
		}
		else if(b->par)
		{
			b = b->par;
		}
		else
		{
			done = true;
			std::cout<<"Parent not found cid,bid,fid "<<c->id<<" "<<bid<<" "<<ppfid<<std::endl;
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

	std::cout<<"Children of "<<c1->id<<" ";
	for(; it1 != c1->ch.end(); it1++)
	{
		std::cout<<(*it1)->id<<" ";
	}
	std::cout<<std::endl;
	std::cout<<"Children of "<<c2->id<<" ";
	for(; it2 != c2->ch.end(); it2++)
	{
		std::cout<<(*it2)->id<<" ";
	}
	std::cout<<std::endl;

	it1 = c1->ch.begin();
	for(; it1 != c1->ch.end(); it1++)
	{
		LNode n1 = g.addNode();
		unsigned int id1 = (*it1)->id;
		ndmap[n1] = id1;
		it2 = c2->ch.begin();
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
			std::cout<<"Edge "<<id1<<" "<<id2<<" wt:"<< edmap[e]<<std::endl;
		}
		create = false;
	}
	
	lemon::MaxWeightedMatching<LGraph,LEdgeMap> mwm(g, edmap);
      	mwm.run();
	float tot = mwm.matchingWeight();
	std::cout<<"Matching "<<tot<<std::endl;
	return tot;
}

void CompMgr::UpSweep(Matrix<float, Dynamic, Dynamic> & U)
{
	for(unsigned int fidx = 0; fidx < fnmap.size(); fidx++)
	{
		if(fidx < fnmap.size() - 1)
		{
			unsigned int fidxsz = fnmap[fidx].size();
			for(unsigned int i = 0; i < fnmap[fidx].size(); i++)
			{
				unsigned int cid = fnmap[fidx][i];
				unsigned int bid = comps[cid]->bid;
				if(bd->BrType(bid, -1))
					SetParent(comps[cid], fidx+1);
			}
		}
		if(fidx == 0) continue;
		for(unsigned int i = 0; i < fnmap[fidx].size(); i++)
		{	
			float maxval = 0;
			unsigned int cid1 = fnmap[fidx][i];
			unsigned int bid1 = comps[cid1]->bid;
			float norm = 1.0 + comps[cid1]->ch.size();

			for(unsigned int j = i+1; j < fnmap[fidx].size(); j++)
			{
				unsigned int cid2 = fnmap[fidx][j];
				unsigned int bid2 = comps[cid2]->bid;
				if(bd->BrType(bid1, -1) && bd->BrType(bid2, -1))
				{
					float fval = Match(comps[cid1], comps[cid2], U);
					U(cid1, cid2) = (U(cid1, cid2) + fval)/norm;
					U(cid2, cid1) = (U(cid2, cid1) + fval)/norm;
					if(fval > maxval) 
						maxval = fval;
					std::cout<<"U of "<<cid1<<" "<<cid2<<": "<<U(cid1,cid2)<<" "<<U(cid2, cid1)<<std::endl;
				}
			}
			U(cid1,cid1) = 1.0;
		}
	}
}
void CompMgr::BuildSimMatrix(Matrix<float, Dynamic, Dynamic> & A)
{
	std::vector<CompNode*>::iterator it = comps.begin();
	unsigned int csz = comps.size();
	fncords = Matrix<float, Dynamic, 1>::Zero(csz, 1);

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
		fncords(c->id,0) = c->fnid*10;
	}

}

void CompMgr::ClusterComps()
{
	unsigned int csz = comps.size();
	Matrix<float, Dynamic, Dynamic> A = Matrix<float, Dynamic, Dynamic>::Zero(csz, csz);
	BuildSimMatrix(A);
	std::cout<<"A matrix:"<<std::endl<<A<<std::endl;
	Matrix<float, Dynamic, Dynamic> U = A;
	UpSweep(U);
	std::cout<<"U matrix:"<<std::endl<<U<<std::endl;
	SelfAdjointEigenSolver<Eigen::Matrix<float, Dynamic, Dynamic> > eigs(U);
	std::cout<<"Eigen Values:"<<std::endl<<eigs.eigenvalues()<<std::endl;
	std::cout<<"Eigen Vectors:"<<std::endl<<eigs.eigenvectors()<<std::endl;
	Matrix<float, Dynamic, Dynamic> V = eigs.eigenvectors();
	Matrix<float, Dynamic, Dynamic> D = eigs.eigenvalues().asDiagonal();
	Matrix<float, Dynamic, Dynamic> VD = V*D;
	std::cout<<"VD matrix:"<<std::endl<<VD<<std::endl;



	unsigned int trunc = 0;
	for(unsigned int i = 0; i < csz; i++)
	{
		unsigned int j = csz-i-1;
		if(eigs.eigenvalues()[j] < 1)
		{
			trunc = i;
			break;
		}
		//std::cout<<"Eigen Vector "<<j<<": "<<std::endl<<eigs.eigenvectors().col(j).transpose()<<std::endl;
		std::cout<<"Eigen Vector "<<j<<": "<<std::endl<<eigs.eigenvalues()[j]*eigs.eigenvectors().col(j).transpose()<<std::endl;
	}

	Matrix<float, Dynamic, Dynamic> eigcords = VD.topRightCorner(csz, trunc);
	symcords = Matrix<float, Dynamic, Dynamic>(eigcords.rows(), eigcords.cols()+1);
	symcords.topLeftCorner(eigcords.rows(),eigcords.cols()) = eigcords;
	symcords.col(eigcords.cols()) = fncords;
	std::cout<<"Cords: "<<std::endl<<symcords<<std::endl;
	cl = new Cluster(symcords);
	std::vector<unsigned int> & cltrs = cl->GetClusters();
	/*
	for(unsigned int i = 0; i < cltrs.size(); i++)
	{
		if(cltrs[i])
			Export(i, cltrs[i]);
	}
	*/
}

void CompMgr::Export(unsigned int cid, unsigned int clid)
{
	CompNode* c = comps[cid];
	unsigned int bid = c->bid;
	std::vector<unsigned int> brmask;
	bd->SetBrMask(bid, brmask);
	std::vector<unsigned int> vmask;
	bd->SetVertMask(clid, cid, vmask, brmask, fvals[c->fnid]);
}
float CompNode::Vote(CompNode* other)
{
	float diff = 0;
	std::vector<float>::iterator it1 = cords.begin();
	std::vector<float>::iterator it2 = other->cords.begin();
	float anorm = 0, bnorm = 0;
	for(; it1 != cords.end(); it1++, it2++)
	{
		anorm += (*it1)*(*it1);
		bnorm += (*it2)*(*it2);
		diff += (*it1 - *it2)*(*it1 - *it2);
	}
	return exp(-(10*diff/std::min(anorm,bnorm)));
}

