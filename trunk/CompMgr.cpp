#include "CompMgr.hpp"
#include "Cluster.hpp"
#include "BD.hpp"
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <lemon/matching.h>

CompMgr::CompMgr(std::vector<BD*> & pbd) : bd(pbd)
{
}
void CompMgr::Init(std::vector<float> & fnvals)
{
	fvals = fnvals;
	fnmap = std::vector<std::vector<unsigned int> > (fnvals.size());
}
void CompMgr::AddComp(CompNode* c)
{
	unsigned int did1 = c->did;
	assert(c->id == comps.size());
	comps.push_back(c);
	unsigned int sz = fnmap[c->fnid].size();
	for(unsigned int i = 0; i < sz; i++)
	{
		unsigned int compidx = fnmap[c->fnid][i];
		CompNode* other = comps[compidx];
		unsigned int did2 = other->did;
		if(bd[did1]->BrType(c->bid, -1) && bd[did2]->BrType(other->bid, -1))
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
	unsigned int did = c->did;
	unsigned int pfid = ppfid;
	unsigned int bid = c->bid;
	bool done = false;
	SymBranch* b = bd[did]->GetBranch(bid);
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
float CompMgr::Match(CompNode* c1, CompNode* c2, unsigned int & norm, Matrix<float, Dynamic, Dynamic> & A)
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
	
	std::cout<<"Size contribution before matching "<<norm<<std::endl;
	std::cout<<"Children of "<<c1->id<<" ";
	for(; it1 != c1->ch.end(); it1++)
	{
		std::cout<<(*it1)->id<<" "<<"("<<(*it1)->csz<<")";
	}
	std::cout<<std::endl;
	std::cout<<"Children of "<<c2->id<<" ";
	for(; it2 != c2->ch.end(); it2++)
	{
		std::cout<<(*it2)->id<<"("<<(*it2)->csz<<")"<<" ";
	}
	std::cout<<std::endl;

	it1 = c1->ch.begin();
	for(; it1 != c1->ch.end(); it1++)
	{
		LNode n1 = g.addNode();
		unsigned int id1 = (*it1)->id;
		ndmap[n1] = id1;
		it2 = c2->ch.begin();
		unsigned int csz1 = comps[id1]->csz;
		norm += csz1;
		for(; it2 != c2->ch.end(); it2++)
		{
			unsigned int id2 = (*it2)->id;
			unsigned int csz2 = comps[id2]->csz;
			if(create)
			{
				LNode n2 = g.addNode();
				ndmap[n2] = id2;
				nodes[id2] = n2;
				norm +=csz2;
			}
			LNode n2 = nodes[id2];
			LEdge e = g.addEdge(n1,n2);
			edmap[e] = (csz1+csz2)*A(id1,id2);
			std::cout<<"Edge "<<id1<<" "<<id2<<" (sz1+sz2)*wt = w: "<<csz1<<"+"<<csz2<<" * "<<
				A(id1,id2)<<" = "<<edmap[e]<<std::endl;
		}
		create = false;
	}
	
	lemon::MaxWeightedMatching<LGraph,LEdgeMap> mwm(g, edmap);
      	mwm.run();
	float tot = mwm.matchingWeight();
	std::cout<<"Matching "<<tot<<" Size "<<norm<<std::endl;
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
				unsigned int did = comps[cid]->did;
				if(bd[did]->BrType(bid, -1))
					SetParent(comps[cid], fidx+1);
			}
		}
		if(fidx == 0) continue;
		for(unsigned int i = 0; i < fnmap[fidx].size(); i++)
		{	
			float maxval = 0;
			unsigned int cid1 = fnmap[fidx][i];
			unsigned int bid1 = comps[cid1]->bid;
			unsigned int did1 = comps[cid1]->did;
			unsigned int csz1 = comps[cid1]->csz;
			unsigned int norm = csz1;
			std::cout<<"size  of "<<cid1<<" "<<csz1<<std::endl;

			for(unsigned int j = i+1; j < fnmap[fidx].size(); j++)
			{
				unsigned int cid2 = fnmap[fidx][j];
				unsigned int bid2 = comps[cid2]->bid;
				unsigned int did2 = comps[cid2]->did;
				unsigned int csz2 = comps[cid2]->csz;
				std::cout<<"size  of "<<cid2<<" "<<csz2<<std::endl;
				norm += csz2;
				if(bd[did1]->BrType(bid1, -1) && bd[did2]->BrType(bid2, -1))
				{
					float fval = Match(comps[cid1], comps[cid2], norm, U);
					std::cout<<"Par Edge "<<cid1<<" "<<cid2<<" (sz1+sz2)*wt = w: "<<csz1<<"+"<<csz2<<" * "<<
						U(cid1,cid2)<<" = "<<(csz1+csz2)*U(cid1, cid2)<<std::endl;
					U(cid1, cid2) = ((csz1+csz2)*U(cid1, cid2) + fval)/norm;
					U(cid2, cid1) = ((csz1+csz2)*U(cid2, cid1) + fval)/norm;
					if(fval > maxval) 
						maxval = fval;
					std::cout<<"U of "<<cid1<<" "<<cid2<<": "<<U(cid1,cid2)<<" "<<U(cid2, cid1)<<std::endl;
				}
				norm = csz1;
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
void CompMgr::FormLrw(Matrix<float, Dynamic, Dynamic> & Lrw, Matrix<float, Dynamic, Dynamic> & U)
{
	Matrix<float, Dynamic, 1> d = U.rowwise().sum();
	unsigned int rows = d.rows();
	std::cout<<"Rows sum :"<<std::endl<<d<<std::endl;
	for(unsigned int i = 0; i < rows; i++)
		//d(i) = 1.0/sqrt(d(i));
		d(i) = 1.0/d(i);
	Matrix<float, Dynamic, Dynamic> D = d.asDiagonal();
	std::cout<<"Inv Rows sum Matrix:"<<std::endl<<D<<std::endl;
	Matrix<float, Dynamic, Dynamic> I = Matrix<float, Dynamic, Dynamic>::Identity(rows, rows);
//	Lrw = I - D*U*D;
	Lrw = I - D*U;
//	Lrw = D - U;
}
void CompMgr::ClusterComps()
{
	unsigned int csz = comps.size();
	Matrix<float, Dynamic, Dynamic> A = Matrix<float, Dynamic, Dynamic>::Zero(csz, csz);
	Matrix<float, Dynamic, Dynamic> Lrw = Matrix<float, Dynamic, Dynamic>::Zero(csz, csz);
	BuildSimMatrix(A);
	std::cout<<"A matrix:"<<std::endl<<A<<std::endl;
	Matrix<float, Dynamic, Dynamic> U = A;
	//UpSweep(U);
	std::cout<<"U matrix:"<<std::endl<<U<<std::endl;
	FormLrw(Lrw, U);
	std::cout<<"Lrw matrix:"<<std::endl<<Lrw<<std::endl;
	SelfAdjointEigenSolver<Eigen::Matrix<float, Dynamic, Dynamic> > eigs(Lrw);
	std::cout<<"Eigen Values:"<<std::endl<<eigs.eigenvalues()<<std::endl;
	std::cout<<"Eigen Vectors:"<<std::endl<<eigs.eigenvectors()<<std::endl;
	Matrix<float, Dynamic, Dynamic> V = eigs.eigenvectors();
	Matrix<float, Dynamic, Dynamic> D = eigs.eigenvalues().asDiagonal();
	Matrix<float, Dynamic, Dynamic> VD = V;
	std::cout<<"VD matrix:"<<std::endl<<VD<<std::endl;



	unsigned int trunc = csz;
	float diff = 0;
	for(unsigned int i = 1; i < csz; i++)
	{
		if(eigs.eigenvalues()[i] - eigs.eigenvalues()[i-1] > diff)
		{
			diff = eigs.eigenvalues()[i] - eigs.eigenvalues()[i-1];
			trunc = i;

		}
		//std::cout<<"Eigen Vector "<<j<<": "<<std::endl<<eigs.eigenvectors().col(j).transpose()<<std::endl;
		std::cout<<"Eigen Vector "<<i<<": "<<std::endl<<eigs.eigenvalues()[i]*eigs.eigenvectors().col(i).transpose()<<std::endl;
	}
	Matrix<float, Dynamic, Dynamic> eigcords = VD.topLeftCorner(csz, trunc);
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
	/*
	bd[did]->SetBrMask(bid, brmask);
	std::vector<unsigned int> vmask;
	bd[did]->SetVertMask(clid, cid, vmask, brmask, fvals[c->fnid]);
	*/
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
	return exp(-(0.1*diff/std::min(anorm,bnorm)));
	
	
}

