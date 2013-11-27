#ifndef COMPMGR_HPP
#define COMPMGR_HPP
#include <vector>
#include <map>
#include <list>
#include <boost/unordered_map.hpp>
#include <eigen3/Eigen/Dense>
#include <lemon/list_graph.h>
#include <lemon/concepts/graph.h>
#include <lemon/concepts/maps.h>

struct CompNode
{
	CompNode(unsigned int cid, unsigned int brid, unsigned int fid) : id(cid), bid(brid), fnid(fid) {};
	unsigned int id;
	CompNode* par;
	unsigned int bid;
	std::vector<float> cords;
	unsigned int fnid;
	unsigned int did;
	unsigned int csz;
	std::list<CompNode*> ch;
	boost::unordered_map<unsigned int, float> votes;
	float Vote(CompNode* other);
};
using namespace Eigen;
class CompMgr
{
	public:
		CompMgr(std::vector<float> & fnvals, std::vector<class BD*> & pbd);
		void AddComp(CompNode* c);
		void ClusterComps();
		void Export(unsigned int clid, unsigned int cid);
	private:
		float Match(CompNode* c1, CompNode* c2, unsigned int & norm, Matrix<float, Dynamic, Dynamic> & A);
		void UpSweep(Matrix<float, Dynamic, Dynamic> & A);
		void SetParent(CompNode* c, unsigned int pfid);
		void BuildSimMatrix(Matrix<float, Dynamic, Dynamic> & A);
		void FormLrw(Matrix<float, Dynamic, Dynamic> & Lrw, Matrix<float, Dynamic, Dynamic> & U);
		std::vector<CompNode*> comps;
		std::vector<std::vector<unsigned int> > fnmap;
		std::vector<BD*> & bd;
		class Cluster* cl;
		Matrix<float, Dynamic, Dynamic> symcords;
		Matrix<float, Dynamic, 1> fncords;
		std::vector<float> & fvals;
};


#endif
