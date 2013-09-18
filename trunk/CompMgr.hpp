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
	unsigned int csz;
	std::list<CompNode*> ch;
	boost::unordered_map<unsigned int, float> votes;
	float Vote(CompNode* other);
};
using namespace Eigen;
class CompMgr
{
	public:
		CompMgr(unsigned int fsz, class BD* pbd);
		void AddComp(CompNode* c);
		void ClusterComps();
		void Group(unsigned int);
	private:
		float Match(CompNode* c1, CompNode* c2, Matrix<float, Dynamic, Dynamic> & A);
		void UpSweep(Matrix<float, Dynamic, Dynamic> & A);
		void SetParent(CompNode* c, unsigned int pfid);
		void BuildSimMatrix(Matrix<float, Dynamic, Dynamic> & A);
		std::vector<CompNode*> comps;
		std::vector<std::vector<unsigned int> > fnmap;
		BD* bd;
		class Cluster* cl;
		Matrix<float, Dynamic, Dynamic> symcords;


};


#endif
