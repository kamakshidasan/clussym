#ifndef COMPMGR_HPP
#define COMPMGR_HPP
#include <vector>
#include <map>
#include <list>
#include <boost/unordered_map.hpp>
#include <eigen3/Eigen/Dense>

struct CompNode
{
	CompNode(unsigned int cid, unsigned int brid, unsigned int fid) : id(cid), bid(brid), fnid(fid) {};
	unsigned int id;
	CompNode* par;
	unsigned int bid;
	std::vector<float> cords;
	unsigned int fnid;
	std::list<CompNode*> ch;
	boost::unordered_map<unsigned int, float> votes;
	float Vote(CompNode* other);
};
using namespace Eigen;
class CompMgr
{
	public:
		CompMgr(unsigned int fsz);
		void AddComp(CompNode* c);
		void Cluster();
	private:
		void BuildSimMatrix(Matrix<float, Dynamic, Dynamic> & A);
		std::vector<CompNode*> comps;
		std::vector<std::vector<unsigned int> > fnmap;
};


#endif
