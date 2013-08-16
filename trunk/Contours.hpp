#ifndef CONTOURS_HPP
#define CONTOURS_HPP

#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include "Elements.hpp"
#include <vector>
#include <vtkAppendPolyData.h>
#include <vtkExtractSelection.h>
#include <vtkContourFilter.h>
#include <boost/unordered_map.hpp>

struct CompNode
{
	CompNode(unsigned int cid, unsigned int brid, float fval) : id(cid), bid(brid), fn(fval) {};
	unsigned int id;
	CompNode* par;
	unsigned int bid;
	std::vector<float> cords;
	float fn;
	std::list<CompNode*> ch;	
};

class Contours
{
	public:
		Contours(const char* fname, vtkSmartPointer<vtkAppendPolyData> allcontours);
		~Contours();
		void ExtractSymmetry();
	private:
		void ComputeBD(vtkSmartPointer<vtkStructuredPoints> vtkstrpts);
		class BD* bd;
		void GenCompCords(CompNode* c, vtkSmartPointer<vtkPolyData> contour);
		void ProcessIsoSurface(float isoval, float prevf, vtkSmartPointer<vtkContourFilter> ctr);
		void GenerateIsoSpace(std::vector<float> & fvals);
		int FindBranchId(vtkSmartPointer<vtkPolyData> contour);
		void SetChildComps(CompNode* c, float curf, float prevf);
		vtkStructuredPointsReader* reader;
		vtkSmartPointer<vtkAppendPolyData> allcts;
		std::vector<Vertex> verts;
		boost::unordered_map<unsigned int, CompNode*> topcomps;
		unsigned int cid;
};
#endif
