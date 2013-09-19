#ifndef CONTOURS_HPP
#define CONTOURS_HPP

#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include "Elements.hpp"
#include "CompMgr.hpp"
#include <vector>
#include <vtkAppendPolyData.h>
#include <vtkExtractSelection.h>
#include <vtkContourFilter.h>
#include <boost/unordered_map.hpp>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataWriter.h>

class Contours
{
	public:
		Contours(const char* fname, vtkSmartPointer<vtkAppendPolyData> allcontours);
		~Contours();
		void ExtractSymmetry();
	private:
		void ComputeBD(vtkSmartPointer<vtkStructuredPoints> vtkstrpts);
		class BD* bd;
		class CompMgr* compmgr;
		void GenCompCords(CompNode* c, vtkSmartPointer<vtkPolyData> contour);
		void ProcessIsoSurface(unsigned int fid, unsigned int prev, vtkSmartPointer<vtkContourFilter> ctr);
		void GenerateIsoSpace();
		int FindBranchId(vtkSmartPointer<vtkPolyData> contour);
		void SetChildComps(CompNode* c, float curf, float prevf);
		vtkStructuredPointsReader* reader;
		vtkSmartPointer<vtkAppendPolyData> allcts;
		std::vector<Vertex> verts;
		boost::unordered_map<unsigned int, CompNode*> topcomps;
		unsigned int cid;
		std::vector<float> fvals;
		vtkSmartPointer<vtkTriangleFilter> trifil;
		vtkSmartPointer<vtkPolyDataWriter> writer;
};
#endif
