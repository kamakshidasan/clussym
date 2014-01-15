#ifndef CONTOURS_HPP
#define CONTOURS_HPP

#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkStructuredPoints.h>
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
		Contours(const char* fname1, const char* fname2, vtkSmartPointer<vtkAppendPolyData> allcontours);
		~Contours();
		void ExtractSymmetry(unsigned int inv, unsigned int dcnt);
	private:
		//void Preprocess(vtkSmartPointer<vtkUnstructuredGrid> tgrid, unsigned int inv, unsigned int did);
		void Preprocess(vtkSmartPointer<vtkStructuredPoints> & tgrid, unsigned int inv, unsigned int did);
		std::vector<class BD*> bd;
		class CompMgr* compmgr;
		void GenCompCords(CompNode* c, vtkSmartPointer<vtkPolyData> contour);
		void ProcessIsoSurface(unsigned int fid, unsigned int prev, vtkSmartPointer<vtkContourFilter> ctr, unsigned int did);
		void GenerateIsoSpace(unsigned int did);
		int FindBranchId(vtkSmartPointer<vtkPolyData> contour, float isoval, unsigned int did);
		void SetChildComps(CompNode* c, float curf, float prevf);
		//std::vector<vtkSmartPointer<vtkUnstructuredGrid> > tgrid;
		std::vector<vtkSmartPointer<vtkStructuredPoints> > tgrid;
		vtkSmartPointer<vtkAppendPolyData> allcts;
		std::vector<std::vector<Vertex> > verts;
		boost::unordered_map<unsigned int, CompNode*> topcomps;
		unsigned int cid;
		std::vector<float> fvals;
		vtkSmartPointer<vtkTriangleFilter> trifil;
		vtkSmartPointer<vtkPolyDataWriter> writer;
};
#endif
