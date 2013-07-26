#ifndef CONTOURS_HPP
#define CONTOURS_HPP

#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include "Elements.hpp"
#include <vector>
#include <vtkAppendPolyData.h>
#include <vtkExtractSelection.h>
#include <vtkContourFilter.h>

class Contours
{
	public:
		Contours(const char* fname, vtkSmartPointer<vtkAppendPolyData> allcontours);
		~Contours();
		void ExtractSymmetry();
	private:
		void ComputeBD(vtkSmartPointer<vtkStructuredPoints> vtkstrpts);
		class BD* bd;
		void ProcessContour(vtkSmartPointer<vtkPolyData> contour);
		void ProcessIsoSurface(float isoval, vtkSmartPointer<vtkContourFilter> ctr);
		void GenerateIsoSpace(std::vector<float> & fvals);
		void FindBrId(vtkSmartPointer<vtkPolyData> contour);
		vtkStructuredPointsReader* reader;
		vtkSmartPointer<vtkAppendPolyData> allcts;
		std::vector<Vertex> verts;
};
#endif
