#ifndef CONTOURS_HPP
#define CONTOURS_HPP

#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include "Elements.hpp"
#include <vector>
#include <vtkAppendPolyData.h>
#include <vtkExtractSelection.h>

class Contours
{
	public:
		Contours(const char* fname, vtkSmartPointer<vtkAppendPolyData> & allcontours);
		~Contours();
		void ExtractSymmetry();
	private:
		void ComputeBD(vtkSmartPointer<vtkStructuredPoints> vtkstrpts);
		void GenerateContour(vtkSmartPointer<vtkExtractSelection>  extr, class SymBranch* curbr, float isoval);
		void GenerateIsoSpace(SymBranch* curbr, std::vector<float> & fvals, vtkSmartPointer<vtkExtractSelection> extr);
		void GenerateGridIds(SymBranch* curbr, std::vector<unsigned int> & cidarray, float isoval);
		class BD* bd;
		vtkStructuredPointsReader* reader;
		vtkSmartPointer<vtkAppendPolyData> & allcts;
		std::vector<Vertex> verts;
};
#endif
