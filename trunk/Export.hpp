#ifndef EXPORT_HPP
#define EXPORT_HPP

#include <vtkStructuredPoints.h>
#include <vtkSmartPointer.h>
#include <vtkStructuredPointsWriter.h>
#include <vector>
#include "Elements.hpp"

class Export
{
	public:
		Export(vtkSmartPointer<vtkStructuredPoints> vtkstrpts);
		void CreateVolume();
		void AddToVolume(unsigned int);
		void FlushVolume(const std::string & fn, unsigned int maskf);
	private:
		vtkSmartPointer<vtkStructuredPoints> newstrpts;
		vtkSmartPointer<vtkStructuredPoints> orgstrpts;
		vtkSmartPointer<vtkStructuredPointsWriter> strpwriter;
		std::vector<unsigned int> mask;
		class BD* bd;
};
#endif
