#include <assert.h>

#include <vector>

#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

#include <boost/numeric/ublas/vector.hpp>

#include "LB.hpp"

void LB::FillMatrix(vtkPolyData* mesh)
{
	namespace ublas = boost::numeric::ublas;

	vtkCellArray* polys = mesh->GetPolys();
	unsigned int n = polys->GetNumberOfCells();
	polys->InitTraversal();
	vtkPoints* pts = mesh->GetPoints();
	vtkIdType npts, *cpts;
	while(polys->GetNextCell(npts, cpts))
	{
		assert(npts == 3);
		double xyz[3][3];
		pts->GetPoint(cpts[0],xyz[0]);
		pts->GetPoint(cpts[1],xyz[1]);
		pts->GetPoint(cpts[2],xyz[2]);

		std::vector<double> v01(3);
		v01[0] = xyz[1][0] - xyz[0][0];
		v01[1] = xyz[1][1] - xyz[0][1];
		v01[2] = xyz[1][2] - xyz[0][2];
		
		std::vector<double> v02(3);
		v02[0] = xyz[2][0] - xyz[0][0];
		v02[1] = xyz[2][1] - xyz[0][1];
		v02[2] = xyz[2][2] - xyz[0][2];

		inner_prod(v01,v02);
		
//		std:vector<double> v01(3);
//		v01(0) = xyz[1][0] - xyz[0][0];
//		v01(1) = xyz[1][1] - xyz[0][1];
//		v01(2) = xyz[1][2] - xyz[0][2];
	}
}

double LB::Cotangent()
{
	double v1v2cos 	= dot(v1,v2);
	double v1sqr 	= dot(v1,v1);
	double v2sqr	= dot(v2v2);
	double v1v2sin 	= sqrt(v1sqr*v2sqr - v1v2cos*v1v2cos);
	double cot 	= v1v2cos / v1v2sin;
	return cot;
}
