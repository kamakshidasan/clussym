#include <assert.h>

#include <vector>

#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

#include <eigen3/Eigen/Dense>

#include "LB.hpp"
#include "Utils.hpp"

void LB::GetCots(vtkIdType *cpts, vtkPoints *pts, double cot[])
{
	double xyz[3][3];
	pts->GetPoint(cpts[0],xyz[0]);
	pts->GetPoint(cpts[1],xyz[1]);
	pts->GetPoint(cpts[2],xyz[2]);

	double v01[3], v02[3];
	Sub(xyz[1], xyz[0], v01);
	Sub(xyz[2], xyz[0], v02);
	cot[0] = Cotangent(v01, v02);	

	double v10[3], v12[3];
	v10[0] = -v01[0];
	v10[1] = -v01[1];
	v10[2] = -v01[2];
	Sub(xyz[2], xyz[1], v12);
	cot[1] = Cotangent(v10, v12);	

	double v20[3], v21[3];
	v20[0] = -v02[0];
	v20[1] = -v02[1];
	v20[2] = -v02[2];
	v21[0] = -v12[0];
	v21[1] = -v12[1];
	v21[2] = -v12[2];
	cot[2] = Cotangent(v20, v21);	
}

bool LB::DuplicateCell(vtkIdType* cpts)
{
	unsigned int v[3];
	v[0] = cpts[0];
	v[1] = cpts[1];
	v[2] = cpts[2];
	Order(v[0], v[1]);	
	Order(v[1], v[2]);	
	Order(v[0], v[1]);

	Tri tri(v);
	return true;
}
void LB::FillMatrix(vtkPolyData* mesh)
{
	using namespace Eigen;
	vtkCellArray *polys = mesh->GetPolys();
	unsigned int n = polys->GetNumberOfCells();
	polys->InitTraversal();
	vtkPoints *pts = mesh->GetPoints();
	unsigned int numpts = pts->GetNumberOfPoints();
	vtkIdType npts, *cpts;
	Matrix<float, Dynamic, Dynamic> A = Matrix<float, Dynamic, Dynamic>::Zero(numpts, numpts);
	std::cout<< A <<std::endl;
	while(polys->GetNextCell(npts, cpts))
	{
		assert(npts == 3);
		if(!DuplicateCell(cpts))
		{
			printf("Points %d %d %d\n", cpts[0],cpts[1],cpts[2]);
			double cot[3];
			GetCots(cpts, pts, cot);
			A(cpts[0], cpts[1]) += 0.5*cot[2];
			A(cpts[1], cpts[0]) += 0.5*cot[2];
			A(cpts[0], cpts[2]) += 0.5*cot[1];
			A(cpts[2], cpts[0]) += 0.5*cot[1];
			A(cpts[1], cpts[2]) += 0.5*cot[0];
			A(cpts[2], cpts[1]) += 0.5*cot[0];
			std::cout<< A <<std::endl;
		}
	}
	SelfAdjointEigenSolver<Eigen::Matrix<float, Dynamic, Dynamic> > eigs(A);
}

double LB::Cotangent(double v1[], double v2[])
{
	double v1v2cos 	= Dot(v1,v2);
	double v1sqr 	= Dot(v1,v1);
	double v2sqr	= Dot(v2,v2);
	double v1v2sin 	= sqrt(v1sqr*v2sqr - v1v2cos*v1v2cos);
	double cot 	= v1v2cos / v1v2sin;
	return cot;
}
