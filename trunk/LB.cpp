#include <assert.h>

#include <vector>

#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

#include <eigen3/Eigen/Dense>

#include <arlsmat.h>
#include <arlssym.h>

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
std::size_t hash_value(const Tri & tri) 
{
	std::size_t h  = 0;
	boost::hash_combine(h, tri.node[0]);
	boost::hash_combine(h, tri.node[1]);
	boost::hash_combine(h, tri.node[2]);
	return h;
}

std::size_t hash_value(const Edge & ed) 
{
	std::size_t h  = 0;
	boost::hash_combine(h, ed.node[0]);
	boost::hash_combine(h, ed.node[1]);
	return h;
}


bool LB::DuplicateEdge(vtkIdType cpts0, vtkIdType cpts1, float val, boost::unordered_map<Edge, unsigned int> & edmap)
{
	unsigned int v[2];
	v[0] = cpts0;
	v[1] = cpts1;
	Order(v[0], v[1]);
	Edge ed(v);
	bool dup = true;
	if(edmap.find(ed) == edmap.end() && val != 0.0)
	{
		edmap[ed] = 1;
		dup = false;
	}

	return dup;

}
bool LB::DuplicateCell(vtkIdType* cpts, boost::unordered_map<Tri, unsigned int> & trimap)
{
	unsigned int v[3];
	v[0] = cpts[0];
	v[1] = cpts[1];
	v[2] = cpts[2];
	Order(v[0], v[1]);	
	Order(v[1], v[2]);	
	Order(v[0], v[1]);

	Tri tri(v);
	bool dup = true;
	if(trimap.find(tri) == trimap.end())
	{
		trimap[tri] = 1;
		dup = false;
	}

	return dup;
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
	//std::cout<< A <<std::endl;
	boost::unordered_map<Tri, unsigned int> trimap;
	boost::unordered_map<Edge, unsigned int> edmap;
	unsigned int nnz = 0;
	while(polys->GetNextCell(npts, cpts))
	{
		assert(npts == 3);
		if(!DuplicateCell(cpts, trimap))
		{
			printf("Points %d %d %d\n", cpts[0],cpts[1],cpts[2]);
			double cot[3];
			GetCots(cpts, pts, cot);

			A(cpts[0], cpts[1]) += 0.5*cot[2];
			A(cpts[1], cpts[0]) += 0.5*cot[2];
			if(!DuplicateEdge(cpts[0], cpts[1], A(cpts[0], cpts[1]), edmap))
			{
				printf("New Edge %d %d (%d)\n", cpts[0], cpts[1], nnz);
				nnz++;
			}

			A(cpts[0], cpts[2]) += 0.5*cot[1];
			A(cpts[2], cpts[0]) += 0.5*cot[1];
			if(!DuplicateEdge(cpts[0], cpts[2], A(cpts[0], cpts[2]), edmap))
			{
				printf("New Edge %d %d (%d)\n", cpts[0], cpts[2], nnz);
				nnz++;
			}

			A(cpts[1], cpts[2]) += 0.5*cot[0];
			A(cpts[2], cpts[1]) += 0.5*cot[0];
			if(!DuplicateEdge(cpts[1], cpts[2], A(cpts[1], cpts[2]), edmap))
			{
				printf("New Edge %d %d (%d)\n", cpts[1], cpts[2], nnz);
				nnz++;
			}

			A(cpts[0], cpts[0]) += 0.5*(cot[2] + cot[1]);
			A(cpts[1], cpts[1]) += 0.5*(cot[2] + cot[0]);
			A(cpts[2], cpts[2]) += 0.5*(cot[1] + cot[0]);
			
			if(!DuplicateEdge(cpts[0], cpts[0], A(cpts[0], cpts[0]), edmap))
			{
				nnz++;
			}
			if(!DuplicateEdge(cpts[1], cpts[1], A(cpts[1], cpts[1]), edmap))
			{
				nnz++;
			}
			if(!DuplicateEdge(cpts[2], cpts[2], A(cpts[2], cpts[2]), edmap))
			{
				nnz++;
			}

		}
		else
		{
			printf("Removed Duplicate\n");
		}
	}
	std::cout<< A <<std::endl;
	float *Acsc = new float[nnz];
	int *irow = new int[nnz];
	int *pcol = new int[numpts+1];
	float *eigv = new float[nnz];
	pcol[0] = 0;
	int idx = 0;
	for(int i = 0; i < numpts; i++)
	{
		for(int j = i; j < numpts; j++)
		{
			if(A(i,j) != 0.0)
			{
				irow[idx] = j;
				Acsc[idx++] = A(i,j);
			}
		}
		pcol[i+1] = idx;
	}
	ARluSymMatrix<float> Am(numpts, nnz, Acsc, irow, pcol); 
	SelfAdjointEigenSolver<Eigen::Matrix<float, Dynamic, Dynamic> > eigs(A);
	std::cout<<"Eigen Values:"<<eigs.eigenvalues().transpose()<<std::endl;
	ARluSymStdEig<float> dprob(2, Am, "SM");
	dprob.FindEigenvectors();
	std::cout<<dprob.Eigenvalue(0)<<" "<<dprob.Eigenvalue(1)<<std::endl;
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
