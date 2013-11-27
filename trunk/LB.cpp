#include <sys/time.h>

#include <assert.h>

#include <vector>

#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

#include <eigen3/Eigen/Dense>

#include <arlsmat.h>
#include <arlssym.h>

//#include <flens/flens.cxx>

#include "LB.hpp"
#include "Utils.hpp"

#include <boost/foreach.hpp>
#undef NDEBUG
void LB::GetCotsLensArea(vtkIdType *cpts, vtkPoints *pts, double cot[], double len[], double & a)
{
	double xyz[3][3];
	double area[3];
	pts->GetPoint(cpts[0],xyz[0]);
	pts->GetPoint(cpts[1],xyz[1]);
	pts->GetPoint(cpts[2],xyz[2]);

//	printf("p1 %lf %lf %lf\n", xyz[0][0], xyz[0][1], xyz[0][2]);
//	printf("p2 %lf %lf %lf\n", xyz[1][0], xyz[1][1], xyz[1][2]);
//	printf("p3 %lf %lf %lf\n", xyz[2][0], xyz[2][1], xyz[2][2]);
	double v01[3], v02[3];
	Sub(xyz[1], xyz[0], v01);
	Sub(xyz[2], xyz[0], v02);

	len[2] = Dot(v01, v01);
	cot[0] = Cotangent(v01, v02, area[0]);	

	double v10[3], v12[3];
	v10[0] = -v01[0];
	v10[1] = -v01[1];
	v10[2] = -v01[2];
	Sub(xyz[2], xyz[1], v12);

	len[0] = Dot(v12, v12);
	cot[1] = Cotangent(v10, v12, area[1]);	

	double v20[3], v21[3];
	v20[0] = -v02[0];
	v20[1] = -v02[1];
	v20[2] = -v02[2];
	v21[0] = -v12[0];
	v21[1] = -v12[1];
	v21[2] = -v12[2];

	len[1] = Dot(v20, v20);
	cot[2] = Cotangent(v20, v21, area[2]);

//	printf("Area %lf %lf %lf\n",area[0],area[1],area[2]);
	assert(fabs(area[0] - area[1]) < 0.000001);
	assert(fabs(area[1] - area[2]) < 0.000001);
	assert(fabs(area[0] - area[2]) < 0.000001);


	a = area[0];
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


bool LB::DuplicateEdge(vtkIdType cpts0, vtkIdType cpts1, boost::unordered_map<Edge, unsigned int> & edmap)
{
	unsigned int v[2];
	v[0] = cpts0;
	v[1] = cpts1;
	Order(v[0], v[1]);
	Edge ed(v);
	bool dup = true;
	if(edmap.find(ed) == edmap.end())
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


void LB::GetVorArea(double cot[3], double len[3], double triarea, double vorarea[3])
{
	if(cot[0] > 0 && cot[1] > 0 && cot[2] > 0)
	{
		vorarea[0] = (cot[1]*len[1]+cot[2]*len[2])/8;
		vorarea[1] = (cot[2]*len[2]+cot[0]*len[0])/8;
		vorarea[2] = (cot[0]*len[0]+cot[1]*len[1])/8;
	}
	else
	{
		if(cot[0] < 0)
		{
			vorarea[0] = triarea/2;
			vorarea[1] = triarea/4;
			vorarea[2] = triarea/4;
		}
		else if(cot[1] < 0)
		{
			vorarea[0] = triarea/4;
			vorarea[1] = triarea/2;
			vorarea[2] = triarea/4;
		}
		else if(cot[2] < 0)
		{
			vorarea[0] = triarea/4;
			vorarea[1] = triarea/4;
			vorarea[2] = triarea/2;

		}
	}

	//printf("VorArea %f %f %f\n", vorarea[0], vorarea[1], vorarea[2]);
}

void LB::GetEigen(vtkPolyData* mesh, std::vector<float> & cords)
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
	std::vector<double> ndvorarea(numpts, 0.0);
	//flens::GeMatrix<flens::FullStorage<float, flens::ColMajor> > //flM(numpts, numpts), VL(numpts, numpts), VR(numpts, numpts);
	unsigned int nnz = 0;
	double sumtriarea = 0.0;
	while(polys->GetNextCell(npts, cpts))
	{
		assert(npts == 3);
		if(!DuplicateCell(cpts, trimap))
		{
//			printf("3 %d %d %d\n", cpts[0],cpts[1],cpts[2]);
			double cot[3], len[3], vorarea[3], triarea;
			GetCotsLensArea(cpts, pts, cot, len, triarea);
			
			//assert(cot[0] != 0.0 && cot[1] != 0.0 && cot[2] != 0.0);
			assert(len[0] != 0.0 && len[1] != 0.0 && len[2] != 0.0);

			GetVorArea(cot, len, triarea, vorarea);

			sumtriarea += triarea;
			//assert(fabs(vorarea[0]+vorarea[1]+vorarea[2]-triarea < 0.001));
			ndvorarea[cpts[0]] += vorarea[0];
			ndvorarea[cpts[1]] += vorarea[1];
			ndvorarea[cpts[2]] += vorarea[2];
	
			//printf("VorArea for (%d %d %d) %f %f %f\n", cpts[0], cpts[1], cpts[2], ndvorarea[cpts[0]], ndvorarea[cpts[1]], ndvorarea[cpts[2]]);

			A(cpts[0], cpts[1]) -= cot[2];//2;
			A(cpts[1], cpts[0]) -= cot[2];//2;
			//flM(cpts[0]+1, cpts[1]+1) += 0.5*cot[2];
			//flM(cpts[1]+1, cpts[0]+1) += 0.5*cot[2];
			if(!DuplicateEdge(cpts[0], cpts[1], edmap))
			{
//				printf("New Edge %d %d (%d)\n", cpts[0], cpts[1], nnz);
				nnz++;
			}

			//flM(cpts[0]+1, cpts[2]+1) += 0.5*cot[1];
			//flM(cpts[2]+1, cpts[0]+1) += 0.5*cot[1];
			A(cpts[0], cpts[2]) -= cot[1];//2;
			A(cpts[2], cpts[0]) -= cot[1];//2;
			if(!DuplicateEdge(cpts[0], cpts[2], edmap))
			{
//				printf("New Edge %d %d (%d)\n", cpts[0], cpts[2], nnz);
				nnz++;
			}

			//flM(cpts[1]+1, cpts[2]+1) += 0.5*cot[0];
			//flM(cpts[2]+1, cpts[1]+1) += 0.5*cot[0];
			A(cpts[1], cpts[2]) -= cot[0];//2;
			A(cpts[2], cpts[1]) -= cot[0];//2;
			if(!DuplicateEdge(cpts[1], cpts[2], edmap))
			{
//				printf("New Edge %d %d (%d)\n", cpts[1], cpts[2], nnz);
				nnz++;
			}

			//flM(cpts[0]+1, cpts[0]+1) += 0.5*(cot[2] + cot[1]);
			//flM(cpts[1]+1, cpts[1]+1) += 0.5*(cot[2] + cot[0]);
			//flM(cpts[2]+1, cpts[2]+1) += 0.5*(cot[1] + cot[0]);
			A(cpts[0], cpts[0]) += (cot[2] + cot[1]);//2;
			A(cpts[1], cpts[1]) += (cot[2] + cot[0]);//2;
			A(cpts[2], cpts[2]) += (cot[1] + cot[0]);//2;
			
			if(!DuplicateEdge(cpts[0], cpts[0], edmap))
			{
				nnz++;
			}
			if(!DuplicateEdge(cpts[1], cpts[1], edmap))
			{
				nnz++;
			}
			if(!DuplicateEdge(cpts[2], cpts[2], edmap))
			{
				nnz++;
			}

		}
		else
		{
//			printf("Removed Duplicate\n");
		}
	}
	typedef std::pair<Edge, unsigned int> edmaptype;
	unsigned int nnzcnt = 0;
	double sumvorarea = 0.0;
	for(unsigned int i = 0; i < numpts; i++)
	{
/*		for(unsigned int j = i; j < numpts; j++)
		{
			if(A(i,j) != 0.0)
			{
				double aij = A(i,j);
				double aji = A(j,i);
				A(i,j) = (aij/(2*ndvorarea[i]) + aji/(2*ndvorarea[j]))/2;
				A(j,i) = (aij/(2*ndvorarea[i]) + aji/(2*ndvorarea[j]))/2;
			}
		}
		sumvorarea += ndvorarea[i];*/
		double xyz[3];
	//	pts->GetPoint(i,xyz);
	//	printf("%f %f %f\n", xyz[0],xyz[1],xyz[2]);
	}
	//assert(fabs(sumvorarea - sumtriarea) < 0.00001);
//	printf("Tot area %f\n", sumtriarea);
	BOOST_FOREACH(edmaptype i, edmap) 
	{
		nnzcnt++;
//		float aij = A(i.first.node[0], i.first.node[1]);
//		float aji = A(i.first.node[1], i.first.node[0]);

		A(i.first.node[0], i.first.node[1]) /= 2*sqrt(ndvorarea[i.first.node[0]]*ndvorarea[i.first.node[1]]);
		if(i.first.node[0] != i.first.node[1])
			A(i.first.node[1], i.first.node[0]) /= 2*sqrt(ndvorarea[i.first.node[1]]*ndvorarea[i.first.node[0]]);
//		aij = A(i.first.node[0], i.first.node[1]);
//		aji = A(i.first.node[1], i.first.node[0]);
		assert(fabs(A(i.first.node[0], i.first.node[1]) - A(i.first.node[1], i.first.node[0])) < 0.000001);
	}
	//assert(nnzcnt == nnz);
//	std::cout<< A <<std::endl;
	struct timeval timeval_start, timeval_end;
	gettimeofday(&timeval_start, NULL);
	SelfAdjointEigenSolver<Eigen::Matrix<float, Dynamic, Dynamic> > eigs(A, EigenvaluesOnly) ;
//	EigenSolver<Eigen::Matrix<float, Dynamic, Dynamic> > eigs(A, EigenvaluesOnly) ;
	gettimeofday(&timeval_end, NULL);
	double time_start = timeval_start.tv_sec + (double) timeval_start.tv_usec/1000000;
	double time_end= timeval_end.tv_sec + (double) timeval_end.tv_usec/1000000;

	printf("Eigen Values :\n");
	for(unsigned int i = 0; i < 10; i++)
	{
		cords.push_back(1.0/(eigs.eigenvalues()[i+1]));
		printf("%f ", cords[i]);
	}
	printf("\n");
//	std::cout<<"Eigen EigenVals: "<<eigs.eigenvalues().transpose()<<std::endl;
//	std::cout<<"Eigen Time: "<<time_end - time_start<<std::endl;
	/*float *Acsc = new float[nnz];
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
	gettimeofday(&timeval_start, NULL);
	ARluSymMatrix<float> Am(numpts, nnz, Acsc, irow, pcol); 
	ARluSymStdEig<float> dprob(10, Am, "SM");
	int nconv = dprob.FindEigenvectors();
	gettimeofday(&timeval_end, NULL);
	time_start = timeval_start.tv_sec + (double) timeval_start.tv_usec/1000000;
	time_end= timeval_end.tv_sec + (double) timeval_end.tv_usec/1000000;
	std::cout<<"Arpack EigenVals:";
	for(unsigned i = 0; i < nconv; i++)
	{
		std::cout<<dprob.Eigenvalue(i)<<" ";
	}
	std::cout<<std::endl;
	std::cout<<"Arpack Time: "<<time_end - time_start<<std::endl;*/
//	flens::DenseVector<flens::Array<float> > wr(numpts), wi(numpts), work;
//	flens::lapack::ev(true, true, //flM, wr, wi, VL, VR, work);

//	std::cout << "wr = " << wr << std::endl;
}

double LB::Cotangent(double v1[], double v2[], double & area)
{
//	printf("v1 %lf %lf %lf\n", v1[0], v1[1], v1[2]);
//	printf("v2 %lf %lf %lf\n", v2[0], v2[1], v2[2]);
	double v1v2cos 	= Dot(v1,v2);
	double v1sqr 	= Dot(v1,v1);
	double v2sqr	= Dot(v2,v2);
	double v1v2sin 	= sqrt(v1sqr*v2sqr - v1v2cos*v1v2cos);
	double cot 	= v1v2cos / v1v2sin;
	area = v1v2sin/2;
//	printf("v1v2cos %lf v1sqr %lf v2sqr %lf v1v2sin %lf area %lf  cot %lf\n", v1v2cos, v1sqr, v2sqr, v1v2sin, area, cot);
	return cot;
}
