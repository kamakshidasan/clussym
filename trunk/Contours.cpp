#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangleFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkPLYWriter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkDecimatePro.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkThreshold.h>
#include <vtkVolume.h>
#include <vtkLinearSubdivisionFilter.h>
#include <vtkSelection.h>
#include <vtkSelectionNode.h>
#include <vtkExtractSelection.h>
#include <vtkInformation.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>

#include <deque>

#include "BD.hpp"
#include "LB.hpp"
#include "Cluster.hpp"
#include "Remesh.hpp"
#include "Contours.hpp"


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;


void vtktoPointList(PointList & pl, vtkPolyData* mesh)
{
	vtkSmartPointer<vtkFloatArray> normals = 
		vtkFloatArray::SafeDownCast(mesh->GetPointData()->GetArray("Normals"));
	
	for(vtkIdType i = 0; i < mesh->GetNumberOfPoints(); i++) 
	{ 
		double p[3]; 
		mesh->GetPoint(i, p);
		Point pt(p[0], p[1], p[2]);
//		printf("%lf %lf %lf\n", p[0],p[1],p[2]);
		if(normals)
		{
			normals->GetTuple(i, p);
		//	printf("%lf %lf %lf\n", p[0],p[1],p[2]);
			Vector n(p[0], p[1], p[2]);
			pl.push_back(Point_with_normal(pt,n));
		}
		else
		{
			printf("No normals\n");
			exit(0);
		}
	} 

}
void Contours::FindBrId(vtkSmartPointer<vtkPolyData> contour)
{
	vtkSmartPointer<vtkCellArray> cells = contour->GetPolys();
	cells->InitTraversal();
	vtkIdType* cellpts;
	vtkIdType  numpts;

	while(cells->GetNextCell(numpts, cellpts))
	{
		double centre[3], a[3], b[3], c[3];
		
		contour->GetPoint(cellpts[0], a);
		contour->GetPoint(cellpts[1], b);
		contour->GetPoint(cellpts[2], c);

		centre[0] = (a[0]+b[0]+c[0])/3.0;
		centre[1] = (a[1]+b[1]+c[1])/3.0;
		centre[2] = (a[2]+b[2]+c[2])/3.0;

		 double p[3] = {1.0,1.0,1.0}; 
		 int subid; 
		 double pcoords[3] = {0,0,0}; 
		 double weights[8]; 
		 contour->FindCell(centre,0,0,0.001,subid,pcoords,weights);
	}
}
void Contours::ProcessContour(vtkSmartPointer<vtkPolyData> contour)
{
	vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
	vtkSmartPointer<vtkPolyDataNormals> gennorms = vtkSmartPointer<vtkPolyDataNormals>::New();
	vtkSmartPointer<vtkCleanPolyData> cleanpolydata = vtkSmartPointer<vtkCleanPolyData>::New();
	if(contour->GetNumberOfPoints() > 20)
	{
		PointList pl;
		gennorms->SetInput(contour);
		gennorms->ComputePointNormalsOn();
		gennorms->ComputeCellNormalsOff();
		gennorms->SplittingOff();
		gennorms->NonManifoldTraversalOff();
		gennorms->Update();

		contour = gennorms->GetOutput();
		vtktoPointList(pl, contour);
		vtkPolyData* mesh = Remesh(pl);

		cleanpolydata->SetInput(mesh);
		cleanpolydata->Update();
		contour = cleanpolydata->GetOutput();
		FindBrId(contour);
		unsigned int ntri = contour->GetNumberOfPolys();
		if(ntri > 20000)
		{
			decimate->SetInputConnection(contour->GetProducerPort());
			float target = 1 - 20000.0/ntri;
			if(target > 0.8) target = 0.8;
			decimate->SetTargetReduction(target);
			decimate->Update();
			contour = decimate->GetOutput();				
			std::cout <<"After decimate "<<contour->GetNumberOfCells()
				  <<" "<<contour->GetNumberOfPolys()
				  <<" "<<contour->GetNumberOfPoints()<<std::endl;
		}

		std::vector<std::vector<double> > surfcords;
		LB lb;
		lb.GetEigen(contour, surfcords);
		allcts->AddInputConnection(contour->GetProducerPort());

		mesh->Delete();
	}
}
void Contours::ProcessIsoSurface(float isoval, vtkSmartPointer<vtkContourFilter> ctr)
{
	vtkSmartPointer<vtkPolyDataConnectivityFilter> confilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	vtkSmartPointer<vtkCleanPolyData> cleanpolydata = vtkSmartPointer<vtkCleanPolyData>::New();

	confilter->SetInputConnection(ctr->GetOutputPort());
	confilter->SetExtractionModeToAllRegions();
	confilter->Update();

	unsigned int nreg = confilter->GetNumberOfExtractedRegions();

	for(unsigned int r = 0; r < nreg; r++)
	{
		std::cout<<" fn "<<isoval<<" nreg "<<r<<std::endl;
		confilter->InitializeSpecifiedRegionList();
		confilter->AddSpecifiedRegion(r);
		cleanpolydata->SetInputConnection(confilter->GetOutputPort());
		cleanpolydata->Update();
		vtkSmartPointer<vtkPolyData> polydata = cleanpolydata->GetOutput();
		ProcessContour(polydata);
	}
}

void Contours::GenerateIsoSpace(std::vector<float> & fvals)
{
	vtkSmartPointer<vtkContourFilter> ctr =	vtkSmartPointer<vtkContourFilter>::New();
	ctr->SetInputConnection(reader->GetOutputPort());
	ctr->ComputeNormalsOff();
	ctr->ComputeGradientsOff();
	ctr->ComputeScalarsOff();

	std::vector<float>::iterator fit = fvals.begin();
	float isoval = *fit;
	for(; fit != fvals.end(); fit++)
	{
		ctr->SetValue(0, isoval);
		ctr->Update();
		ProcessIsoSurface(isoval, ctr);
	}
}

Contours::Contours(const char* fname, vtkSmartPointer<vtkAppendPolyData> allcontours)
			: allcts(allcontours)
{
	reader = vtkStructuredPointsReader::New();
	reader->SetFileName(fname);
	reader->Update();
}

Contours::~Contours()
{
	reader->Delete();
}


void Contours::ComputeBD(vtkSmartPointer<vtkStructuredPoints> vtkstrpts)
{
	vtkSmartPointer<vtkDataArray> ps = vtkstrpts->GetPointData()->GetScalars();
	for(unsigned int i = 0; i < vtkstrpts->GetNumberOfPoints(); i++)
	{
		double p[3];
		double w = ps->GetComponent(i,0);
		vtkstrpts->GetPoint(i,p);
		verts.push_back(Vertex(p,w));
	}
	
	bd = new BD(verts);
	bd->BuildBD();
}

void Contours::ExtractSymmetry()
{
	double range[2];
	vtkSmartPointer<vtkStructuredPoints> vtkstrpts = reader->GetOutput();
	vtkstrpts->GetScalarRange(range);

	std::vector<float> fvals;
	for(unsigned int i = 1; i < 10; i++)
	{
		float isoval = range[0] + i*(range[1] - range[0])/10.0;
		//fvals.push_back(isoval);
	}

	fvals.push_back(0.987);
	ComputeBD(vtkstrpts);

	vtkSmartPointer<vtkIntArray> bidarray = vtkIntArray::New();
	bidarray->SetName("bids");
	std::vector<int> & bmap = bd->GetVertMap();
	bidarray->SetArray(&bmap[0], bmap.size(), 1);
	vtkstrpts->GetPointData()->AddArray(bidarray);
	vtkstrpts->Update();

	vtkSmartPointer<vtkStructuredPointsWriter> strpwriter =
		vtkSmartPointer<vtkStructuredPointsWriter>::New();
	strpwriter->SetInput(vtkstrpts);
	strpwriter->SetFileName("t.vtk");
	strpwriter->Update();

	GenerateIsoSpace(fvals);
}

