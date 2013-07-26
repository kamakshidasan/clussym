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
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>
#include <vtkIdTypeArray.h>
#include <vtkPointData.h>
#include <vtkThreshold.h>
#include <vtkVolume.h>
#include <vtkContourFilter.h>
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
/*
void ComputeBD(vtkStructuredPoints *vtkstrpts, std::vector<Vertex> & vlist)
{
	unsigned int n = vtkstrpts->GetNumberOfPoints();
	vtkDataArray* pscl = vtkstrpts->GetPointData()->GetScalars(); 
	for(unsigned int i = 0; i < n; i++)
	{
		double p[3], f;
		vtkstrpts->GetPoint(i, p);
		f = pscl->GetComponent(i, 0);
		Vertex v(p,f);
		vlist.push_back(v);	
	}
	BD bd(vlist);
	 
	vtkSmartPointer<vtkIntArray> brlbl =
		vtkSmartPointer<vtkIntArray>::New();
	brlbl->SetName("Branch");
	brlbl->SetNumberOfComponents(1);
	brlbl->SetNumberOfValues(vlist.size());
	bd.BuildBD(brlbl);
	
	vtkstrpts->GetPointData()->SetScalars(brlbl);
}
*/

void Contours::GenerateContour(vtkSmartPointer<vtkExtractSelection> extr, SymBranch* curbr, float isoval)
{

	vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
	vtkSmartPointer<vtkThreshold> thr = vtkSmartPointer<vtkThreshold>::New();
	vtkSmartPointer<vtkContourFilter> ctr =	vtkSmartPointer<vtkContourFilter>::New();
	vtkSmartPointer<vtkPolyDataConnectivityFilter> confilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	vtkSmartPointer<vtkCleanPolyData> cleanpolydata = vtkSmartPointer<vtkCleanPolyData>::New();
	vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkSmartPointer<vtkTriangleFilter> trifil = vtkSmartPointer<vtkTriangleFilter>::New();

//	thr->SetInputConnection(extr->GetOutputPort());
//	thr->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "bids");
//	thr->AllScalarsOff();
//	thr->ThresholdBetween(curbr->bid, curbr->bid);

	ctr->SetInputConnection(extr->GetOutputPort());
	ctr->ComputeNormalsOff();
	ctr->ComputeGradientsOff();
	ctr->ComputeScalarsOff();
	ctr->SetValue(0, isoval);
	ctr->Update();

	cleanpolydata->SetInputConnection(ctr->GetOutputPort());
	cleanpolydata->Update();
	vtkPolyData* polydata = cleanpolydata->GetOutput();

	confilter->SetInputConnection(ctr->GetOutputPort());
	confilter->SetExtractionModeToAllRegions();
	confilter->Update();

	unsigned int nreg = confilter->GetNumberOfExtractedRegions();

	std::cout <<"brid "<<curbr->bid<<" fn "<<isoval<<" nreg "<<nreg<<" size  = "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
//	if(nreg == 1)
	{
		if(polydata->GetNumberOfPoints() > 20)
		{
			PointList pl;
			normalGenerator->SetInput(polydata);
			normalGenerator->ComputePointNormalsOn();
			normalGenerator->ComputeCellNormalsOff();
			normalGenerator->SplittingOff();
			normalGenerator->NonManifoldTraversalOff();
			normalGenerator->Update();

			polydata = normalGenerator->GetOutput();
			vtktoPointList(pl, polydata);
			//	vtkSmartPointer<vtkPolyData> newpoly = polydata;
			/*vtkPolyData* newpoly = Remesh(pl);

			cleanpolydata->SetInput(newpoly);
			cleanpolydata->Update();
			polydata = cleanpolydata->GetOutput();

			char fn[100];
			sprintf(fn,"%dct.vtk",1);
			writer->SetFileName(fn);
			trifil->SetInput(polydata);
			writer->SetInputConnection(trifil->GetOutputPort());
			writer->Write();

			unsigned int ntri = polydata->GetNumberOfPolys();
			if(ntri > 20000)
			{
				decimate->SetInputConnection(polydata->GetProducerPort());
				float target = 1 - 20000.0/ntri;
				if(target > 0.8) target = 0.8;
				decimate->SetTargetReduction(target);
				decimate->Update();
				//polydata->ShallowCopy(decimate->GetOutput());
				polydata = decimate->GetOutput();				
					std::cout <<"After decimate "<<polydata->GetNumberOfCells()<<" "
						<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
				}

				std::vector<std::vector<double> > surfcords;
				LB lb;
				lb.GetEigen(polydata, surfcords);*/
				allcts->AddInputConnection(polydata->GetProducerPort());
				/*newpoly->Delete();

			}*/

		}
	}
}

void Contours::IntersectingCells(vtkUnstructuredGrid* ugrid, unsigned int brid, vtkIdTypeArray* cids, float isoval)
{

	unsigned int ncells = ugrid->GetNumberOfCells();
	unsigned int i = 0;
	vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
	vtkFloatArray* pscalars = vtkFloatArray::SafeDownCast(ugrid->GetPointData()->GetScalars());
	vtkIntArray* bidarray = vtkIntArray::SafeDownCast(ugrid->GetPointData()->GetArray("bids"));
	for(; i < ncells; i++)
	{
		ugrid->GetCellPoints(i, ptids);
		unsigned int npts = ptids->GetNumberOfIds();
		assert(npts == 8);
		unsigned int crossing[8];
		unsigned int brids[8];
		bool one = false;
		bool zero = false;
		for(vtkIdType j = 0; j < npts; j++)
		{
			unsigned int id = ptids->GetId(j);
			float fval = pscalars->GetValue(id);
			brids[j] = bidarray->GetValue(id);
			if(fval <= isoval)
			{
				zero = true;
				crossing[j] = 0;
			}
			else
			{
				one = true;
				crossing[j] = 1;
			}
		}
		if(one && zero)
		{
			if(((crossing[0] != crossing[1]) && (brids[0] == brid && brids[1] == brid))||
			   ((crossing[0] != crossing[2]) && (brids[0] == brid && brids[2] == brid))||
			   ((crossing[0] != crossing[4]) && (brids[0] == brid && brids[4] == brid))||
			   ((crossing[3] != crossing[1]) && (brids[3] == brid && brids[1] == brid))||
			   ((crossing[3] != crossing[2]) && (brids[3] == brid && brids[2] == brid))||
			   ((crossing[3] != crossing[7]) && (brids[3] == brid && brids[7] == brid))||
			   ((crossing[5] != crossing[1]) && (brids[5] == brid && brids[1] == brid))||
			   ((crossing[5] != crossing[4]) && (brids[5] == brid && brids[4] == brid))||
			   ((crossing[5] != crossing[7]) && (brids[5] == brid && brids[7] == brid))||
			   ((crossing[6] != crossing[2]) && (brids[6] == brid && brids[2] == brid))||
			   ((crossing[6] != crossing[4]) && (brids[6] == brid && brids[4] == brid))||
			   ((crossing[6] != crossing[7]) && (brids[6] == brid && brids[7] == brid)))
			   {
				cids->InsertNextValue(i);
			   }
		}


	}
/*

	std::list<SymBranch*>::iterator it = curbr->ch.begin();
	assert(!bidarray[curbr->bid]);
	bidarray[curbr->bid] = 1;
	std::cout<<" genbid "<<curbr->bid<<endl;
	for(; it != curbr->ch.end(); it++)
	{
		if((verts[(*it)->ext].w > verts[(*it)->sad].w && verts[(*it)->ext].w < isoval) ||
				(verts[(*it)->ext].w < verts[(*it)->sad].w && verts[(*it)->sad].w < isoval))
		{
			GenerateGridIds(*it, bidarray, isoval);
		}

	}*/
}
void Contours::GenerateIsoSpace(SymBranch* curbr, std::vector<float> & fvals)
{
	unsigned int brid = curbr->bid;
	vtkSmartPointer<vtkThreshold> thr = vtkSmartPointer<vtkThreshold>::New();
	thr->SetInputConnection(reader->GetOutputPort());
	thr->SetInputArrayToProcess(0, 0, 0, vtkDataObject::FIELD_ASSOCIATION_POINTS, "bids");
	thr->AllScalarsOff();
	thr->ThresholdBetween(brid, brid);
	thr->Update();
	
	vtkSmartPointer<vtkExtractSelection> extr = vtkSmartPointer<vtkExtractSelection>::New();
	vtkUnstructuredGrid* ugrid = thr->GetOutput();
	extr->SetInput(0, ugrid);
	vtkSmartPointer<vtkUnstructuredGridWriter> strpwriter =
					vtkSmartPointer<vtkUnstructuredGridWriter>::New();
	strpwriter->SetInput(thr->GetOutput());
	strpwriter->SetFileName("u.vtk");
	strpwriter->Update();

	std::vector<float>::iterator fit = fvals.begin();
	for(; fit != fvals.end(); fit++)
	{
		float isoval = *fit;
		if((isoval > verts[curbr->ext].w && isoval < verts[curbr->sad].w) ||
				(isoval < verts[curbr->ext].w && isoval > verts[curbr->sad].w))
		{
			std::cout<<"For bid "<<brid<<endl;
			vtkSmartPointer<vtkIdTypeArray> cids = vtkSmartPointer<vtkIdTypeArray>::New();
			IntersectingCells(ugrid, brid, cids, isoval);
			cids->SetNumberOfComponents(1);
			vtkSmartPointer<vtkSelection> select = vtkSmartPointer<vtkSelection>::New();
			vtkSmartPointer<vtkSelectionNode> node= vtkSmartPointer<vtkSelectionNode>::New();
			node->SetFieldType(vtkSelectionNode::CELL);
			node->SetContentType(vtkSelectionNode::INDICES);
			node->SetSelectionList(cids);
			select->AddNode(node);
			extr->SetInput(1, select);
			extr->Update();
			vtkSmartPointer<vtkUnstructuredGridWriter> strpwriter =
				vtkSmartPointer<vtkUnstructuredGridWriter>::New();
			strpwriter->SetInput(extr->GetOutput());
			strpwriter->SetFileName("v.vtk");
			strpwriter->Update();
			GenerateContour(extr, curbr, isoval);
		}
	}
	std::list<SymBranch*>::iterator brit = curbr->ch.begin();
	for(; brit != curbr->ch.end(); brit++)
	{
		//		GenerateIsoSpace(*brit, fvals, extr);
	}
}

Contours::Contours(const char* fname, vtkSmartPointer<vtkAppendPolyData> & allcontours)
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

	GenerateIsoSpace(bd->bridsarr[1], fvals);
}

