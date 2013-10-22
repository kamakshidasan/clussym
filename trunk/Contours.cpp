#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPLYWriter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkDecimatePro.h>
#include <vtkQuadricDecimation.h>
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

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Point_with_normal_3.h>
//#include "Remesh.hpp"

#include <deque>

#include "BD.hpp"
#include "LB.hpp"
#include "Cluster.hpp"
#include "Contours.hpp"
#include "CompMgr.hpp"


/*typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
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

}*/
int Contours::FindBranchId(vtkSmartPointer<vtkPolyData> contour)
{
	vtkSmartPointer<vtkCellArray> cells = contour->GetPolys();
	cells->InitTraversal();
	vtkIdType* cellpts;
	vtkIdType  numpts;
	vtkSmartPointer<vtkStructuredPoints> vtkstrpts = reader->GetOutput();
	std::vector<int> & bmap = bd->GetVertMap();
	int bid = -1;
	while(cells->GetNextCell(numpts, cellpts))
	{
		double centre[3], a[3], b[3], c[3];

		assert(numpts == 3);

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
		unsigned int cid = vtkstrpts->FindCell(centre,0,0,0.001,subid,pcoords,weights);
		//printf("centre (%f %f %f) cid %d\n",centre[0], centre[1], centre[2],cid);
		assert(cid >= 0);
		vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
		vtkstrpts->GetCellPoints(cid, ptids);
		assert(8 == ptids->GetNumberOfIds());
		bid = bmap[ptids->GetId(0)];
		bool iddiff = false;
		for(unsigned int i = 1; i < 8; i++)
		{
			unsigned int bid1 = bmap[ptids->GetId(i)];
			if(bid != bid1)
			{
				iddiff = true;
				bid = -1;
				break;
			}
		}
		if(iddiff == false) break;
	}
	return bid;
}
void Contours::GenCompCords(CompNode* c, vtkSmartPointer<vtkPolyData> contour)
{
	vtkSmartPointer<vtkQuadricDecimation> decimate = vtkSmartPointer<vtkQuadricDecimation>::New();
	//vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();
	vtkSmartPointer<vtkPolyDataNormals> gennorms = vtkSmartPointer<vtkPolyDataNormals>::New();
	vtkSmartPointer<vtkCleanPolyData> cleanpolydata = vtkSmartPointer<vtkCleanPolyData>::New();
	/*PointList pl;
	gennorms->SetInput(contour);
	gennorms->ComputePointNormalsOn();
	gennorms->ComputeCellNormalsOff();
	gennorms->SplittingOff();
	gennorms->NonManifoldTraversalOff();
	gennorms->Update();

	contour = gennorms->GetOutput();
	vtktoPointList(pl, contour);
	vtkPolyData* mesh = Remesh(pl);

	cleanpolydata->SetInput(mesh);*/

	cleanpolydata->SetInput(contour);
	cleanpolydata->Update();
	contour = cleanpolydata->GetOutput();
	unsigned int ntri = contour->GetNumberOfPolys();
	if(ntri > 2000)
	{
		decimate->SetInputConnection(contour->GetProducerPort());
		float target = 1 - 2000.0/ntri;
		if(target > 0.95) target = 0.95;
		decimate->SetTargetReduction(target);
		decimate->Update();
		contour = decimate->GetOutput();				
		std::cout <<"After decimate "<<contour->GetNumberOfCells()
			<<" "<<contour->GetNumberOfPolys()
			<<" "<<contour->GetNumberOfPoints()<<std::endl;
	}

	LB lb;
	lb.GetEigen(contour, c->cords);
	allcts->AddInputConnection(contour->GetProducerPort());

	//mesh->Delete();
}

void Contours::SetChildComps(CompNode* c, float curf, float prevf)
{
	SymBranch* b = bd->GetBranch(c->bid);
	boost::unordered_map<unsigned int, CompNode*>::iterator cit;
	cit = topcomps.find(c->bid);
	if(cit != topcomps.end())
	{
		cit->second->par = c;
		c->ch.push_back(cit->second);
	}

	std::list<SymBranch*>::iterator bit = b->ch.begin();
	for(; bit != b->ch.end(); bit++)
	{
		float sadw = verts[(*bit)->sad].w ;
		if(sadw >= prevf && sadw < curf)
		{
			cit = topcomps.find((*bit)->bid);
			assert(cit != topcomps.end());
			cit->second->par = c;
			c->ch.push_back(cit->second);
		}
	}

	topcomps[c->bid] = c;

}

void Contours::ProcessIsoSurface(unsigned int fid, unsigned int prev, vtkSmartPointer<vtkContourFilter> ctr)
{
	vtkSmartPointer<vtkPolyDataConnectivityFilter> confilter = vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	vtkSmartPointer<vtkCleanPolyData> cleanpolydata = vtkSmartPointer<vtkCleanPolyData>::New();

	confilter->SetInputConnection(ctr->GetOutputPort());
	confilter->SetExtractionModeToAllRegions();
	confilter->Update();

	unsigned int nreg = confilter->GetNumberOfExtractedRegions();
	confilter->SetExtractionModeToSpecifiedRegions();

	float isoval = fvals[fid];

	for(unsigned int r = 0; r < nreg; r++)
	{
		confilter->InitializeSpecifiedRegionList();
		confilter->AddSpecifiedRegion(r);
		cleanpolydata->SetInputConnection(confilter->GetOutputPort());
		cleanpolydata->Update();
		vtkSmartPointer<vtkPolyData> polydata = cleanpolydata->GetOutput();

		if(polydata->GetNumberOfPoints() > 20)
		{
			int bid = FindBranchId(polydata);
			if(bid == -1) 
			{
				std::cout<<"Branch Id is -1!!!"<<std::endl;
				continue;
			}

			std::cout<<" fnid "<<fid<<"cid,  nreg "<<cid<<" of "<<nreg<<" bid "<<bid<<std::endl;
			std::cout <<" size  = "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
			CompNode* c = new CompNode(cid++, bid, fid);
			SymBranch* b = bd->GetBranch(c->bid);
			boost::unordered_map<unsigned int, unsigned int>::iterator cit;
			cit = b->comps.find(fid);
			assert(cit == b->comps.end());
			b->comps[fid] = c->id;
			GenCompCords(c, polydata);
			c->csz = polydata->GetNumberOfPoints();
			compmgr->AddComp(c);
			char fn[100];
			sprintf(fn,"%d.vtk",c->id);
			writer->SetFileName(fn);
			trifil->SetInput(polydata);
			writer->SetInputConnection(trifil->GetOutputPort());
			writer->Write();

		}
	}
}

void Contours::GenerateIsoSpace()
{
	vtkSmartPointer<vtkContourFilter> ctr =	vtkSmartPointer<vtkContourFilter>::New();
	ctr->SetInputConnection(reader->GetOutputPort());
	ctr->ComputeNormalsOff();
	ctr->ComputeGradientsOff();
	ctr->ComputeScalarsOff();

	unsigned int i = 0;
	unsigned int prev = 0;
	for(; i < fvals.size(); i++)
	{
		ctr->SetValue(0, fvals[i]);
		ctr->Update();
		ProcessIsoSurface(i, prev, ctr);
		prev = i;
	}
	compmgr->ClusterComps();
}

Contours::Contours(const char* fname, vtkSmartPointer<vtkAppendPolyData> allcontours)
			: allcts(allcontours), cid(0)
{
	reader = vtkStructuredPointsReader::New();
	reader->SetFileName(fname);
	reader->Update();
	trifil = vtkSmartPointer<vtkTriangleFilter>::New();
	writer = vtkSmartPointer<vtkPolyDataWriter>::New();
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

	for(unsigned int i = 1; i < 3; i++)
	{
		float isoval = range[0] + i*(range[1] - range[0])/10.0;
		printf("isoval %f\n",isoval);
		fvals.push_back(isoval);
	}

//	fvals.push_back(-0.024);
//	fvals.push_back(-0.021);
	
	ComputeBD(vtkstrpts);
	compmgr = new CompMgr(fvals.size(), bd);

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

	GenerateIsoSpace();
}

