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
#include "Sampler.hpp"

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
int Contours::FindBranchId(vtkSmartPointer<vtkPolyData> contour, float isoval)
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
		assert(cid >= 0);
		vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
		vtkstrpts->GetCellPoints(cid, ptids);
		assert(8 == ptids->GetNumberOfIds());
		bid = 0;
		int inbid = 0, outbid = 0;
		for(unsigned int i = 0; i < 8; i++)
		{
			unsigned int v = ptids->GetId(i);
			unsigned int bid1 = bmap[v];
			if(bid == 0)
				bid = bid1;
			else if(bid != bid1)
				bid = -1;

			if(verts[v].w < isoval)
			{
				unsigned int inbid1 = bid1;
				if(inbid == 0)
					inbid = inbid1;
				else if(inbid != inbid1)
					inbid = -1;
			}
			if(verts[v].w > isoval)
			{
				unsigned int outbid1 = bid1;
				if(outbid == 0)
					outbid = outbid1;
				else if(outbid != outbid1)
					outbid = -1;
			}

			if(bid == -1 && inbid == -1 && outbid == -1)
			{
				break;
			}
		}
		if(bid > 0) 
			;
		else if(inbid > 1 && outbid > 1)
			continue;
		else if(inbid > 1)
		{
			bid = inbid;
		}
		else if(outbid > 1)
		{
			bid = outbid;
		}
				
		if(bid <= 1)
			bid = 0;
		else
			break;
	}
	if(bid < 1 ) bid = 1;
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
//		if(target > 0.95) target = 0.95;
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

		if(nreg == 1)
		{
			std::cout<<" fnid "<<fid<<" cid "<<cid<<" nreg "<<r<<" of "<<nreg<<std::endl;
			std::cout <<" size  = "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
			std::cout<<" Only 1 contour, not processing"<<std::endl;
		}
		else if(polydata->GetNumberOfPoints() > 20)
		{
			char fn[100];
			int bid = FindBranchId(polydata, isoval);
			std::cout<<" fnid "<<fid<<" cid "<<cid<<" nreg "<<r<<" of "<<nreg<<" bid "<<bid<<std::endl;
			std::cout <<" size  = "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
			if(bid == -1) 
			{
				sprintf(fn,"disbid-%d-%d.vtk",cid, fid);
				writer->SetFileName(fn);
				trifil->SetInput(polydata);
				writer->SetInputConnection(trifil->GetOutputPort());
				writer->Write();
				std::cout<<"Branch Id is -1!!!"<<std::endl;
				continue;
			}

			SymBranch* b = bd->GetBranch(bid);
			boost::unordered_map<unsigned int, unsigned int>::iterator cit;
			cit = b->comps.find(fid);
			if(cit != b->comps.end())
			{
				std::cout<<"Selecting contour "<<cid<<" from bid "<<bid<<" at fnid"<<fid<<" - "<<isoval<<std::endl;
				CompNode* c = new CompNode(cid++, bid, fid);
				b->comps[fid] = c->id;
				GenCompCords(c, polydata);
				c->csz = polydata->GetNumberOfPoints();
				compmgr->AddComp(c);
				sprintf(fn,"%d.vtk",c->id);
			}
			else
			{
				sprintf(fn,"%d-%d-dis.vtk",bid, fid);
				std::cout<<"Discarding contour for bid "<<bid<<" at "<<isoval<<std::endl;
			}
			writer->SetFileName(fn);
			trifil->SetInput(polydata);
			writer->SetInputConnection(trifil->GetOutputPort());
			writer->Write();

		}
		else
		{
			std::cout<<" fnid "<<fid<<" cid "<<cid<<" nreg "<<r<<" of "<<nreg<<std::endl;
			std::cout <<" size  = "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
			std::cout<<" Too few points, not processing"<<std::endl;

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


void Contours::Preprocess(vtkSmartPointer<vtkStructuredPoints> vtkstrpts)
{
	vtkSmartPointer<vtkDataArray> ps = vtkstrpts->GetPointData()->GetScalars();
	for(unsigned int i = 0; i < vtkstrpts->GetNumberOfPoints(); i++)
	{
		double p[3];
		double w = ps->GetComponent(i,0);
		vtkstrpts->GetPoint(i,p);
		verts.push_back(Vertex(p,w));
	}
	int dim[3];
	vtkstrpts->GetDimensions(dim);
	bd = new BD(verts, dim[0], dim[1], dim[2]);
	std::vector<unsigned int> sadidx;
	bd->BuildBD(sadidx);
	Sampler s(bd, sadidx, verts);
	s.Sample(fvals);
}

void Contours::ExtractSymmetry()
{
	double range[2];
	vtkSmartPointer<vtkStructuredPoints> vtkstrpts = reader->GetOutput();
	vtkstrpts->GetScalarRange(range);

	Preprocess(vtkstrpts);

	compmgr = new CompMgr(fvals, bd);
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

