#include <vtkStructuredPointsWriter.h>
#include <vtkStructuredPoints.h>
#include <vtkUnstructuredGridWriter.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkUnstructuredGrid.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataReader.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkPLYWriter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkDecimatePro.h>
#include <vtkGeometryFilter.h>
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
#include <vtkDataSetTriangleFilter.h>
#include <vtkImageGaussianSmooth.h>
#include <vtkImageToStructuredPoints.h>
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

#include <sys/time.h>

struct timeval clus_start, clus_end;
struct timeval ct_start, ct_end;
struct timeval dis_start, dis_end;
double distime;
struct timeval ctr_start, ctr_end;
double ctrtime;
struct timeval tot_start, tot_end;
struct timeval oh_start, oh_end;
double ohtime;
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
int Contours::FindBranchId(vtkSmartPointer<vtkPolyData> contour, float isoval, unsigned int did)
{
	vtkSmartPointer<vtkCellArray> cells = contour->GetPolys();
	cells->InitTraversal();
	vtkIdType* cellpts;
	vtkIdType  numpts;
	//vtkSmartPointer<vtkStructuredPoints> vtkstrpts = reader->GetOutput();
	std::vector<int> & bmap = bd[did]->GetVertMap();
	int bid = -1, ret1 = 0;
	if(cells->GetNextCell(numpts, cellpts))
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
		unsigned int cid = tgrid[did]->FindCell(centre,0,0,0.001,subid,pcoords,weights);
		//unsigned int cid = vtkstrpts->FindCell(centre,0,0,0.001,subid,pcoords,weights);
		assert(cid >= 0);
		vtkSmartPointer<vtkIdList> ptids = vtkSmartPointer<vtkIdList>::New();
		//vtkstrpts->GetCellPoints(cid, ptids);
		tgrid[did]->GetCellPoints(cid, ptids);
		assert(8 == ptids->GetNumberOfIds());
		bid = 0;
		int inbid = 0, outbid = 0;
		for(unsigned int i = 0; i < 8; i++)
		{
			unsigned int v = ptids->GetId(i);
			unsigned int bid1 = bmap[v];
			unsigned int done = 0;
			while(!done)
			{
				if(bid1 == 1) done = 1;
				else
				{
					unsigned int sad = bd[did]->bridsarr[bid1]->sad;
					unsigned int ext = bd[did]->bridsarr[bid1]->ext;
					if(verts[did][sad].w > verts[did][ext].w)
					{
						if(isoval >= verts[did][ext].w && isoval <= verts[did][sad].w)
							done = 1;
						else
							bid1 = bd[did]->bridsarr[bid1]->par->bid;
					}
					else
					{
						if(isoval <= verts[did][ext].w && isoval >= verts[did][sad].w)
							done = 1;
						else
							bid1 = bd[did]->bridsarr[bid1]->par->bid;
					}

				}
			}

			if(bid1 != 1)
			{
				unsigned int sad = bd[did]->bridsarr[bid1]->sad;
				unsigned int ext = bd[did]->bridsarr[bid1]->ext;
				if(verts[did][sad].w > verts[did][ext].w)
				{
					if(verts[did][v].w <= isoval && verts[did][sad].w >= isoval)
					{
						if(bid == 0)
							bid = bid1;
						else if(bid != bid1)
						{
							bid = -1;
							//break;
						}
					}
				}
				else 
				{
					if(verts[did][v].w >= isoval && verts[did][sad].w <= isoval)
					{
						if(bid == 0)
							bid = bid1;
						else if(bid != bid1)
						{
							bid = -1;
							//break;
						}
					}
				}
			}
			else
			{
				ret1 = 1;
			}
			/*if(verts[v].w < isoval)
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
			}*/
		}
		//if(bid > 1) 
		//	break;
		/*else if(inbid > 1 && outbid > 1)
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
			break;*/
	}
	if(bid < 1 ) 
		if(ret1) bid = 1;
		else bid = -1;
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
//		std::cout <<"After decimate "<<contour->GetNumberOfCells()
//			<<" "<<contour->GetNumberOfPolys()
//			<<" "<<contour->GetNumberOfPoints()<<std::endl;
	}

	LB lb;
	lb.GetEigen(contour, c->cords);
	allcts->AddInputConnection(contour->GetProducerPort());

	//mesh->Delete();
}

void Contours::SetChildComps(CompNode* c, float curf, float prevf)
{
	unsigned int did = c->did;
	SymBranch* b = bd[did]->GetBranch(c->bid);
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
		float sadw = verts[did][(*bit)->sad].w ;
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

void Contours::ProcessIsoSurface(unsigned int fid, unsigned int prev, vtkSmartPointer<vtkContourFilter> ctr, unsigned int did)
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
		vtkSmartPointer<vtkPolyData> polydata1 = confilter->GetOutput();
		//std::cout<<"begin cleaning"<<std::endl;
		gettimeofday(&oh_start, NULL);
		cleanpolydata->SetInputConnection(confilter->GetOutputPort());
		cleanpolydata->Update();
		vtkSmartPointer<vtkPolyData> polydata = cleanpolydata->GetOutput();
		gettimeofday(&oh_end, NULL);
		ohtime += (oh_end.tv_sec-oh_start.tv_sec+oh_end.tv_usec/1000000.0-oh_start.tv_usec/1000000.0);
		//std::cout<<"end cleaning"<<std::endl;

		if(nreg == 1)
		{
//			std::cout<<" fnid "<<fid<<" did "<<did<<" cid "<<cid<<" nreg "<<r<<" of "<<nreg<<std::endl;
//			std::cout <<" size  = "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
//			std::cout<<" Only 1 contour, not processing"<<std::endl;
		}
		else if(polydata->GetNumberOfPoints() > 20)
		{
			char fn[100];
			int bid = FindBranchId(polydata, isoval, did);
		//	std::cout<<" fnid "<<fid<<" did "<<did<<" cid "<<cid<<" nreg "<<r<<" of "<<nreg<<" bid "<<bid<<std::endl;
		//	std::cout <<" size  = "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
			if(bid == -1) 
			{
/*				sprintf(fn,"disbid-%d-%d.vtk",cid, fid);
				writer->SetFileName(fn);
				trifil->SetInput(polydata);
				writer->SetInputConnection(trifil->GetOutputPort());
				writer->Write();*/
//				std::cout<<"Branch Id is -1!!!"<<std::endl;
				//assert(0);
				//bid = 1;
				continue;
			}
			SymBranch* b = bd[did]->GetBranch(bid);
			boost::unordered_map<unsigned int, unsigned int>::iterator cit;
			cit = b->comps.find(fid);
			if(cit != b->comps.end() && cit->second == -1)
			//if(cit == b->comps.end())
			{
				
//				std::cout<<"Selecting contour "<<cid<<" from bid "<<bid<<" at fnid"<<fid<<" - "<<isoval<<std::endl;
				CompNode* c = new CompNode(cid++, bid, fid);
				b->comps[fid] = c->id;
				gettimeofday(&dis_start, NULL);
				GenCompCords(c, polydata);
				gettimeofday(&dis_end, NULL);
				distime += (dis_end.tv_sec-dis_start.tv_sec+dis_end.tv_usec/1000000.0-dis_start.tv_usec/1000000.0);
				c->csz = polydata->GetNumberOfPoints();
				c->did = did;
				compmgr->AddComp(c);
				sprintf(fn,"%d-%d.vtk",c->did, c->id);
				writer->SetFileName(fn);
				trifil->SetInput(polydata);
				writer->SetInputConnection(trifil->GetOutputPort());
				writer->Write();
			}
			else
			{
//				sprintf(fn,"%d-%d-dis.vtk",bid, fid);
//				std::cout<<"Discarding contour for bid "<<bid<<" at "<<isoval<<std::endl;
			}
//			std::cout<<" ext "<<verts[did][bd[did]->bridsarr[bid]->ext].w<<" sad "<<verts[did][bd[did]->bridsarr[bid]->sad].w<<std::endl;

		}
		else
		{
/*			std::cout<<" fnid "<<fid<<" did "<<did<<" cid "<<cid<<" nreg "<<r<<" of "<<nreg<<std::endl;
			std::cout <<" size  = "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
			std::cout<<" Too few points, not processing"<<std::endl;*/

		}
	}
}

void Contours::GenerateIsoSpace(unsigned int did)
{
	vtkSmartPointer<vtkContourFilter> ctr =	vtkSmartPointer<vtkContourFilter>::New();
	ctr->SetInput(tgrid[did]);
	//ctr->SetInputConnection(reader->GetOutputPort());
	ctr->ComputeNormalsOff();
	ctr->ComputeGradientsOff();
	ctr->ComputeScalarsOff();

	int i = 0;
	unsigned int prev = 0;
	gettimeofday(&ctr_start, NULL);
	for(; i < fvals.size(); i++)
	{
//		std::cout<<"GenIso iso"<<fvals[i]<<std::endl;
		ctr->SetValue(0, fvals[i]);
		ctr->Update();
		ProcessIsoSurface(i, prev, ctr, did);
		prev = i;
	}
	gettimeofday(&ctr_end, NULL);
}

Contours::Contours(const char* fname1, const char* fname2, vtkSmartPointer<vtkAppendPolyData> allcontours)
			: allcts(allcontours), cid(0)
{
	vtkSmartPointer<vtkStructuredPointsReader> reader1 = vtkStructuredPointsReader::New();
	reader1->SetFileName(fname1);
	reader1->Update();
	trifil = vtkSmartPointer<vtkTriangleFilter>::New();
	writer = vtkSmartPointer<vtkPolyDataWriter>::New();


	//vtkSmartPointer<vtkDataSetTriangleFilter> vtktet1 = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
	//vtktet1->SetInputConnection(reader1->GetOutputPort());
	//tgrid.push_back(vtktet1->GetOutput());
	tgrid.push_back(reader1->GetOutput());
	tgrid[0]->Update();

	if(fname2)
	{
		vtkSmartPointer<vtkStructuredPointsReader> reader2 = vtkStructuredPointsReader::New();
		reader2->SetFileName(fname2);
		reader2->Update();
		//vtkSmartPointer<vtkDataSetTriangleFilter> vtktet2 = vtkSmartPointer<vtkDataSetTriangleFilter>::New();
		//vtktet2->SetInputConnection(reader2->GetOutputPort());
		//tgrid.push_back(vtktet2->GetOutput());
		tgrid.push_back(reader2->GetOutput());
		tgrid[1]->Update();

	}
}

Contours::~Contours()
{
}


//void Contours::Preprocess(vtkSmartPointer<vtkUnstructuredGrid> tgrid, unsigned int inv, unsigned int did)
void Contours::Preprocess(vtkSmartPointer<vtkStructuredPoints> & tgrid, unsigned int inv, unsigned int did, float alpha, float delta)
{
	/*vtkSmartPointer<vtkImageGaussianSmooth> gaussianSmoother =
		vtkSmartPointer<vtkImageGaussianSmooth>::New();

	gaussianSmoother->SetInput(tgrid);
	gaussianSmoother->SetStandardDeviations(1,1,1);
	gaussianSmoother->Update();
	vtkSmartPointer<vtkImageData> img = gaussianSmoother->GetOutput();
	vtkSmartPointer<vtkImageToStructuredPoints> i2sp =
	     vtkSmartPointer<vtkImageToStructuredPoints>::New();
	i2sp->SetInput(img);
	i2sp->Update();
	tgrid = i2sp->GetStructuredPointsOutput();
	*/
	vtkSmartPointer<vtkDataArray> ps = tgrid->GetPointData()->GetScalars();
	for(unsigned int i = 0; i < tgrid->GetNumberOfPoints(); i++)
	{
		double p[3];
		double w = ps->GetComponent(i,0);
		if(inv)
		{
			w = -w;
			ps->SetComponent(i,0, w);
		}
		tgrid->GetPoint(i,p);
		verts[did].push_back(Vertex(p,w));
	}
	int dim[3];
	tgrid->GetDimensions(dim);
	bd.push_back(new BD(verts[did], dim[0], dim[1], dim[2]));
	std::vector<unsigned int> sadidx;
	std::vector<float> tempf;
	gettimeofday(&ct_start, NULL);
	bd[did]->BuildBD(sadidx,tempf, delta);

	Sampler s(bd[did], sadidx, verts[did]);
	if(did == 0)
	{
		s.Sample(fvals, alpha);
	}
	else
	{
		float max = verts[0][bd[0]->bridsarr[1]->ext].w;
		float min = verts[0][bd[0]->bridsarr[1]->sad].w;
		s.PropogateValues(fvals, min, max);
	}
	gettimeofday(&ct_end, NULL);
}

void Contours::ExtractSymmetry(unsigned int inv, float epsd, float alpha, float delta, unsigned int dcnt)
{
	/*vtkSmartPointer<vtkUnstructuredGridReader> prdr1 = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	vtkSmartPointer<vtkUnstructuredGridReader> prdr2 = vtkSmartPointer<vtkUnstructuredGridReader>::New();
	vtkSmartPointer<vtkPolyData> pd1 = vtkSmartPointer<vtkPolyData>::New();
	vtkSmartPointer<vtkPolyData> pd2 = vtkSmartPointer<vtkPolyData>::New();
	prdr1->SetFileName("elcup.vtk");
//	prdr1->SetFileName("isacup.vtk");
	prdr1->Update();
	prdr2->Update();

	vtkGeometryFilter * geometryFilter =  vtkGeometryFilter::New(); 
	geometryFilter->SetInput(prdr1->GetOutput()); 
	geometryFilter->Update(); 

	pd1= geometryFilter->GetOutput(); 
	std::cout << "Output has " << pd1->GetNumberOfPoints() << " points." << std::endl; 

	vtkSmartPointer<vtkTriangleFilter> trifil = vtkSmartPointer<vtkTriangleFilter>::New();
	trifil->SetInput(pd1);
	trifil->Update();
	pd2 = trifil->GetOutput();

	std::cout << "Output has " << pd2->GetNumberOfPoints() << " points." << std::endl; 
	vtkSmartPointer<vtkQuadricDecimation> decimate = vtkSmartPointer<vtkQuadricDecimation>::New();
	decimate->SetInputConnection(pd2->GetProducerPort());
		float target = 0.1;
		decimate->SetTargetReduction(target);
		decimate->Update();
		pd2 = decimate->GetOutput();				
		std::cout <<"After decimate "<<pd2->GetNumberOfCells()
			<<" "<<pd2->GetNumberOfPolys()
			<<" "<<pd2->GetNumberOfPoints()<<std::endl;
	std::vector<float> c1;
	std::vector<float> c2;
	LB lb1;
	lb1.GetEigen(pd2, c1);
	exit(0);*/
	struct timeval timeval_start, timeval_end;
	gettimeofday(&timeval_start, NULL);
	compmgr = new CompMgr(bd);
	verts = std::vector<std::vector<Vertex> >(dcnt);
	for(unsigned int did = 0; did < dcnt; did++)
	{
		Preprocess(tgrid[did], inv, did, alpha, delta);
		if(did == 0)
			compmgr->Init(fvals);
		GenerateIsoSpace(did);
	}
	compmgr->ClusterComps(epsd, bd[0]->bridsarr.size());
	//compmgr->DistanceList();
	gettimeofday(&timeval_end, NULL);
	double time_start = ct_start.tv_sec + (double) ct_start.tv_usec/1000000;
	double time_end= ct_end.tv_sec + (double) ct_end.tv_usec/1000000;
	printf("Time: cttime %f\n", time_end - time_start);
	printf("Time: distime %f\n", distime);
	time_start = clus_start.tv_sec + (double) clus_start.tv_usec/1000000;
	time_end= clus_end.tv_sec + (double) clus_end.tv_usec/1000000;
	printf("Time: clustime %f\n", time_end - time_start);
	printf("Time: ohtime %f\n", ohtime);
	time_start = ctr_start.tv_sec + (double) ctr_start.tv_usec/1000000;
	time_end= ctr_end.tv_sec + (double) ctr_end.tv_usec/1000000;
	printf("Time: ctrtime %f\n", time_end-time_start);
	time_start = timeval_start.tv_sec + (double) timeval_start.tv_usec/1000000;
	time_end= timeval_end.tv_sec + (double) timeval_end.tv_usec/1000000;
	printf("Time: tottime %f\n", time_end-time_start);
	fflush(stdout);
}

