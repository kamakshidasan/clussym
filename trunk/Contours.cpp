#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkCleanPolyData.h>
#include <vtkPolyDataNormals.h>
#include <vtkTriangleFilter.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataWriter.h>
#include <vtkPLYWriter.h>
#include <vtkWindowedSincPolyDataFilter.h>
#include <vtkDelaunay2D.h>
#include <vtkDecimatePro.h>
#include <vtkDoubleArray.h>
#include <vtkPointData.h>
#include <vtkLinearSubdivisionFilter.h>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Point_with_normal_3.h>

#include <deque>

#include "LB.hpp"
#include "Cluster.hpp"
#include "Remesh.hpp"


typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;

void vtktoPointList(PointList & pl, vtkPolyData* mesh)
{
	vtkSmartPointer<vtkDoubleArray> normals = 
		vtkDoubleArray::SafeDownCast(mesh->GetPointData()->GetNormals());
	
	assert(mesh->GetPointData());
	assert(mesh->GetPointData()->GetNumberOfTuples() == mesh->GetNumberOfPoints());

	for(vtkIdType i = 0; i < mesh->GetNumberOfPoints(); i++) 
	{ 
		double p[3]; 
		mesh->GetPoint(i, p);
		Point pt(p[0], p[1], p[2]);
		mesh->GetPointData()->GetTuple(i, p);
		Vector n(p[0], p[1], p[2]);
		pl.push_back(Point_with_normal(pt,n));
	} 

}
int main(int argc, char* argv[])
{
	if (argc < 2)
	{
		std::cerr << "Usage: " << argv[0] 
			<< " InputFile(.vtk)"<< std::endl;
		return EXIT_FAILURE;
	}
	const char* fileName = argv[1];

	// Load data
	vtkSmartPointer<vtkStructuredPointsReader> reader =
		vtkSmartPointer<vtkStructuredPointsReader>::New();
	reader->SetFileName(fileName);
	reader->Update();
	double range[2];
	vtkStructuredPoints* vtkstrpts = reader->GetOutput();
	vtkstrpts->GetScalarRange(range);
	vtkSmartPointer<vtkMarchingCubes> mc =
		vtkSmartPointer<vtkMarchingCubes>::New();
	mc->SetInputConnection(reader->GetOutputPort());
	mc->ComputeNormalsOff();
	mc->ComputeGradientsOff();
	mc->ComputeScalarsOff();

	vtkSmartPointer<vtkPolyDataConnectivityFilter> confilter =
		vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
	// Create a 3D model using marching cubes
	vtkSmartPointer<vtkDecimatePro> decimate = vtkSmartPointer<vtkDecimatePro>::New();

	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
	vtkSmartPointer<vtkPLYWriter> PLYwriter = vtkSmartPointer<vtkPLYWriter>::New();
	vtkSmartPointer<vtkDelaunay2D> delaunay = vtkSmartPointer<vtkDelaunay2D>::New();

	  vtkSmartPointer<vtkWindowedSincPolyDataFilter> smoother =
		      vtkSmartPointer<vtkWindowedSincPolyDataFilter>::New();
	  unsigned int smoothingIterations = 50;
	  double passBand = 0.001;
	  double featureAngle = 120.0;
	std::vector<std::vector<double> > surfcords;
	float val[] = {.84,.86,.88,.9};
	unsigned int nsurf = 0;
	for(unsigned int s = 0; s < 4; s++)
	{
		float isoval = range[0] + val[s]*(range[1] - range[0])/9.0;
		mc->SetValue(0, val[s]);
		//mc->SetValue(0, isoval);
		confilter->SetInputConnection(mc->GetOutputPort());
		confilter->SetExtractionModeToAllRegions();
		confilter->Update();
		unsigned int nreg = confilter->GetNumberOfExtractedRegions();
		std::cout << "Regions = "<<nreg<<std::endl;
		confilter->SetExtractionModeToSpecifiedRegions();
//		unsigned int r = 0;
		for(unsigned int r = 0; r < nreg; r++)
		{
			confilter->InitializeSpecifiedRegionList();
			confilter->AddSpecifiedRegion(r);
			vtkSmartPointer<vtkCleanPolyData> cleanpolydata = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanpolydata->SetInputConnection(confilter->GetOutputPort());
			cleanpolydata->Update();
			vtkSmartPointer<vtkPolyData> polydata = cleanpolydata->GetOutput();


/*			smoother->SetInputConnection(cleanpolydata->GetOutputPort());
			smoother->SetNumberOfIterations(smoothingIterations);
			smoother->BoundarySmoothingOff();
			smoother->FeatureEdgeSmoothingOn();
			smoother->SetFeatureAngle(featureAngle);
			smoother->SetPassBand(passBand);
			smoother->NonManifoldSmoothingOn();
			smoother->NormalizeCoordinatesOn();
			smoother->Update();

//			delaunay->SetInput(polydata);
//			delaunay->SetSource(polydata);
//			delaunay->Update();
			polydata = smoother->GetOutput();	*/
			std::cout <<"Surface "<<s<<" Region "<<r<<" size  = "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
			unsigned int ntri = polydata->GetNumberOfPolys();
			if(ntri > 300)
			{
				decimate->SetInputConnection(polydata->GetProducerPort());
				float target = 1 - 300.0/ntri;
				decimate->SetTargetReduction(target);
				decimate->Update();
//				polydata->ShallowCopy(decimate->GetOutput());
				polydata = decimate->GetOutput();				
				std::cout <<"After decimate "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
			}
			if(polydata->GetNumberOfPoints() > 20)
			{
				vtkSmartPointer<vtkLinearSubdivisionFilter> subdivisionFilter = 
					vtkSmartPointer<vtkLinearSubdivisionFilter>::New();
				subdivisionFilter->SetNumberOfSubdivisions(2);
				subdivisionFilter->SetInputConnection(polydata->GetProducerPort());
				polydata = subdivisionFilter->GetOutput();
				PointList pl;
				vtkSmartPointer<vtkPolyDataNormals> normalGenerator = vtkSmartPointer<vtkPolyDataNormals>::New();
				normalGenerator->SetInput(polydata);
				normalGenerator->ComputePointNormalsOn();
				normalGenerator->ComputeCellNormalsOff();
				normalGenerator->SplittingOff();
				normalGenerator->NonManifoldTraversalOff();
				normalGenerator->Update();

				polydata = normalGenerator->GetOutput();
				vtktoPointList(pl, polydata);
				vtkPolyData* newpoly = Remesh(pl, nsurf);
				//PointListtovtk();
				LB lb;
				lb.GetEigen(newpoly, surfcords);
				char fn[100];
				sprintf(fn,"%d.vtk",nsurf);
				writer->SetFileName(fn);
				vtkSmartPointer<vtkTriangleFilter> trifil = vtkSmartPointer<vtkTriangleFilter>::New();
				trifil->SetInput(newpoly);
				writer->SetInputConnection(trifil->GetOutputPort());
				writer->Write();
				nsurf++;
				newpoly->Delete();
			}
		}

	}
	//  mc->GenerateValues(10, -0.9615, 2.644);  

	// Create a mapper
	Cluster cl(surfcords);
	cl.GetClusters();
	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputConnection(mc->GetOutputPort());

	mapper->ScalarVisibilityOff();    // utilize actor's property I set

	// Visualize
	vtkSmartPointer<vtkActor> actor =
		vtkSmartPointer<vtkActor>::New();
	actor->GetProperty()->SetColor(1,1,1);
	actor->GetProperty()->SetInterpolationToFlat();
	actor->GetProperty()->SetEdgeColor(1,0,0);
	actor->GetProperty()->EdgeVisibilityOn();
	actor->SetMapper(mapper);

	vtkSmartPointer<vtkRenderer> renderer =
		vtkSmartPointer<vtkRenderer>::New();
	renderer->AddActor(actor);

	vtkSmartPointer<vtkRenderWindow> renwin =
		vtkSmartPointer<vtkRenderWindow>::New();
	renwin->AddRenderer(renderer);

	vtkSmartPointer<vtkRenderWindowInteractor> iren =
		vtkSmartPointer<vtkRenderWindowInteractor>::New();
	iren->SetRenderWindow(renwin);
	iren->Initialize();
	iren->Start();

	return EXIT_SUCCESS;
}
