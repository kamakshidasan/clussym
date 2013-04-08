#include <vtkSmartPointer.h>
#include <vtkStructuredPointsReader.h>
#include <vtkStructuredPoints.h>
#include <vtkMarchingCubes.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkPolyDataMapper.h>
#include <vtkCleanPolyData.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataWriter.h>

#include "LB.hpp"
#include "Cluster.hpp"

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
	vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
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
		std::cout <<isoval<< " Regions = "<<nreg<<std::endl;
		confilter->SetExtractionModeToSpecifiedRegions();
//		unsigned int r = 0;
		for(unsigned int r = 0; r < nreg; r++)
		{
			confilter->InitializeSpecifiedRegionList();
			confilter->AddSpecifiedRegion(r);
			vtkSmartPointer<vtkCleanPolyData> cleanpolydata = vtkSmartPointer<vtkCleanPolyData>::New();
			cleanpolydata->SetInputConnection(confilter->GetOutputPort());
			cleanpolydata->Update();
			vtkPolyData* polydata = cleanpolydata->GetOutput();
			std::cout <<"Surface "<<s<<" Region "<<r<<" size  = "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
			if(polydata->GetNumberOfPoints() > 20)
			{
				LB lb;
				lb.GetEigen(polydata, surfcords);
				char fn[100];
				sprintf(fn,"%d.vtk",nsurf);
				writer->SetFileName(fn);
				writer->SetInputConnection(cleanpolydata->GetOutputPort());
				writer->Write();
				nsurf++;
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
