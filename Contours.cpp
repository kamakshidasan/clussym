#include "vtkSmartPointer.h"
#include "vtkStructuredPointsReader.h"
#include "vtkStructuredPoints.h"
#include "vtkMarchingCubes.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkPolyDataMapper.h"
#include "vtkCleanPolyData.h"
#include "vtkActor.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
 
int main(int argc, char* argv[])
{
  if (argc < 3)
    {
    std::cerr << "Usage: " << argv[0] 
              << " InputFile(.vtk) Threshold" << std::endl;
    return EXIT_FAILURE;
    }
  const char* fileName = argv[1];
  float threshold = atof(argv[2]);
  int extractLargest = 1;
  //if (argc == 4)
    {
    extractLargest = 0;//atoi(argv[3]);
    }

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
  mc->ComputeNormalsOn();
  mc->ComputeGradientsOn();

  vtkSmartPointer<vtkPolyDataConnectivityFilter> confilter =
	  vtkSmartPointer<vtkPolyDataConnectivityFilter>::New();
  // Create a 3D model using marching cubes
  unsigned int s = 3;
  //for(unsigned int s = 0; s < 10; s++)
  {
  	mc->SetValue(0, range[0] + s*(range[1] - range[0])/10.0);

	confilter->SetInputConnection(mc->GetOutputPort());
	vtkPolyData* currentComponent = NULL;
	confilter->SetExtractionModeToAllRegions();
	confilter->Update();
	unsigned int nreg = confilter->GetNumberOfExtractedRegions();
	std::cout << "Regions = "<<nreg<<std::endl;
	confilter->SetExtractionModeToSpecifiedRegions();
	unsigned int i = 26;
	//for(unsigned int i = 0; i < nreg; i++)
	{
		confilter->InitializeSpecifiedRegionList();
		confilter->AddSpecifiedRegion(i);
		confilter->Update();
		vtkSmartPointer<vtkCleanPolyData> cleanpolydata = vtkSmartPointer<vtkCleanPolyData>::New();
		cleanpolydata->SetInputConnection(confilter->GetOutputPort());
		cleanpolydata->Update();
		vtkPolyData* polydata = cleanpolydata->GetOutput();
		std::cout << "Surface "<<s<<" Region "<<i<<" size  = "<<polydata->GetNumberOfCells()<<" "<<polydata->GetNumberOfPolys()<<" "<<polydata->GetNumberOfPoints()<<std::endl;
	}
	
  }
  //  mc->GenerateValues(10, -0.9615, 2.644);  
 
  // Create a mapper
  vtkSmartPointer<vtkPolyDataMapper> mapper =
	  vtkSmartPointer<vtkPolyDataMapper>::New();
  mapper->SetInputConnection(confilter->GetOutputPort());

  mapper->ScalarVisibilityOff();    // utilize actor's property I set
 
  // Visualize
  vtkSmartPointer<vtkActor> actor =
    vtkSmartPointer<vtkActor>::New();
  actor->GetProperty()->SetColor(1,1,1);
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
