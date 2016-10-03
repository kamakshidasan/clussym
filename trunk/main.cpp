#include "Contours.hpp"
#include <vtkPolyData.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyDataMapper.h>
#include <vtkAppendPolyData.h>

int main(int argc, char* argv[])
{
	if (argc < 3)
	{
		std::cerr << "Usage: " << argv[0] 
			<< " InputFile(.vtk) 0/1 (min/max)"<< std::endl;
		return EXIT_FAILURE;
	}

	vtkSmartPointer<vtkAppendPolyData> allcontours = vtkSmartPointer<vtkAppendPolyData>::New();
	Contours ct(argv[1],  argv[6] ? argv[6] : 0, allcontours);

	ct.ExtractSymmetry(atoi(argv[2]), atof(argv[3]), atof(argv[4]), atof(argv[5]), argv[6] ? 2 : 1);

	vtkSmartPointer<vtkPolyDataMapper> mapper =
		vtkSmartPointer<vtkPolyDataMapper>::New();
	allcontours->Update();
	mapper->SetInputData(allcontours->GetOutput());


	mapper->ScalarVisibilityOff();    // utilize actor's property I set

	// Visualize
	/*vtkSmartPointer<vtkActor> actor =
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
	*/
	return EXIT_SUCCESS;
}
