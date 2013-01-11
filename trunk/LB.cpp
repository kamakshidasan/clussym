#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkSmartPointer.h>

#include "LB.hpp"


void LB::FillMatrix(vtkPolyData* mesh)
{
	vtkCellArray* polys = mesh->GetPolys();
	unsigned int n = polys->GetNumberOfCells();
	polys->InitTraversal();
	vtkIdType npts, *pts;
	while(polys->GetNextCell(npts, pts))
	{
	}
}
