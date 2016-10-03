#include "Export.hpp"
#include "BD.hpp"
#include <vtkPointData.h>
#include <vtkDoubleArray.h>
#include <vtkFloatArray.h>

Export::Export(vtkSmartPointer<vtkStructuredPoints> vtkstrpts) : orgstrpts(vtkstrpts)
{
	strpwriter = vtkSmartPointer<vtkStructuredPointsWriter>::New();
	newstrpts = vtkSmartPointer<vtkStructuredPoints>::New();
	newstrpts->DeepCopy(orgstrpts);
}
void Export::CreateVolume()
{
	mask = std::vector<unsigned int>(orgstrpts->GetNumberOfPoints(), 0);
	newstrpts->GetPointData()->SetScalars(orgstrpts->GetPointData()->GetScalars());
}

void Export::AddToVolume(unsigned int cid)
{
	bd->AppendExportMask(cid, mask);
}

void Export::FlushVolume(const std::string & fn, unsigned int maskf)
{
	double range[2];
	double val = range[maskf];
	orgstrpts->GetScalarRange(range);
	vtkSmartPointer<vtkDataArray> f = newstrpts->GetPointData()->GetScalars();
	for(unsigned int i = 0; i < newstrpts->GetNumberOfPoints(); i++)
	{
		if(mask[i])
			f->SetComponent(i,0,val);
	}

	//newstrpts->GetPointData()->SetScalars(f);
	strpwriter->SetInputData(newstrpts);
	strpwriter->SetFileName(fn.c_str());
}
