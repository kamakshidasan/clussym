#ifndef LB_HPP 
#define LB_HPP

class LB
{
public:
	void FillMatrix(vtkPolyData*);
private:
	double Cotangent(double v1[], double v2[]);
	void GetCots(vtkIdType *cpts, vtkPoints *pts, double cot[]);
	bool DuplicateCell(vtkIdType* cpts);
};

#endif
