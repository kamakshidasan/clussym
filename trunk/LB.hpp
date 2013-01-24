#ifndef LB_HPP 
#define LB_HPP

#include <boost/unordered_map.hpp>
#include "Elements.hpp"
class LB
{
public:
	void FillMatrix(vtkPolyData*);
private:

	boost::unordered_map<Tri, unsigned int> trimap;


	double Cotangent(double v1[], double v2[]);
	void GetCots(vtkIdType *cpts, vtkPoints *pts, double cot[]);
	bool DuplicateCell(vtkIdType* cpts);
};

#endif
