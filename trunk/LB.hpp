#ifndef LB_HPP 
#define LB_HPP

#include <boost/unordered_map.hpp>

#include "Elements.hpp"
class LB
{
public:
	void FillMatrix(vtkPolyData*);
private:
	double Cotangent(double v1[], double v2[]);
	void GetCots(vtkIdType *cpts, vtkPoints *pts, double cot[]);
	bool DuplicateCell(vtkIdType* cpts, boost::unordered_map<Tri, unsigned int> & trimap);
	bool DuplicateEdge(vtkIdType cpts0, vtkIdType cpts1, float val, boost::unordered_map<Edge, unsigned int> & edmap);
};

#endif
