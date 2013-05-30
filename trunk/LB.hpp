#ifndef LB_HPP 
#define LB_HPP

#include <boost/unordered_map.hpp>

#include "Elements.hpp"
class LB
{
public:
	void GetEigen(vtkPolyData*, std::vector<std::vector<double> > & surfcords);
private:
	double Cotangent(double v1[], double v2[], double & area);
	void GetVorArea(double cot[3], double len[3], double triarea, double vorarea[3]);
	void GetCotsLensArea(vtkIdType *cpts, vtkPoints *pts, double cot[], double len[], double & a);
	bool DuplicateCell(vtkIdType* cpts, boost::unordered_map<Tri, unsigned int> & trimap);
	bool DuplicateEdge(vtkIdType cpts0, vtkIdType cpts1, boost::unordered_map<Edge, unsigned int> & edmap);
};

#endif
