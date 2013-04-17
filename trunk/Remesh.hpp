#ifndef REMESH_HPP
#define REMESH_HPP

#include <vector>

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef CGAL::Point_with_normal_3<Kernel> Point_with_normal;
typedef std::vector<Point_with_normal> PointList;

vtkPolyData* Remesh(PointList & points, int id);

#endif
