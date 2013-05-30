cmake_minimum_required(VERSION 2.8)
 
PROJECT(Contours)
add_definitions(-DCGAL_EIGEN3_ENABLED)

find_package(VTK REQUIRED)
find_package(CGAL REQUIRED)
include(${VTK_USE_FILE})
include(${CGAL_USE_FILE})


add_executable(Contours Contours LB Decompose Cluster Remesh Utils BD)

include_directories("/usr/include/eigen3/")
include_directories("/home/ibis/Downloads/arpack++/include")
include_directories("/home/ibis/Downloads/FLENS-2012-10-01/")

target_link_libraries(Contours ${VTK_LIBRARIES} tourtre arpack superlu ANN ) 

