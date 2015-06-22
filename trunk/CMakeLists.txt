cmake_minimum_required(VERSION 2.8)
 
PROJECT(Contours)
#add_definitions(-DCGAL_EIGEN3_ENABLED)

find_package(VTK REQUIRED)
find_package(CGAL REQUIRED)
include(${VTK_USE_FILE})
include(${CGAL_USE_FILE})


add_executable(sym main Contours LB Decompose Cluster Utils BD CompMgr Export Sampler)

include_directories("/usr/include/eigen3/")

target_link_libraries(sym ${VTK_LIBRARIES} tourtre arpack superlu ANN ) 

