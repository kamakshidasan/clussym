cmake_minimum_required(VERSION 2.8)
 
PROJECT(Contours)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})
add_definitions(-std=c++11)
add_executable(Contours Contours LB Decompose)

include_directories("/home/ibis/Downloads/arpack++/include")
include_directories("/home/ibis/Downloads/FLENS-2012-10-01/")

target_link_libraries(Contours ${VTK_LIBRARIES} arpack superlu) 

