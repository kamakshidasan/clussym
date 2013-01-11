cmake_minimum_required(VERSION 2.8)
 
PROJECT(Contours)

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

add_executable(Contours Contours LB)

target_link_libraries(Contours ${VTK_LIBRARIES})

