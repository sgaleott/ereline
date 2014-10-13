cmake_minimum_required (VERSION 2.6)
project (only_smooth_project)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -I/usr/common/usg/gsl/1.13/include -L/usr/common/usg/gsl/1.13/lib") 

add_library(misc misc.cpp)
add_library(gain_table gain_table.cpp misc)
add_library(smooth_interface SHARED gain_table misc smooth_interface.cpp)
target_link_libraries(smooth_interface gsl gslcblas)
add_executable(only_smooth only_smooth.cpp)

target_link_libraries(only_smooth smooth_interface gain_table misc gsl gslcblas)
