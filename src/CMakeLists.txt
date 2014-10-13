cmake_minimum_required (VERSION 2.6)
project (only_smooth_project)

find_package(Boost COMPONENTS system filesystem REQUIRED)

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11") 

add_library(misc misc.cpp)
add_library(gain_table gain_table.cpp misc)
add_library(smooth_interface SHARED gain_table misc smooth_interface.cpp)
target_link_libraries(smooth_interface ${Boost_LIBRARIES} gsl gslcblas)
add_executable(only_smooth only_smooth.cpp)

target_link_libraries(only_smooth smooth_interface gain_table misc ${Boost_LIBRARIES} gsl gslcblas)
