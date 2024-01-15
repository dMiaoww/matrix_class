cmake_minimum_required(VERSION 2.8.3)
project(matrix)
add_definitions(-std=c++11)


find_package(Eigen3 REQUIRED)
include_directories("/usr/include/eigen3")

add_subdirectory(class1)
add_subdirectory(class2)
add_subdirectory(class3)