cmake_minimum_required(VERSION 2.8.3)
project(matrix)
add_definitions(-std=c++11)


find_package(Eigen3 REQUIRED)
include_directories("/usr/include/eigen3")

add_subdirectory(class1)
add_subdirectory(class2)
add_subdirectory(class3)
add_subdirectory(class4)
add_subdirectory(class5)

add_executable(class6 class6.cpp)
target_link_libraries(class6 ${EIGEN_LIBRARIES})

add_executable(class8 class8.cpp)
target_link_libraries(class8 ${EIGEN_LIBRARIES})

add_executable(class9 class9.cpp)
target_link_libraries(class9 ${EIGEN_LIBRARIES})