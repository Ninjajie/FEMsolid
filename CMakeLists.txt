cmake_minimum_required(VERSION 3.9)
project(FEMsolid)

set(CMAKE_CXX_STANDARD 11)
# Find eigen in system directory and treat it as an interface library
find_package(Eigen3 REQUIRED)
add_library(eigen INTERFACE IMPORTED)

add_executable(FEMsolid
        main.cpp)
target_include_directories(FEMsolid SYSTEM PUBLIC ${EIGEN3_INCLUDE_DIR})
target_link_libraries(FEMsolid eigen)