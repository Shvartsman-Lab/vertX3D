cmake_minimum_required(VERSION 3.5.1)
project(vertX3D)

set(CMAKE_CXX_FLAGS "-std=c++11 -w")

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(include)

# Add the executable
add_executable(vertX3D_spreading_flat "examples/spreading_flat.cpp")
add_executable(vertX3D_spreading_sphere "examples/spreading_sphere.cpp")
add_executable(vertX3D_homeostasis "examples/homeostasis.cpp")
add_executable(vertX3D_fluidization "examples/fluidization.cpp")
add_executable(vertX3D_sphere "examples/sphere.cpp")
