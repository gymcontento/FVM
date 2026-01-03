cmake_minimum_required(VERSION 3.10.0)

set(CMAKE_CXX_STANDARD 17)

project(demo)

include_directories(include)

add_executable(demo src/main.cpp src/structure_mesh.cpp)