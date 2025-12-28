cmake_minimum_required(VERSION 3.10.0)

set(CMAKE_CXX_STANDARD 11)

project(1d_convection_diffusion)

add_executable(1d_convection_diffusion main.cpp 1d_cd.cpp)