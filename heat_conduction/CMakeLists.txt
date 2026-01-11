cmake_minimum_required(VERSION 3.10.0)

set(CMAKE_CXX_STANDARD 17)

project(demo)

include_directories(include)

add_executable(demo src/main.cpp 
                    src/structure_mesh.cpp 
                    src/postprocess.cpp
                    src/solver_settings.cpp
                    src/material_settings.cpp
                    src/boundarysettings.cpp
                    src/kernel.cpp
                    src/assemblesystem.cpp
                    src/solver.cpp)