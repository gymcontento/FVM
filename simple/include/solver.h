/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-06 23:33:07
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-11 00:34:49
 * FilePath: \FVM\heat_conduction\include\solver.h
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#ifndef SOLVER_H
#define SOLVER_H
#include <vector>
#include "structure_mesh.h"
#include "postprocess.h"
#include "solver_settings.h"

class Solver
{
public:
    std::vector<float> CalScalarNorm2(StructureMesh& mesh, SolverSettings& solversettings, std::string var);
    void GaussSiedelMethod(int& nlitera_nums, StructureMesh& mesh, PostProcess& postprocess, SolverSettings& solversettings);
    void JacobiIteraMethod(int& nlitera_nums, StructureMesh& mesh, PostProcess& postprocess, SolverSettings& solversettings);
    bool stop_sim{false};
};

#endif
