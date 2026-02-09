/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-06 23:33:07
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-02-09 20:57:35
 * FilePath: \simple\include\linearsolver.h
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

class LinearSolver
{
public:
    std::vector<float> CalScalarNorm2(StructureMesh& mesh, SolverSettings& solversettings, std::string var);
    void GaussSiedelMethod(int& nlitera_nums,std::vector<std::vector<std::vector<std::vector<float>>>>& cellcoef,
                            StructureMesh& mesh, PostProcess& postprocess, 
                            SolverSettings& solversettings, int& iter_num, int& total_iter_num, float& rc, float& torlence);
    void JacobiIteraMethod(int& nlitera_nums,std::vector<std::vector<std::vector<std::vector<float>>>>& cc,
                            StructureMesh& mesh, PostProcess& postprocess, 
                            SolverSettings& solversettings, int& iter_num, int& total_iter_num, float& rc, float& torlence);
    bool stop_sim{false};
};

#endif
