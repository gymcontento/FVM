/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-06 23:33:07
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-06 23:39:43
 * FilePath: \FVM\heat_conduction\include\solver.h
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#ifndef SOLVER_H
#define SOLVER_H
#include <vector>
#include "structure_mesh.h"

class Solver
{
public:
    std::vector<float> CalScalarNorm2(StructureMesh& mesh, 
                                    const float& u0, const float& u,
                                    std::string& var);
    void ScalarGaussSiedel();
};

#endif