/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-02-05 19:27:03
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-02-09 20:44:26
 * FilePath: \simple\src\solver.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "solver.h"
// #include "linearsolver.h"
#include "assemblesystem.h"

bool Solver::CollocatedSegregated(int &nlitera_nums, StructureMesh &mesh, 
    SolverSettings &solversettings, MaterialSettings &material, 
    BoundarySettings &boundary, PostProcess &postprocess)
{
    AssembleSystem assemblesystem;
    if(solversettings.eqn_type == solversettings.eqn_type_conduct_flow ||
        solversettings.eqn_type == solversettings.eqn_type_flow)
    {
        //计算momemntum linking coefs in all directions（所有方向上的动量连接系数）, 不考虑边界条件
        assemblesystem.MomentumCoefs(mesh, solversettings, material);
    }
    
    return false;
}
