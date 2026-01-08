/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-07 23:19:21
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-07 23:47:09
 * FilePath: \FVM\heat_conduction\include\assemblesystem.h
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#ifndef ASSEMBLESYSTEM_H
#define ASSEMBLESYSTEM_H
#include "structure_mesh.h"
#include "solver_settings.h"
#include "material_settings.h"
#include "boundarysettings.h"

class AssembleSystem
{
public:
    void ConductionCoefs(StructureMesh& mesh, SolverSettings& solver,
                        MaterialSettings& material);
    void ConductionCoefsBoud(StructureMesh& mesh, MaterialSettings& material,
                            BoundarySettings& boundary);
};

#endif