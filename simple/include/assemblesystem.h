/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-07 23:19:21
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-02-06 23:28:19
 * FilePath: \simple\include\assemblesystem.h
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
    void MomentumCoefs(StructureMesh& mesh, SolverSettings& solversettings,
                        MaterialSettings& material);
    void ConductionCoefs(StructureMesh& mesh, SolverSettings& solversettings,
                        MaterialSettings& material);
    void ConductionCoefsBoud(StructureMesh& mesh, MaterialSettings& material,
                            BoundarySettings& boundary);
};

#endif
