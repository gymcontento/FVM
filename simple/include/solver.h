#ifndef SOLVER_H
#define SOLVER_H
#include "structure_mesh.h"
#include "solver_settings.h"
#include "material_settings.h"
#include "boundarysettings.h"
#include "postprocess.h"

class Solver
{
private:

public:
    bool CollocatedSegregated(int& nlitera_nums, StructureMesh &mesh, 
                            SolverSettings &solversettings, MaterialSettings &material, 
                            BoundarySettings &boundary, PostProcess& postprocess);

    bool stop_sim{false};
};

#endif