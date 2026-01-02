#include "structure_mesh.h"
#include <vector>

int main()
{
    int dim = 2;
    std::vector<int> nums{2,2,2};
    StructureMesh sm;
    sm.CreateMesh(dim,nums); 
}