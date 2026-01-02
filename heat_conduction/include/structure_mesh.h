#ifndef STRUCTURE_MESH_H
#define STRUCTURE_MESH_H
#include <vector>
#include <string>
#include <iostream>

class StructureMesh
{
public:
    virtual ~StructureMesh() = default;
    void CreateMesh(int& dimensions, std::vector<int>& nums)
    {
        // Number of Spatial dimensions
        dim = dimensions;

        // Number of cells in each direction
        cellnum[0] = nums[0];
        cellnum[1] = nums[1];
        cellnum[2] = nums[2];

        // Number of nodes in each direction
        nodenum[0] = nums[0] + 1;
        nodenum[1] = nums[1] + 1;
        nodenum[2] = nums[2] + 1;

        std::cout << "Number of cells in x, y, z = " << cellnum[0] << ", "
            << cellnum[1] << ", " << cellnum[2] << std::endl;
        std::cout << "Number of nodes in x, y, z = " << nodenum[0] << ", "
            << nodenum[1] << ", " << nodenum[2] << std::endl;
    }
private:
    int dim;
    std::vector<int> cellnum{0,0,0};
    std::vector<int> nodenum{0,0,0};
};

#endif