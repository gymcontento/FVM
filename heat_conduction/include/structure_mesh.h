#ifndef STRUCTURE_MESH_H
#define STRUCTURE_MESH_H
#include <vector>
#include <string>
#include <iostream>

class StructureMesh
{
public:
    virtual ~StructureMesh() = default;

    void CreateMesh(int& dimensions, std::vector<int>& nums);
    void CreateCoordinates(std::vector<std::vector<float>>& irange);
    void CreateFieldMeshData();
    void SetInitialT(float& T);
    void CreateCoeffMeshData();
    void CreateSimulationData();

    std::vector<int> cellnum_export() {return cellnum;}
    std::vector<int> nodenum_export() {return nodenum;}

    std::vector<float> xnodes_export() {return x_nodes;}
    std::vector<float> ynodes_export() {return y_nodes;}
    std::vector<float> znodes_export() {return z_nodes;}

    std::vector<float> xcells_export() {return x_cells;}
    std::vector<float> ycells_export() {return y_cells;}
    std::vector<float> zcells_export() {return z_cells;}

    std::vector<std::vector<std::vector<float>>> tfield_export() {return t;}

private:
    int dim;
    std::vector<int> cellnum{0,0,0};
    std::vector<int> nodenum{0,0,0};
    std::vector<std::vector<float>> range = {
        {0.0f, 0.0f}, {0.0f, 0.0f}, {0.0f, 0.0f}
    };

    std::vector<float> x_nodes;
    std::vector<float> y_nodes;
    std::vector<float> z_nodes;

    std::vector<float> x_cells;
    std::vector<float> y_cells;
    std::vector<float> z_cells;

    float dx;
    float dy;
    float dz;

    std::vector<std::vector<std::vector<float>>> t;
    std::vector<std::vector<std::vector<float>>> t0;

    int numcoef;
    int id_ap, id_ae, id_aw, id_an, id_as, id_aT, id_aB, id_bsrc;
    std::vector<std::vector<std::vector<std::vector<int>>>> cellcoeff;
};

#endif