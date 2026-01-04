#ifndef POSTPROCESS_H
#define POSTPROCESS_H
#include <string>
#include "structure_mesh.h"

class PostProcess
{
private:
    std::vector<int> cellnum;
    std::vector<int> nodenum;
    std::vector<float> x_nodes;
    std::vector<float> y_nodes;
    std::vector<float> z_nodes;

    std::vector<float> x_cells;
    std::vector<float> y_cells;
    std::vector<float> z_cells;

    std::vector<std::vector<std::vector<float>>> t;
public:
    void WriteVTKCollocated_temp(
        StructureMesh& mesh,std::string filename = "post_temp.vtk");
    

};

#endif