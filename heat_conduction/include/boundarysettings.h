#ifndef BOUNDARYSETTINGS_H
#define BOUNDARYSETTINGS_H
#include <vector>
#include "structure_mesh.h"

class BoundarySettings
{
private:
    // None, Wall, Inlet, Outlet
    int btype{0};

    //初始化各个数据
    //每个网格的面，二维有四个面，三维有六个面，为这些面编号
    int fid_e = 0;
    int fid_w = 1;
    int fid_n = 2;
    int fid_s = 3;
    int fid_t = 4;
    int fid_b = 5;

    //物理边界条件类型，进出口和壁面以及空。这里的空指的是网格在计算域内部，所有的边都不和边界边相接触
    int bc_none   = 0;
    int bc_wall   = 1;
    int bc_inlet  = 2;
    int bc_outlet = 3;
 
    //为每个边界创建一个索引，上下 左右 前后
    std::vector<std::vector<std::vector<std::vector<int>>>> bcid; //四维数据，存储每个单元的的边界id
    int bcid_none   = 0; //这里表示非边界
    int bcid_xmin   = 1;
    int bcid_xmax   = 2;
    int bcid_ymin   = 3;
    int bcid_ymax   = 4;
    int bcid_zmin   = 5;
    int bcid_zmax   = 6;
    int num_bcs = 7;     //这样类型的个数

    //一个单元的边界条件类型和取值
    std::vector<int> bc_type_ocell;
    std::vector<int> bc_temp_type_ocell;

    //温度边界条件取值
    float t_e{0.0f};
    float t_w{0.0f};
    float t_n{0.0f};
    float t_s{0.0f};
    float t_t{0.0f};
    float t_b{0.0f};

    //温度边界条件类型
    int t_bc_type_dirichlet = 0;
    int t_bc_type_neumann = 1;

    //温度边界条件情况
    int t_btype{0};
    float dirichlet_t{0.0f};
    float neumann_t{0.0f};
public:
    void SetBcFacesForCells(StructureMesh& mesh);


};

#endif