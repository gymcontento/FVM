/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-03 00:42:12
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-22 09:39:52
 * FilePath: \simple\include\structure_mesh.h
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
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
    void SetInitialVelocity(const std::vector<float>& initialvelocity);
    void SetInitialFaceVelocity(const std::vector<float>& facevelocity);
    void SetInitialPressure(const float& initialpressure);
    void SetSourceTerm(const std::vector<float>& sourceterm);
    void SetHeatSource(const std::vector<float>& sourceterm);
    
    // void SetInitialUF(float& u);
    // void SetInitialVF(float& v);
    // void SetInitialWF(float& w);
    void CreateCoeffMeshData();

    //Display Info
    const int& dim_export() {return dim;}
    const std::vector<int>& cellnum_export() {return cellnum;}
    const std::vector<int>& nodenum_export() {return nodenum;}
    const std::vector<std::vector<float>>& range_export() {return range;}

    const std::vector<float>& xnodes_export() {return x_nodes;}
    const std::vector<float>& ynodes_export() {return y_nodes;}
    const std::vector<float>& znodes_export() {return z_nodes;}

    const std::vector<float>& xcells_export() {return x_cells;}
    const std::vector<float>& ycells_export() {return y_cells;}
    const std::vector<float>& zcells_export() {return z_cells;}

    const float& dx_export() {return dx;};
    const float& dy_export() {return dy;};
    const float& dz_export() {return dz;};

    const std::vector<std::vector<std::vector<float>>>& tfield_export() {return t;}
    const std::vector<std::vector<std::vector<float>>>& uffield_export() {return uf;}
    const std::vector<std::vector<std::vector<float>>>& vffield_export() {return vf;}
    const std::vector<std::vector<std::vector<float>>>& wffield_export() {return wf;}

    std::vector<std::vector<std::vector<std::vector<float>>>> cellcoeff_u;
    std::vector<std::vector<std::vector<std::vector<float>>>> cellcoeff_v;
    std::vector<std::vector<std::vector<std::vector<float>>>> cellcoeff_w;
    std::vector<std::vector<std::vector<std::vector<float>>>> cellcoeff_p;
    std::vector<std::vector<std::vector<std::vector<float>>>> cellcoeff_t;
    int id_aP, id_aE, id_aW, id_aN, id_aS, id_aT, id_aB, id_bsrc;
    int numcoef;
    int dim;
    //温度场和速度场（单元界面）
    //速度场、压力场和温度场、保守力场
    std::vector<std::vector<std::vector<float>>> t;
    std::vector<std::vector<std::vector<float>>> u;
    std::vector<std::vector<std::vector<float>>> v;
    std::vector<std::vector<std::vector<float>>> w;
    std::vector<std::vector<std::vector<float>>> p;
    std::vector<std::vector<std::vector<float>>> x_source;
    std::vector<std::vector<std::vector<float>>> y_source;
    std::vector<std::vector<std::vector<float>>> z_source;  
    //旧速度、压力和温度场
    std::vector<std::vector<std::vector<float>>> t0;
    std::vector<std::vector<std::vector<float>>> u0;
    std::vector<std::vector<std::vector<float>>> v0;
    std::vector<std::vector<std::vector<float>>> w0;
    std::vector<std::vector<std::vector<float>>> p0;
    //面上速度和压力校正场
    std::vector<std::vector<std::vector<float>>> pp;
    std::vector<std::vector<std::vector<float>>> uf;
    std::vector<std::vector<std::vector<float>>> vf;
    std::vector<std::vector<std::vector<float>>> wf;
    //动量和温度方程的源项
    std::vector<std::vector<std::vector<float>>> ru;
    std::vector<std::vector<std::vector<float>>> rv;
    std::vector<std::vector<std::vector<float>>> rw;
    std::vector<std::vector<std::vector<float>>> rp;
    std::vector<std::vector<std::vector<float>>> rt;
private:
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
};

#endif