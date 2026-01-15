/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-05 23:35:05
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-15 14:40:41
 * FilePath: \FVM\convection_diffusion\src\boundarysettings.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "boundarysettings.h"
#include <cmath>

void BoundarySettings::SetBcFacesForCells(StructureMesh &mesh)
{
    const int& dim = mesh.dim_export();
    const std::vector<std::vector<float>>& range = mesh.range_export();
    const std::vector<int>& cellnum = mesh.cellnum_export();
    const std::vector<float>& x_nodes = mesh.xnodes_export();
    const std::vector<float>& y_nodes = mesh.ynodes_export();
    const std::vector<float>& z_nodes = mesh.znodes_export();

    //存储单元边界条件id，若为内部点，则设置为0
    //顺序为 e w n s t b
    bcid = std::vector<std::vector<std::vector<std::vector<int>>>>(
        cellnum[0], std::vector<std::vector<std::vector<int>>>(cellnum[1],
        std::vector<std::vector<int>>(cellnum[2], std::vector<int>(2 * dim, 0))));

    //设置检查容差
    float eps = 1.0e-6f;
    // float x0, y0, z0;

    //判断每个网格的边界是否是在外边界上，并且记录到bcid中
    for(int k=0; k < cellnum[2]; ++k)
    {
        for(int j=0; j < cellnum[1]; ++j)
        {
            for(int i=0; i < cellnum[0]; ++i)
            {
                // 每个网格的东边界
                if (std::abs(x_nodes[i+1]-range[0][1]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                {   
                    bcid[i][j][k][fid_e] = bcid_xmax;
                }
                // 每个网格的西边界
                if (std::abs(x_nodes[i]-range[0][0]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                {   
                    bcid[i][j][k][fid_w] = bcid_xmin;
                }
                // 每个网格的北边界
                if (std::abs(y_nodes[j+1]-range[1][1]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                {   
                    bcid[i][j][k][fid_n] = bcid_ymax;
                }
                // 每个网格的南边界
                if (std::abs(y_nodes[j]-range[1][0]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                {   
                    bcid[i][j][k][fid_s] = bcid_ymin;
                }
                
                //如果维度为三维，要多加两个边界判断
                    // if dim == 3:
                   // 每个网格的底边界
                if(dim == 3)
                {
                    //每个网格的底边界
                    if (std::abs(z_nodes[k]-range[2][0]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                    {   
                        bcid[i][j][k][fid_b] = bcid_zmax;
                    }
                    //每个网格的顶边界
                    if (std::abs(z_nodes[k+1]-range[2][1]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                    {   
                        bcid[i][j][k][fid_t] = bcid_zmax;
                    }
                }
            }
        }
    }

    // std::cout << "---------- Check BC Interface Type ----------" << "\n";
    // for(int i = 0; i<4; ++i)
    // {
    //     std::cout << bcid[0][2][0][i] << "\t";
    // }
    // std::cout << "\n";
    // for(int i = 0; i<4; ++i)
    // {
    //     std::cout << bcid[1][2][0][i] << "\t";
    // }
    // std::cout << "\n";
    // for(int i = 0; i<4; ++i)
    // {
    //     std::cout << bcid[2][2][0][i] << "\t";
    // }
    // std::cout << "\n";
    // for(int i = 0; i<4; ++i)
    // {
    //     std::cout << bcid[3][2][0][i] << "\t";
    // }
    // std::cout << "\n";
}

void BoundarySettings::SetPhysicalBoundary(StructureMesh &mesh,
                                            const std::string &bc_xmin, 
                                            const std::string &bc_xmax, 
                                            const std::string &bc_ymin, 
                                            const std::string &bc_ymax,
                                            const std::string& bc_zmin,
                                            const std::string& bc_zmax)
{
    const int& dim = mesh.dim_export();
    bc_physics = std::vector<int>(2*dim + 1, 0); //+1是包括内部的点，[0]表示内部的点

    //内部边界条件
    bc_physics[bcid_none] = bc_none;

    //Xmin边界条件
    if(bc_xmin == "inlet")
    {
        bc_physics[bcid_xmin] = bc_inlet;
    }else if(bc_xmin == "outlet")
    {
        bc_physics[bcid_xmin] = bc_outlet;
    }else if(bc_xmin == "wall")
    {
        bc_physics[bcid_xmin] = bc_wall;
    }

    //Xmax边界条件
    if(bc_xmax == "inlet")
    {
        bc_physics[bcid_xmax] = bc_inlet;
    }else if(bc_xmax == "outlet")
    {
        bc_physics[bcid_xmax] = bc_outlet;
    }else if(bc_xmax == "wall")
    {
        bc_physics[bcid_xmax] = bc_wall;
    }

    //Ymin边界条件
    if(bc_ymin == "inlet")
    {
        bc_physics[bcid_ymin] = bc_inlet;
    }else if(bc_ymin == "outlet")
    {
        bc_physics[bcid_ymin] = bc_outlet;
    }else if(bc_ymin == "wall")
    {
        bc_physics[bcid_ymin] = bc_wall;
    }

    //Ymax边界条件
    if(bc_ymax == "inlet")
    {
        bc_physics[bcid_ymax] = bc_inlet;
    }else if(bc_ymax == "outlet")
    {
        bc_physics[bcid_ymax] = bc_outlet;
    }else if(bc_ymax == "wall")
    {
        bc_physics[bcid_ymax] = bc_wall;
    }

    if(dim == 3)
    {
        //Ymin边界条件
        if(bc_zmin == "inlet")
        {
            bc_physics[bcid_zmin] = bc_inlet;
        }else if(bc_zmin == "outlet")
        {
            bc_physics[bcid_zmin] = bc_outlet;
        }else if(bc_zmin == "wall")
        {
            bc_physics[bcid_zmin] = bc_wall;
        }

        //Ymin边界条件
        if(bc_zmax == "inlet")
        {
            bc_physics[bcid_zmax] = bc_inlet;
        }else if(bc_zmax == "outlet")
        {
            bc_physics[bcid_zmax] = bc_outlet;
        }else if(bc_zmax == "wall")
        {
            bc_physics[bcid_zmax] = bc_wall;
        }
    }
}

void BoundarySettings::SetNumericalBoundary(StructureMesh& mesh,
                                        std::string& xmin_type, float& xmin_value,
                                        std::string& xmax_type, float& xmax_value,
                                        std::string& ymin_type, float& ymin_value,
                                        std::string& ymax_type, float& ymax_value,
                                        std::string& zmin_type, float& zmin_value,
                                        std::string& zmax_type, float& zmax_value)
{
    // 在 SetThermalBoundary 开始处初始化大小
    const int& dim = mesh.dim_export();
    bc_temp = std::vector<temp>(2 * dim + 1); 

    //internal edges / faces
    bc_temp[bcid_none].bc_t_type = bc_none;
    bc_temp[bcid_none].t_dirichlet = 0.0f;
    bc_temp[bcid_none].t_neumann = 0.0f;

    //Xmin边界条件
    if(xmin_type == "constant")
    {
        bc_temp[bcid_xmin].bc_t_type = t_bc_type_dirichlet;
        bc_temp[bcid_xmin].t_dirichlet = xmin_value;
    }else if(xmin_type == "heat_flux")
    {
        bc_temp[bcid_xmin].bc_t_type = t_bc_type_neumann;
        bc_temp[bcid_xmin].t_neumann = xmin_value;
    }

    //Xmax边界条件
    if(xmax_type == "constant")
    {
        bc_temp[bcid_xmax].bc_t_type = t_bc_type_dirichlet;
        bc_temp[bcid_xmax].t_dirichlet = xmax_value;
    }else if(xmax_type == "heat_flux")
    {
        bc_temp[bcid_xmax].bc_t_type = t_bc_type_neumann;
        bc_temp[bcid_xmax].t_neumann = xmax_value;
    }

    //Ymin边界条件
    if(ymin_type == "constant")
    {
        bc_temp[bcid_ymin].bc_t_type = t_bc_type_dirichlet;
        bc_temp[bcid_ymin].t_dirichlet = ymin_value;
    }else if(ymin_type == "heat_flux")
    {
        bc_temp[bcid_ymin].bc_t_type = t_bc_type_neumann;
        bc_temp[bcid_ymin].t_neumann = ymin_value;
    }

    //Ymax边界条件
    if(ymax_type == "constant")
    {
        bc_temp[bcid_ymax].bc_t_type = t_bc_type_dirichlet;
        bc_temp[bcid_ymax].t_dirichlet = ymax_value;
    }else if(ymax_type == "heat_flux")
    {
        bc_temp[bcid_ymax].bc_t_type = t_bc_type_neumann;
        bc_temp[bcid_ymax].t_neumann = ymax_value;
    }

    if(dim == 3)
    {
        //Zmin边界条件
        if(zmin_type == "constant")
        {
            bc_temp[bcid_zmin].bc_t_type = t_bc_type_dirichlet;
            bc_temp[bcid_zmin].t_dirichlet = zmin_value;
        }else if(zmin_type == "heat_flux")
        {
            bc_temp[bcid_zmin].bc_t_type = t_bc_type_neumann;
            bc_temp[bcid_zmin].t_neumann = zmin_value;
        }

        //Zmax边界条件
        if(zmax_type == "constant")
        {
            bc_temp[bcid_zmax].bc_t_type = t_bc_type_dirichlet;
            bc_temp[bcid_zmax].t_dirichlet = zmax_value;
        }else if(zmax_type == "heat_flux")
        {
            bc_temp[bcid_zmax].bc_t_type = t_bc_type_neumann;
            bc_temp[bcid_zmax].t_neumann = zmax_value;
        }
    }   
}

void BoundarySettings::DisplayFluidType()
{
    std::cout << "---------- Fluid BC Type ----------" << std::endl;
    for(auto& i : bc_physics)
    {
        std::cout << i << "\t";
    }
    std::cout << "\n";
}
