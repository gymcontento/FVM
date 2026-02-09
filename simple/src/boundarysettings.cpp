/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-05 23:35:05
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-02-03 14:55:08
 * FilePath: \simple\src\boundarysettings.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "boundarysettings.h"
#include <cmath>
#include <iomanip>

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
                    bcid[i][j][k][fid_b] = bcid_zmin;
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
        //Zmin边界条件
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

        //Zmax边界条件
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

void BoundarySettings::SetNumericalBoundary(StructureMesh& mesh, std::vector<std::vector<std::string>>& numerical_bc_type,
                            std::vector<std::vector<float>>& numerical_bc_value) 
{
    // 在 SetThermalBoundary 开始处初始化大小
    const int& dim = mesh.dim_export();
    bc_vel = std::vector<vel>(2 * dim + 1); 

    //internal edges / faces
    bc_vel[bcid_none].u_type = bc_none;
    bc_vel[bcid_none].v_type = bc_none;
    bc_vel[bcid_none].u_dirichlet = 0.0f;
    bc_vel[bcid_none].v_dirichlet = 0.0f;
    bc_vel[bcid_none].u_neumann = 0.0f;
    bc_vel[bcid_none].v_neumann = 0.0f;
    if(dim == 3){
        bc_vel[bcid_none].w_type = bc_none;
        bc_vel[bcid_none].w_dirichlet = 0.0f;
        bc_vel[bcid_none].w_neumann = 0.0f;
    }

    //Xmin边界条件
    //设置 u 方向速度的边界条件
    if(numerical_bc_type[0][0] == "constant")
    {
        bc_vel[bcid_xmin].u_type = bc_type_dirichlet;
        bc_vel[bcid_xmin].u_dirichlet = numerical_bc_value[0][0];
    }else if(numerical_bc_type[0][0] == "flux")
    {
        bc_vel[bcid_xmin].u_type = bc_type_neumann;
        bc_vel[bcid_xmin].u_neumann = numerical_bc_value[0][0];
    }
    //设置 v 方向速度的边界条件
    if(numerical_bc_type[0][1] == "constant")
    {
        bc_vel[bcid_xmin].v_type = bc_type_dirichlet;
        bc_vel[bcid_xmin].v_dirichlet = numerical_bc_value[0][1];
    }else if(numerical_bc_type[0][1] == "flux")
    {
        bc_vel[bcid_xmin].v_type = bc_type_neumann;
        bc_vel[bcid_xmin].v_neumann = numerical_bc_value[0][1];
    }
    //设置 w 方向速度的边界条件
    if(dim == 3){
        if(numerical_bc_type[0][2] == "constant")
        {
            bc_vel[bcid_xmin].w_type = bc_type_dirichlet;
            bc_vel[bcid_xmin].w_dirichlet = numerical_bc_value[0][2];
        }else if(numerical_bc_type[0][2] == "flux")
        {
            bc_vel[bcid_xmin].w_type = bc_type_neumann;
            bc_vel[bcid_xmin].w_neumann = numerical_bc_value[0][2];
        } 
    }

    //Xmax边界条件
    //设置 u 方向速度的边界条件
    if(numerical_bc_type[1][0] == "constant")
    {
        bc_vel[bcid_xmax].u_type = bc_type_dirichlet;
        bc_vel[bcid_xmax].u_dirichlet = numerical_bc_value[1][0];
    }else if(numerical_bc_type[1][0] == "flux")
    {
        bc_vel[bcid_xmax].u_type = bc_type_neumann;
        bc_vel[bcid_xmax].u_neumann = numerical_bc_value[1][0];
    }
    //设置 v 方向速度的边界条件
    if(numerical_bc_type[1][1] == "constant")
    {
        bc_vel[bcid_xmax].v_type = bc_type_dirichlet;
        bc_vel[bcid_xmax].v_dirichlet = numerical_bc_value[1][1];
    }else if(numerical_bc_type[1][1] == "flux")
    {
        bc_vel[bcid_xmax].v_type = bc_type_neumann;
        bc_vel[bcid_xmax].v_neumann = numerical_bc_value[1][1];
    }
    //设置 w 方向速度的边界条件
    if(dim == 3){
        if(numerical_bc_type[1][2] == "constant")
        {
            bc_vel[bcid_xmax].w_type = bc_type_dirichlet;
            bc_vel[bcid_xmax].w_dirichlet = numerical_bc_value[1][2];
        }else if(numerical_bc_type[1][2] == "flux")
        {
            bc_vel[bcid_xmax].w_type = bc_type_neumann;
            bc_vel[bcid_xmax].w_neumann = numerical_bc_value[1][2];
        } 
    }

    //Ymin边界条件
    //设置 u 方向速度的边界条件
    if(numerical_bc_type[2][0] == "constant")
    {
        bc_vel[bcid_ymin].u_type = bc_type_dirichlet;
        bc_vel[bcid_ymin].u_dirichlet = numerical_bc_value[2][0];
    }else if(numerical_bc_type[2][0] == "flux")
    {
        bc_vel[bcid_ymin].u_type = bc_type_neumann;
        bc_vel[bcid_ymin].u_neumann = numerical_bc_value[2][0];
    }
    //设置 v 方向速度的边界条件
    if(numerical_bc_type[2][1] == "constant")
    {
        bc_vel[bcid_ymin].v_type = bc_type_dirichlet;
        bc_vel[bcid_ymin].v_dirichlet = numerical_bc_value[2][1];
    }else if(numerical_bc_type[2][1] == "flux")
    {
        bc_vel[bcid_ymin].v_type = bc_type_neumann;
        bc_vel[bcid_ymin].v_neumann = numerical_bc_value[2][1];
    }
    //设置 w 方向速度的边界条件
    if(dim == 3){
        if(numerical_bc_type[2][2] == "constant")
        {
            bc_vel[bcid_ymin].w_type = bc_type_dirichlet;
            bc_vel[bcid_ymin].w_dirichlet = numerical_bc_value[2][2];
        }else if(numerical_bc_type[2][2] == "flux")
        {
            bc_vel[bcid_ymin].w_type = bc_type_neumann;
            bc_vel[bcid_ymin].w_neumann = numerical_bc_value[2][2];
        } 
    }

    //Ymax边界条件
    //设置 u 方向速度的边界条件
    if(numerical_bc_type[3][0] == "constant")
    {
        bc_vel[bcid_ymax].u_type = bc_type_dirichlet;
        bc_vel[bcid_ymax].u_dirichlet = numerical_bc_value[3][0];
    }else if(numerical_bc_type[3][0] == "flux")
    {
        bc_vel[bcid_ymax].u_type = bc_type_neumann;
        bc_vel[bcid_ymax].u_neumann = numerical_bc_value[3][0];
    }
    //设置 v 方向速度的边界条件
    if(numerical_bc_type[3][1] == "constant")
    {
        bc_vel[bcid_ymax].v_type = bc_type_dirichlet;
        bc_vel[bcid_ymax].v_dirichlet = numerical_bc_value[3][1];
    }else if(numerical_bc_type[3][1] == "flux")
    {
        bc_vel[bcid_ymax].v_type = bc_type_neumann;
        bc_vel[bcid_ymax].v_neumann = numerical_bc_value[3][1];
    }
    //设置 w 方向速度的边界条件
    if(dim == 3){
        if(numerical_bc_type[3][2] == "constant")
        {
            bc_vel[bcid_ymax].w_type = bc_type_dirichlet;
            bc_vel[bcid_ymax].w_dirichlet = numerical_bc_value[3][2];
        }else if(numerical_bc_type[3][2] == "flux")
        {
            bc_vel[bcid_ymax].w_type = bc_type_neumann;
            bc_vel[bcid_ymax].w_neumann = numerical_bc_value[3][2];
        } 
    }

    if(dim == 3){
        //Zmin边界条件
        //设置 u 方向速度的边界条件
        if(numerical_bc_type[4][0] == "constant")
        {
            bc_vel[bcid_zmin].u_type = bc_type_dirichlet;
            bc_vel[bcid_zmin].u_dirichlet = numerical_bc_value[4][0];
        }else if(numerical_bc_type[4][0] == "flux")
        {
            bc_vel[bcid_zmin].u_type = bc_type_neumann;
            bc_vel[bcid_zmin].u_neumann = numerical_bc_value[4][0];
        }
        //设置 v 方向速度的边界条件
        if(numerical_bc_type[4][1] == "constant")
        {
            bc_vel[bcid_zmin].v_type = bc_type_dirichlet;
            bc_vel[bcid_zmin].v_dirichlet = numerical_bc_value[4][1];
        }else if(numerical_bc_type[4][1] == "flux")
        {
            bc_vel[bcid_zmin].v_type = bc_type_neumann;
            bc_vel[bcid_zmin].v_neumann = numerical_bc_value[4][1];
        }
        //设置 w 方向速度的边界条件
        if(dim == 3){
            if(numerical_bc_type[4][2] == "constant")
            {
                bc_vel[bcid_zmin].w_type = bc_type_dirichlet;
                bc_vel[bcid_zmin].w_dirichlet = numerical_bc_value[4][2];
            }else if(numerical_bc_type[4][2] == "flux")
            {
                bc_vel[bcid_zmin].w_type = bc_type_neumann;
                bc_vel[bcid_zmin].w_neumann = numerical_bc_value[4][2];
            } 
        }

        //Zmax边界条件
        //设置 u 方向速度的边界条件
        if(numerical_bc_type[5][0] == "constant")
        {
            bc_vel[bcid_zmax].u_type = bc_type_dirichlet;
            bc_vel[bcid_zmax].u_dirichlet = numerical_bc_value[5][0];
        }else if(numerical_bc_type[5][0] == "flux")
        {
            bc_vel[bcid_zmax].u_type = bc_type_neumann;
            bc_vel[bcid_zmax].u_neumann = numerical_bc_value[5][0];
        }
        //设置 v 方向速度的边界条件
        if(numerical_bc_type[5][1] == "constant")
        {
            bc_vel[bcid_zmax].v_type = bc_type_dirichlet;
            bc_vel[bcid_zmax].v_dirichlet = numerical_bc_value[5][1];
        }else if(numerical_bc_type[5][1] == "flux")
        {
            bc_vel[bcid_zmax].v_type = bc_type_neumann;
            bc_vel[bcid_zmax].v_neumann = numerical_bc_value[5][1];
        }
        //设置 w 方向速度的边界条件
        if(dim == 3){
            if(numerical_bc_type[5][2] == "constant")
            {
                bc_vel[bcid_zmax].w_type = bc_type_dirichlet;
                bc_vel[bcid_zmax].w_dirichlet = numerical_bc_value[5][2];
            }else if(numerical_bc_type[5][2] == "flux")
            {
                bc_vel[bcid_zmax].w_type = bc_type_neumann;
                bc_vel[bcid_zmax].w_neumann = numerical_bc_value[5][2];
            } 
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

void BoundarySettings::DisplayNumericalBcType()
{
    std::cout << "---------- Numerical BC Type ----------" << std::endl;
    
    // 方向名称
    std::vector<std::string> dir_names = {"u", "v", "w"};
    
    // 定义边界ID和名称的对应关系（与bcid定义一致）
    std::vector<std::pair<int, std::string>> bc_list = {
        {bcid_xmin, "Xmin"},
        {bcid_xmax, "Xmax"},
        {bcid_ymin, "Ymin"},
        {bcid_ymax, "Ymax"},
        {bcid_zmin, "Zmin"},
        {bcid_zmax, "Zmax"}
    };
    
    for(const auto& bc_pair : bc_list)
    {
        int bc_id = bc_pair.first;
        const std::string& bc_name = bc_pair.second;
        
        // 跳过超出bc_vel范围的边界（2D情况下跳过Zmin, Zmax）
        if(bc_id >= static_cast<int>(bc_vel.size()))
            continue;
        
        std::cout << bc_name << ": ";
        // u方向
        if(bc_vel[bc_id].u_type == bc_type_dirichlet)
            std::cout << dir_names[0] << "=" << bc_vel[bc_id].u_dirichlet << "(constant) ";
        else if(bc_vel[bc_id].u_type == bc_type_neumann)
            std::cout << dir_names[0] << "=" << bc_vel[bc_id].u_neumann << "(flux) ";
        else
            std::cout << dir_names[0] << "=(none) ";
        
        // v方向
        if(bc_vel[bc_id].v_type == bc_type_dirichlet)
            std::cout << dir_names[1] << "=" << bc_vel[bc_id].v_dirichlet << "(constant) ";
        else if(bc_vel[bc_id].v_type == bc_type_neumann)
            std::cout << dir_names[1] << "=" << bc_vel[bc_id].v_neumann << "(flux) ";
        else
            std::cout << dir_names[1] << "=(none) ";
        
        // w方向（如果是3D）
        if(bc_vel[bc_id].w_type == bc_type_dirichlet)
            std::cout << dir_names[2] << "=" << bc_vel[bc_id].w_dirichlet << "(constant)";
        else if(bc_vel[bc_id].w_type == bc_type_neumann)
            std::cout << dir_names[2] << "=" << bc_vel[bc_id].w_neumann << "(flux)";
        else
            std::cout << dir_names[2] << "=(none)";
        
        std::cout << std::endl;
    }
}

void BoundarySettings::DisplayNumericalBcValue()
{
    std::cout << "---------- Numerical BC Value ----------" << std::endl;
    std::cout << "bc_vel size: " << bc_vel.size() << std::endl;
    
    // 设置浮点数输出精度
    std::cout << std::fixed << std::setprecision(6);
    
    // 方向名称
    std::vector<std::string> dir_names = {"u", "v", "w"};
    
    // 定义边界ID和名称的对应关系（与bcid定义一致）
    std::vector<std::pair<int, std::string>> bc_list = {
        {bcid_xmin, "Xmin"},
        {bcid_xmax, "Xmax"},
        {bcid_ymin, "Ymin"},
        {bcid_ymax, "Ymax"},
        {bcid_zmin, "Zmin"},
        {bcid_zmax, "Zmax"}
    };
    
    for(const auto& bc_pair : bc_list)
    {
        int bc_id = bc_pair.first;
        const std::string& bc_name = bc_pair.second;
        
        // 跳过超出bc_vel范围的边界（2D情况下跳过Zmin, Zmax）
        if(bc_id >= static_cast<int>(bc_vel.size())) {
            std::cout << bc_name << ": SKIPPED (id=" << bc_id << " >= size=" << bc_vel.size() << ")" << std::endl;
            continue;
        }
        
        std::cout << bc_name << " (array index=" << bc_id << "):" << std::endl;
        
        // u方向数值
        std::cout << "  " << dir_names[0] << ": ";
        std::cout << "Dirichlet=" << bc_vel[bc_id].u_dirichlet << ", ";
        std::cout << "Neumann=" << bc_vel[bc_id].u_neumann << std::endl;
        
        // v方向数值
        std::cout << "  " << dir_names[1] << ": ";
        std::cout << "Dirichlet=" << bc_vel[bc_id].v_dirichlet << ", ";
        std::cout << "Neumann=" << bc_vel[bc_id].v_neumann << std::endl;
        
        // w方向数值
        std::cout << "  " << dir_names[2] << ": ";
        std::cout << "Dirichlet=" << bc_vel[bc_id].w_dirichlet << ", ";
        std::cout << "Neumann=" << bc_vel[bc_id].w_neumann << std::endl;
    }
    
    // 恢复默认输出格式
    std::cout << std::defaultfloat;
}
