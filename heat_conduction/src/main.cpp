/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-03 00:35:43
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-08 15:02:48
 * FilePath: \FVM\heat_conduction\src\main.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "structure_mesh.h"
#include "solver_settings.h"
#include "postprocess.h"
#include "material_settings.h"
#include "boundarysettings.h"
#include "kernel.h"
#include "assemblesystem.h"
#include <ctime>
#include <vector>
#include <iomanip>
#include <map>


int main(){
    clock_t start = clock();

    //设置网格参数
    int dim = 2;
    std::vector<int> cellnums{4,3,1};
    std::vector<std::vector<float>> range = {
        {0.0f, 0.833f}, {0.0f, 0.833f}, {0.0f, 1.0f}};
    float initalT = 373.0f;
    
    StructureMesh case1;
    case1.CreateMesh(dim, cellnums);
    case1.CreateCoordinates(range);
    case1.CreateFieldMeshData();
    case1.SetInitialT(initalT);
    case1.CreateCoeffMeshData();

    //求解参数设置
    int nlin_steps = 2000;          //非线性求解次数
    float nlinsol_residual = 1.e-6f;        // 非线性迭代的收敛误差判定值
    int t_iter_nums = 100;           // 温度线性求解器迭代次数
    float relax_coef = 0.75f;          // 松弛因子
    float linsol_residual_t = 0.1f;             // 线性求解器迭代的收敛误差判定值
    SolverSettings solver;
    solver.SetIterNums(nlin_steps);
    solver.SetThermalSolverParam(t_iter_nums, relax_coef,
         linsol_residual_t, nlinsol_residual);
    solver.CheckSolverSettings();
    
    //材料参数设置
    float density{1.0f};
    float viscosity{0.001f};
    float conductivity{81.0f};
    float specficheat{1.0f};
    MaterialSettings material;
    material.SetInitialDenstVist(case1,density,viscosity);
    material.SetCondtSpcHt(conductivity, specficheat);
    material.CheckDensityViscosity(case1);

    //边界条件设置
    //fluid: inlet, wall, outlet
    std::vector<std::string> physical_bc_type{"inlet","wall","inlet","wall","wall","inlet"};
    //temp: constant heat_flux
    std::vector<std::string> numerical_bc_type{"constant","constant","constant","constant","constant","constant"};
    std::vector<float> numerical_bc_value{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
    BoundarySettings boundary;
    boundary.SetBcFacesForCells(case1);
    boundary.SetPhysicalBoundary(case1, physical_bc_type[0], physical_bc_type[1], physical_bc_type[2],
                            physical_bc_type[3], physical_bc_type[4], physical_bc_type[5]);
    boundary.SetNumericalBoundary(case1, numerical_bc_type[0], numerical_bc_value[0],
                                numerical_bc_type[1], numerical_bc_value[1],
                                numerical_bc_type[2], numerical_bc_value[2],
                                numerical_bc_type[3], numerical_bc_value[3],
                                numerical_bc_type[4], numerical_bc_value[4],
                                numerical_bc_type[5], numerical_bc_value[5]);
    boundary.DisplayFluidType();                            

    //后处理
    PostProcess temp1;
    temp1.WriteVTKCollocated_temp(case1);

    //输出设置
    int res_freq = 2;

    //开始仿真
    // for(int i=1; i < nlin_steps+1; ++i)
    // {
        
    // }
    // {
    //     if(i % 2 == 0 || i == 1 || i == nlin_steps)
    //     {
    //         std::cout << "\n";
    //         std::cout << "----------------------------" << "\n";
    //         std::cout << "Begin iter = " << i << "\n";
    //         std::cout << "----------------------------" << "\n";
    //     }
    //     case1.
    // }

    AssembleSystem linearsystem;
    linearsystem.ConductionCoefs(case1,solver,material);

    clock_t end = clock();
    double duration = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << 
        std::fixed << std::setprecision(10) 
        << duration << " seconds" << std::endl;
}