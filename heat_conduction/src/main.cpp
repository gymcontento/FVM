/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-03 00:35:43
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-04 23:38:43
 * FilePath: \FVM\heat_conduction\src\main.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "structure_mesh.h"
#include "solver_settings.h"
#include "postprocess.h"
#include "material_settings.h"
#include <ctime>
#include <vector>
#include <iomanip>

int main(){
    clock_t start = clock();

    int dim = 2;
    std::vector<int> numcell{4,3,2};
    std::vector<std::vector<float>> range = {
        {0.0f, 1.0f}, {0.0f, 1.0f}, {0.0f, 1.0f}};
    float initalT = 100.23f;
    
    //网格参数设置
    StructureMesh case1;
    case1.CreateMesh(dim, numcell);
    case1.CreateCoordinates(range);
    case1.CreateFieldMeshData();
    case1.SetInitialT(initalT);
    case1.CreateCoeffMeshData();

    //求解参数设置
    int nlin_steps = 2000;          //求解系数设置
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

    // PostProcess temp1;
    // temp1.WriteVTKCollocated_temp(case1);

    clock_t end = clock();
    double duration = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << 
        std::fixed << std::setprecision(10) 
        << duration << " seconds" << std::endl;
}