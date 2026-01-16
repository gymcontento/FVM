/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-03 00:35:43
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-17 00:41:19
 * FilePath: \convection_diffusion\src\main.cpp
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
#include "solver.h"
#include <ctime>
#include <vector>
#include <iomanip>
#include <map>
#include <fstream>
#include <chrono> 


int main(){
    clock_t start = clock();

    //设置网格参数
    int dim = 2;
    std::vector<int> cellnums{64,1,1};
    std::vector<std::vector<float>> range = {
        {0.0f, 0.833f}, {0.0f, 1.0f}, {0.0f, 1.0f}};
    range[0][0]  = -0.5f/float(cellnums[0]-1);
    range[0][1] = 1.0f + 0.5f/float(cellnums[0]-1);

    float heatsource = 0.0f;
    float initialT = 0.5f;
    float initialUF = 1.0f; //uf
    float initialVF = 0.0f;
    float initialWF = 0.0f;
    
    StructureMesh case1;
    case1.CreateMesh(dim, cellnums);
    case1.CreateCoordinates(range);
    case1.CreateFieldMeshData();
    case1.SetHeatSource(heatsource);
    case1.SetInitialT(initialT);
    case1.SetInitialUF(initialUF);
    case1.SetInitialVF(initialVF);
    case1.SetInitialWF(initialWF);
    case1.CreateCoeffMeshData();

    //求解参数设置
    std::string eqn_type = "conduction";
    std::string conv_scheme = "cd"; //默认upwind "uw" "cd" "CD" "pl" "PL" "sou" "SOU"
    int nlin_steps = 200;          //非线性求解次数
    float nlinsol_residual = 1.e-8f;        // 非线性迭代的收敛误差判定值
    int t_iter_nums = 10000;           // 温度线性求解器迭代次数
    float relax_coef = 0.75f;          // 松弛因子
    float linsol_residual_t = 0.01f;             // 线性求解器迭代的收敛误差判定值
    SolverSettings solversettings;
    solversettings.SetEqnType(eqn_type);
    solversettings.SetConvScheme(conv_scheme);
    solversettings.SetIterNums(nlin_steps);
    solversettings.SetThermalSolverParam(t_iter_nums, relax_coef,
         linsol_residual_t, nlinsol_residual);
    solversettings.CheckSolverSettings();
    std::cout << "------------" << "\n";
    std::cout << "Conv Scheme: " << solversettings.conv_scheme << "\n";
    
    //材料参数设置
    float density{1.0f};
    float viscosity{1.0f};
    float conductivity{1000000.0f};
    float specficheat{0.0f};
    MaterialSettings material;
    material.SetInitialDenstVist(case1,density,viscosity);
    material.SetCondtSpcHt(conductivity, specficheat);
    material.CheckDensityViscosity(case1);

    //边界条件设置
    //fluid: inlet, wall, outlet
    //xmin xmax ymin ymax zmin zmax
    std::vector<std::string> physical_bc_type{"wall","wall","wall","wall","wall","inlet"};
    //temp: constant heat_flux
    std::vector<std::string> numerical_bc_type{"constant","constant","constant","constant","constant","constant"};
    std::vector<float> numerical_bc_value{373.0f, 373.0f, 373.0f, 293.0f, 0.0f, 0.0f};
    //benchmark所用边界条件
    if(cellnums[1] == 1){
        physical_bc_type[0] = "wall"; numerical_bc_type[0] = "constant"; numerical_bc_value[0] = 0.0f;
        physical_bc_type[1] = "wall"; numerical_bc_type[1] = "constant"; numerical_bc_value[1] = 1.0f;
        physical_bc_type[2] = "wall"; numerical_bc_type[2] = "heat_flux"; numerical_bc_value[2] = 0.0f;
        physical_bc_type[3] = "wall"; numerical_bc_type[3] = "heat_flux"; numerical_bc_value[3] = 0.0f;
        physical_bc_type[4] = "wall"; numerical_bc_type[4] = "heat_flux"; numerical_bc_value[4] = 0.0f;
        physical_bc_type[5] = "wall"; numerical_bc_type[5] = "heat_flux"; numerical_bc_value[5] = 0.0f;
    }else if(cellnums[0] == 1){
        physical_bc_type[0] = "wall"; numerical_bc_type[0] = "heat_flux"; numerical_bc_value[0] = 0.0f;
        physical_bc_type[1] = "wall"; numerical_bc_type[1] = "heat_flux"; numerical_bc_value[1] = 0.0f;
        physical_bc_type[2] = "wall"; numerical_bc_type[2] = "constant"; numerical_bc_value[2] = 0.0f;
        physical_bc_type[3] = "wall"; numerical_bc_type[3] = "constant"; numerical_bc_value[3] = 1.0f;
        physical_bc_type[4] = "wall"; numerical_bc_type[4] = "heat_flux"; numerical_bc_value[4] = 0.0f;
        physical_bc_type[5] = "wall"; numerical_bc_type[5] = "heat_flux"; numerical_bc_value[5] = 0.0f;        
    }else if(cellnums[2] == 1){
        if(dim == 3){
            physical_bc_type[0] = "wall"; numerical_bc_type[0] = "heat_flux"; numerical_bc_value[0] = 0.0f;
            physical_bc_type[1] = "wall"; numerical_bc_type[1] = "heat_flux"; numerical_bc_value[1] = 0.0f;
            physical_bc_type[2] = "wall"; numerical_bc_type[2] = "heat_flux"; numerical_bc_value[2] = 0.0f;
            physical_bc_type[3] = "wall"; numerical_bc_type[3] = "heat_flux"; numerical_bc_value[3] = 0.0f;
            physical_bc_type[4] = "wall"; numerical_bc_type[4] = "constant"; numerical_bc_value[4] = 0.0f;
            physical_bc_type[5] = "wall"; numerical_bc_type[5] = "constant"; numerical_bc_value[5] = 1.0f;
        }
    }

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

    //后处理 残差输出频率，vtk文件输出频率
    int res_freq = 2;
    int out_freq = 1000;
    PostProcess post;
    post.SetOuputFreq(res_freq, out_freq);
    post.WriteVTKCollocated_temp(case1);
    //输出文件
    std::string& nonlinsol_fname = post.nonlinsol_fname;
    std::string& linsol_fname = post.linsol_fname;

    // Open nonlinear residual file (写入模式，会覆盖已有文件)
    std::ofstream nonlinsol_fid(nonlinsol_fname);
    nonlinsol_fid << "#it, walltime, l2_t/l2_max_t\n";
    nonlinsol_fid.close();

    // Open linear residual file (写入模式，会覆盖已有文件)
    std::ofstream linsol_fid(linsol_fname);
    linsol_fid << "#it_nl, it, tot_it, norm, init, max, rel\n";
    linsol_fid.close();


    //开始仿真
    AssembleSystem linearsystem;
    Solver solver;

    //求解 
    //外层循环为时间循环，每次都会重新组装矩阵，因为系数和源项有可能随着场量的变化而变化（因此是非线性）
    //内存循环组装线性系统，并进行求解
    for(int i=1; i < nlin_steps+1; ++i)
    {
        if(i % 2 == 0 || i == 1 || i == nlin_steps)
        {
            std::cout << "\n";
            std::cout << "----------------------------" << "\n";
            std::cout << "Begin iter = " << i << "\n";
            std::cout << "----------------------------" << "\n";
        }

        if(solversettings.eqn_type == solversettings.eqn_type_conduct){
            case1.t0 = case1.t; //t0为了计算范数做准备
            //计算温度方程coefs
            linearsystem.ConductionCoefs(case1,solversettings,material);
            //对网格的边界类型做处理，修改coefs
            linearsystem.ConductionCoefsBoud(case1, material, boundary);
            
            //进行雅克比迭代
            // solver.JacobiIteraMethod(i, case1, post, solversettings);
            solver.GaussSiedelMethod(i, case1, post, solversettings);

            //将原来的值更新，这个结果与非线性迭代的收敛误差判定值对比
            std::vector<float> l2_para = solver.CalScalarNorm2(case1, solversettings, "temp");
            solversettings.l2_t = l2_para[0];
            solversettings.l2_max_t = l2_para[1];
    
            //如果满足条件就令stop_sim为真，程序将在主文件那里跳出迭代循环
            if(solversettings.l2_t / solversettings.l2_max_t < solversettings.nlinsol_residual){
                solver.stop_sim = true;
                std::cout << std::endl;
                std::cout << "----------------------------" << std::endl;
                std::cout << "Final iter = " << i << std::endl;
                std::cout << "it, l2_t/l2_max_t " << i << " " << solversettings.l2_t / solversettings.l2_max_t << std::endl;
                std::cout << "----------------------------" << std::endl;
            }
        }
        
        if(i % post.OutfreqExport() == 0 || i == nlin_steps || solver.stop_sim){
            post.WriteVTKCollocated_temp(case1);
        }

        if (i == 1 || i % post.ResfreqExport() == 0 || i == nlin_steps || solver.stop_sim) {
            // 计算运行时间（从clock_begin到现在）
            clock_t current = clock();
            double duration_time = double(current - start) / CLOCKS_PER_SEC;
            
            // 输出到控制台
            std::cout << "it, walltime, l2_t/l2_max_t " << i << " " << duration_time << " " << solversettings.l2_t / solversettings.l2_max_t << std::endl;
            
            // 以追加模式打开文件并写入
            std::ofstream nonlinsol_fid(post.nonlinsol_fname, std::ios::app);
            if (nonlinsol_fid.is_open()) {
                nonlinsol_fid << i << " " << duration_time << " " << solversettings.l2_t / solversettings.l2_max_t << "\n";
                nonlinsol_fid.close();
            } else {
                std::cerr << "Error: Unable to open file " << post.nonlinsol_fname << std::endl;
            }
        }

        if(solver.stop_sim == true){
            break;
        }
    }

    //针对一维特殊情况做的输出
    // post.WriteVTKCollocated_temp_Pe_L_center(case1, material, solversettings);
    post.WriteVTKCollocated_temp_Pe_L(case1, material);
    
    clock_t end = clock();
    double duration = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << 
        std::fixed << std::setprecision(10) 
        << duration << " seconds" << std::endl;
    
    return 0;
}
