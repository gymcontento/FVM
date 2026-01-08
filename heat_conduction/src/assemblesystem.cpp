/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-07 23:45:38
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-08 23:10:51
 * FilePath: \FVM\heat_conduction\src\assemblesystem.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "assemblesystem.h"
#include "kernel.h"

void AssembleSystem::ConductionCoefs(StructureMesh &mesh, SolverSettings &solver, MaterialSettings &material)
{
    Kernel inport_kernel;
    //获取面积和体积
    const int& dim = mesh.dim_export();
    const auto& dx = mesh.dx_export();
    const auto& dy = mesh.dy_export();
    const auto& dz = mesh.dz_export();
    std::cout << "dx: " << dx << " dy: " << dy << " dz: " << dz << "\n";
    std::vector<float> area_vol = inport_kernel.CalcuAreaVolume(dx, dy, dz);

    //材料参数
    //初始化源和源系数
    const auto& sC = material.HeatSourceExport();
    float sP = 0.0f;
    float bsrc; //C表示常数，P表示节点，b表示右端 
    //密度
    const auto& density = material.DensityExport(); 
    float rho;
    //热参数
    const auto& conductivity = material.ConductivityExport();
    const auto& specficheat = material.SpecficheatExport();
    float right_cond, left_cond;
    //速度场
    const auto& uf = mesh.uffield_export();
    const auto& vf = mesh.vffield_export();
    const auto& wf = mesh.wffield_export();
    float ul, ur, vl, vr, wl, wr;
    //温度场
    const auto& t_field = mesh.tfield_export();

    //求解参数
    float idt = 1.0f / solver.dt;

    // Coefficient storage positions in four dimensional array
    int& id_aP = mesh.id_aP;
    int& id_aE = mesh.id_aE;
    int& id_aW = mesh.id_aW;
    int& id_aN = mesh.id_aN;
    int& id_aS = mesh.id_aS;
    int& id_aT = mesh.id_aT;
    int& id_aB = mesh.id_aB;
    int& id_bsrc = mesh.id_bsrc;

    std::vector<int> cellnum = mesh.cellnum_export();

    //对每一个网格进行循环，计算系数aP, aE, aW, aN, aS, aT, aB, bsrc 
    float aE, aW, aN, aS, aT, aB, aP, aP_t;
    for(int k=0; k < cellnum[2]; ++k)
    {
        for(int j=0; j < cellnum[1]; ++j)
        {
            for(int i=0; i < cellnum[0]; ++i)
            {
                //初始化系数，令系数为零
                aE = 0.0f; aW = 0.0f;
                aN = 0.0f; aS = 0.0f;
                aT = 0.0f; aB = 0.0f;
                aP = 0.0f;

                //初始化源和源系数为零（这里的热源没有引入，所以可以不考虑）
                sP = 0.0f;
                bsrc = 0.0f;

                //计算每个网格中的系数，aP, aE, aW, aN, aS, aT, aB, bsrc
                //东方向的面(密度场使用两个单元的插值，速度使用界面速度，界面速度之前是两个单元速度的平均)
                if(i == cellnum[0]-1)
                {
                    rho = density[i][j][k];          
                }else{
                    rho = 0.5f * (density[i][j][k] + density[i+1][j][k]);
                }
                
                right_cond = conductivity;//如果有需要 也可以弄成场分布
                left_cond = conductivity;
                ul = uf[i+1][j][k];
                ur = uf[i+1][j][k];
                aE = inport_kernel.CalCoefA(area_vol[0], dx, ul, ur, right_cond, left_cond, rho, -1.0f); 

                //西方向的面(密度场插值，边界上不需要)
                if(i == 0)
                {
                    rho = density[i][j][k];          
                }else{
                    rho = 0.5f * (density[i][j][k] + density[i-1][j][k]);
                }
                
                right_cond = conductivity;//如果有需要 也可以弄成场分布
                left_cond = conductivity;
                ul = uf[i][j][k];
                ur = uf[i][j][k];
                aW = inport_kernel.CalCoefA(area_vol[0], dx, ul, ur, right_cond, left_cond, rho, 1.0f); 
                
                //北方向的面(密度场插值，边界上不需要)
                if(j == cellnum[1]-1)
                {
                    rho = density[i][j][k];          
                }else{
                    rho = 0.5f * (density[i][j][k] + density[i][j+1][k]);
                }
                
                right_cond = conductivity;//如果有需要 也可以弄成场分布
                left_cond = conductivity;
                vl = vf[i][j+1][k];
                vr = vf[i][j+1][k];
                aN = inport_kernel.CalCoefA(area_vol[1], dx, vl, vr, right_cond, left_cond, rho, -1.0f); 

                //南方向的面(密度场插值，边界上不需要)
                if(j == 0)
                {
                    rho = density[i][j][k];          
                }else{
                    rho = 0.5f * (density[i][j][k] + density[i][j-1][k]);
                }
                
                right_cond = conductivity;//如果有需要 也可以弄成场分布
                left_cond = conductivity;
                vl = vf[i][j][k];
                vr = vf[i][j][k];
                aS = inport_kernel.CalCoefA(area_vol[1], dx, vl, vr, right_cond, left_cond, rho, 1.0f); 

                if(dim == 3)
                {
                    //顶面(密度场插值，边界上不需要)
                    if(k == cellnum[2]-1)
                    {
                        rho = density[i][j][k];          
                    }else{
                        rho = 0.5f * (density[i][j][k] + density[i][j][k+1]);
                    }
                    
                    right_cond = conductivity;//如果有需要 也可以弄成场分布
                    left_cond = conductivity;
                    wl = wf[i][j][k+1];
                    wr = wf[i][j][k+1];
                    aT = inport_kernel.CalCoefA(area_vol[2], dx, wl, wr, right_cond, left_cond, rho, -1.0f); 

                    //底面(密度场插值，边界上不需要)
                    if(k == 0)
                    {
                        rho = density[i][j][k];          
                    }else{
                        rho = 0.5f * (density[i][j][k] + density[i][j][k-1]);
                    }
                    
                    right_cond = conductivity;//如果有需要 也可以弄成场分布
                    left_cond = conductivity;
                    wl = wf[i][j][k];
                    wr = wf[i][j][k];
                    aB = inport_kernel.CalCoefA(area_vol[2], dx, wl, wr, right_cond, left_cond, rho, 1.0f); 
                }
                // 对density进行更新
                rho = density[i][j][k];
                aP_t = rho * specficheat * area_vol[3] * idt;
                aP = aP + aE + aW + aN + aS + aT + aB + aP_t - sP * area_vol[3];

                //存储系数
                mesh.cellcoeff[i][j][k][id_aP] = aP;
                mesh.cellcoeff[i][j][k][id_aE] = -aE;
                mesh.cellcoeff[i][j][k][id_aW] = -aW;
                mesh.cellcoeff[i][j][k][id_aN] = -aN;
                mesh.cellcoeff[i][j][k][id_aS] = -aS;
                if(dim == 3)
                {
                    mesh.cellcoeff[i][j][k][id_aT] = -aT;
                    mesh.cellcoeff[i][j][k][id_aB] = -aB;
                }
                bsrc = sC * area_vol[3] + aP_t * t_field[i][j][k];
                mesh.cellcoeff[i][j][k][id_bsrc] = bsrc;

                std::cout << "aP_t: " << aP_t << "\t"
                        << " aP: " << aP << "\t" 
                        << " bsrc: " << bsrc << "\n";
            }
            std::cout << "\n";
        }
    }

    return;
}

void AssembleSystem::ConductionCoefsBoud(StructureMesh &mesh, MaterialSettings &material, BoundarySettings &boundary)
{
    // linearsystem系统中的a_coef
    int& id_aP = mesh.id_aP;
    int& id_aE = mesh.id_aE;
    int& id_aW = mesh.id_aW;
    int& id_aN = mesh.id_aN;
    int& id_aS = mesh.id_aS;
    int& id_aT = mesh.id_aT;
    int& id_aB = mesh.id_aB;
    int& id_bsrc = mesh.id_bsrc;

    //网格点的坐标
    const auto& x_cells = mesh.xcells_export();
    const auto& y_cells = mesh.ycells_export();
    const auto& z_cells = mesh.zcells_export();

    //边界类型(物理边界)
    const auto& bc_physics = boundary.BcPhysicsExport();

    //定义与fluidboundary对象相关的局部变量
    const auto& cellnum = mesh.cellnum_export();
    for(int k=0; k < cellnum[2]; ++k)
    {
        for(int j=0; j < cellnum[1]; ++j)
        {
            for(int i=0; i < cellnum[0]; ++i)
            {
                //确定单元的界面类型


                //东面
                 
            
            }
        }
    }
}
