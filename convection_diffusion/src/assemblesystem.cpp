/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-07 23:45:38
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-15 20:52:21
 * FilePath: \convection_diffusion\src\assemblesystem.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "assemblesystem.h"
#include "kernel.h"
#include <fstream>
#include <iomanip>

void AssembleSystem::ConductionCoefs(StructureMesh &mesh, SolverSettings &solversettings, MaterialSettings &material)
{
    Kernel inport_kernel;
    //获取面积和体积
    const int& dim = mesh.dim_export();
    const auto& dx = mesh.dx_export();
    const auto& dy = mesh.dy_export();
    const auto& dz = mesh.dz_export();
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
    float idt = 1.0f / solversettings.dt;
    int conv_scheme = solversettings.conv_scheme;

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
                aE = inport_kernel.CalCoefA(conv_scheme, area_vol[0], dx, ul, ur, right_cond, left_cond, rho, -1.0f); 

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
                aW = inport_kernel.CalCoefA(conv_scheme, area_vol[0], dx, ul, ur, right_cond, left_cond, rho, 1.0f); 
                
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
                aN = inport_kernel.CalCoefA(conv_scheme, area_vol[1], dy, vl, vr, right_cond, left_cond, rho, -1.0f);

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
                aS = inport_kernel.CalCoefA(conv_scheme, area_vol[1], dy, vl, vr, right_cond, left_cond, rho, 1.0f);

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
                    aT = inport_kernel.CalCoefA(conv_scheme, area_vol[2], dz, wl, wr, right_cond, left_cond, rho, -1.0f); 

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
                    aB = inport_kernel.CalCoefA(conv_scheme, area_vol[2], dz, wl, wr, right_cond, left_cond, rho, 1.0f); 
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
            }
        }
    }

    return;
}

void AssembleSystem::ConductionCoefsBoud(StructureMesh &mesh, MaterialSettings &material, BoundarySettings &boundary)
{
    //维度
    const int& dim = mesh.dim_export();
    
    // linearsystem系统中的a_coef
    int& id_aP = mesh.id_aP;
    int& id_aE = mesh.id_aE;
    int& id_aW = mesh.id_aW;
    int& id_aN = mesh.id_aN;
    int& id_aS = mesh.id_aS;
    int& id_aT = mesh.id_aT;
    int& id_aB = mesh.id_aB;
    int& id_bsrc = mesh.id_bsrc;

    //网格点的坐标 和边界坐标
    const auto& x_cells = mesh.xcells_export();
    const auto& y_cells = mesh.ycells_export();
    const auto& z_cells = mesh.zcells_export();
    const auto& x_nodes = mesh.xnodes_export();
    const auto& y_nodes = mesh.ynodes_export();
    const auto& z_nodes = mesh.znodes_export();

    //单元的面编号
    int& face_id_e = boundary.fid_e;
    int& face_id_w = boundary.fid_w;
    int& face_id_n = boundary.fid_n;
    int& face_id_s = boundary.fid_s;
    int& face_id_t = boundary.fid_t;
    int& face_id_b = boundary.fid_b;

    //边界类型(物理边界)
    int bc_none   = 0;
    int bc_wall   = 1;
    int bc_inlet  = 2;
    int bc_outlet = 3;
    const auto& bc_physics_type = boundary.BcPhysicsExport();
    int bc_physics_type_e, bc_physics_type_w, bc_physics_type_n, bc_physics_type_s, bc_physics_type_t, bc_physics_type_b;
    //温度边界类型
    const auto& conductivity = material.ConductivityExport();
    const auto& specfiheat = material.SpecficheatExport();
    const auto& bc_temp = boundary.BcTempExport();
    int t_bc_type_dirichlet = 0;
    int t_bc_type_neumann = 1;
    // float bc_dirichlet_value_e, bc_dirichlet_value_w, bc_dirichlet_value_n, bc_dirichlet_value_s, bc_dirichlet_value_t, bc_dirichlet_value_b;

    
    //定义与fluidboundary对象相关的局部变量
    const auto& cellnum = mesh.cellnum_export();
    const auto& bcid = boundary.BcidExport();
    int bcid_e, bcid_w, bcid_n, bcid_s, bcid_t, bcid_b;

    for(int k=0; k < cellnum[2]; ++k)
    {
        for(int j=0; j < cellnum[1]; ++j)
        {
            for(int i=0; i < cellnum[0]; ++i)
            {
                //确定单元的界面类型，是否是边界，如果是，是什么类型
                bcid_e = bcid[i][j][k][face_id_e];
                bcid_w = bcid[i][j][k][face_id_w];
                bcid_n = bcid[i][j][k][face_id_n];
                bcid_s = bcid[i][j][k][face_id_s];
                if(dim == 3)
                {
                    bcid_t = bcid[i][j][k][face_id_t];
                    bcid_b = bcid[i][j][k][face_id_b];
                }

                //设定界面类型，边界和不是边界
                bc_physics_type_e = bc_physics_type[bcid_e];
                bc_physics_type_w = bc_physics_type[bcid_w];
                bc_physics_type_n = bc_physics_type[bcid_n];
                bc_physics_type_s = bc_physics_type[bcid_s];
                if(dim == 3)
                {
                    bc_physics_type_t = bc_physics_type[bcid_t];
                    bc_physics_type_b = bc_physics_type[bcid_b];
                }

                //dirichlet和neumann边界条件的值
                //东方向
                if(bc_physics_type_e == bc_inlet){
                    //计算的系数a都是为正的，cellcoeff中的除aP的系数都是已经加上负号，因此这里要再次加上负号变成正的
                    mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aE];
                    mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aE] * bc_temp[bcid_e].t_dirichlet;
                    mesh.cellcoeff[i][j][k][id_aE] = 0.0f;
                } else if(bc_physics_type_e == bc_outlet){
                    mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aE];
                    mesh.cellcoeff[i][j][k][id_aE] = 0.0f;
                }else if(bc_physics_type_e == bc_wall){
                    if(bc_temp[bcid_e].bc_t_type == t_bc_type_dirichlet){
                        mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aE];
                        mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aE] * bc_temp[bcid_e].t_dirichlet;
                        mesh.cellcoeff[i][j][k][id_aE] = 0.0f;
                    }else if(bc_temp[bcid_e].bc_t_type == t_bc_type_neumann){
                        mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aE];
                        //通量的负号和cellcoef的负号，负负得正
                        //第二类边界条件一般是给出通量q=-k\frac{\partial u}{\partial x}，因此这里计算的时候也需要消去k 即除以con
                        mesh.cellcoeff[i][j][k][id_bsrc] += mesh.cellcoeff[i][j][k][id_aE] * 
                                                        bc_temp[bcid_e].t_neumann * 2.0f * (x_nodes[i+1] - x_cells[i]) / conductivity;
                        mesh.cellcoeff[i][j][k][id_aE] = 0.0f;
                    }
                }
            
                //西方向
                if(bc_physics_type_w == bc_inlet){
                    //计算的系数a都是为正的，cellcoeff中的除aP的系数都是已经加上负号，因此这里要再次加上负号变成正的
                    mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aW];
                    mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aW] * bc_temp[bcid_w].t_dirichlet;
                    mesh.cellcoeff[i][j][k][id_aW] = 0.0f;
                } else if(bc_physics_type_w == bc_outlet){
                    mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aW];
                    mesh.cellcoeff[i][j][k][id_aW] = 0.0f;
                }else if(bc_physics_type_w == bc_wall){
                    if(bc_temp[bcid_w].bc_t_type == t_bc_type_dirichlet){
                        mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aW];
                        mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aW] * bc_temp[bcid_w].t_dirichlet;
                        mesh.cellcoeff[i][j][k][id_aW] = 0.0f;
                    }else if(bc_temp[bcid_w].bc_t_type == t_bc_type_neumann){
                        mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aW];
                        //通量的负号和cellcoef的负号，负负得正
                        //第二类边界条件一般是给出通量q=-k\frac{\partial u}{\partial x}，因此这里计算的时候也需要消去k 即除以con
                        mesh.cellcoeff[i][j][k][id_bsrc] += mesh.cellcoeff[i][j][k][id_aW] * 
                                                        bc_temp[bcid_w].t_neumann * 2.0f * (x_nodes[i] - x_cells[i]) / conductivity;
                        mesh.cellcoeff[i][j][k][id_aW] = 0.0f;
                    }
                }

                //北方向
                if(bc_physics_type_n == bc_inlet){
                    //计算的系数a都是为正的，cellcoeff中的除aP的系数都是已经加上负号，因此这里要再次加上负号变成正的
                    mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aN];
                    mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aN] * bc_temp[bcid_n].t_dirichlet;
                    mesh.cellcoeff[i][j][k][id_aN] = 0.0f;
                } else if(bc_physics_type_n == bc_outlet){
                    mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aN];
                    mesh.cellcoeff[i][j][k][id_aN] = 0.0f;
                }else if(bc_physics_type_n == bc_wall){
                    if(bc_temp[bcid_n].bc_t_type == t_bc_type_dirichlet){
                        mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aN];
                        mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aN] * bc_temp[bcid_n].t_dirichlet;
                        mesh.cellcoeff[i][j][k][id_aN] = 0.0f;
                    }else if(bc_temp[bcid_n].bc_t_type == t_bc_type_neumann){
                        mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aN];
                        //通量的负号和cellcoef的负号，负负得正
                        //第二类边界条件一般是给出通量q=-k\frac{\partial u}{\partial x}，因此这里计算的时候也需要消去k 即除以con
                        mesh.cellcoeff[i][j][k][id_bsrc] += mesh.cellcoeff[i][j][k][id_aN] * 
                                                        bc_temp[bcid_n].t_neumann * 2.0f * (y_nodes[j+1] - y_cells[j]) / conductivity;
                        mesh.cellcoeff[i][j][k][id_aN] = 0.0f;
                    }
                }

                //南方向
                if(bc_physics_type_s == bc_inlet){
                    //计算的系数a都是为正的，cellcoeff中的除aP的系数都是已经加上负号，因此这里要再次加上负号变成正的
                    mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aS];
                    mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aS] * bc_temp[bcid_s].t_dirichlet;
                    mesh.cellcoeff[i][j][k][id_aS] = 0.0f;
                } else if(bc_physics_type_s == bc_outlet){
                    mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aS];
                    mesh.cellcoeff[i][j][k][id_aS] = 0.0f;
                }else if(bc_physics_type_s == bc_wall){
                    if(bc_temp[bcid_s].bc_t_type == t_bc_type_dirichlet){
                        mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aS];
                        mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aS] * bc_temp[bcid_s].t_dirichlet;
                        mesh.cellcoeff[i][j][k][id_aS] = 0.0f;
                    }else if(bc_temp[bcid_s].bc_t_type == t_bc_type_neumann){
                        mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aS];
                        //通量的负号和cellcoef的负号，负负得正
                        //第二类边界条件一般是给出通量q=-k\frac{\partial u}{\partial x}，因此这里计算的时候也需要消去k 即除以con
                        mesh.cellcoeff[i][j][k][id_bsrc] += mesh.cellcoeff[i][j][k][id_aS] * 
                                                        bc_temp[bcid_s].t_neumann * 2.0f * (y_nodes[j] - y_cells[j]) / conductivity;
                        mesh.cellcoeff[i][j][k][id_aS] = 0.0f;
                    }
                }
                
                if(dim == 3)
                {
                    //顶面
                    if(bc_physics_type_t == bc_inlet){
                        //计算的系数a都是为正的，cellcoeff中的除aP的系数都是已经加上负号，因此这里要再次加上负号变成正的
                        mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aT];
                        mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aT] * bc_temp[bcid_t].t_dirichlet;
                        mesh.cellcoeff[i][j][k][id_aT] = 0.0f;
                    } else if(bc_physics_type_t == bc_outlet){
                        mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aT];
                        mesh.cellcoeff[i][j][k][id_aT] = 0.0f;
                    }else if(bc_physics_type_t == bc_wall){
                        if(bc_temp[bcid_t].bc_t_type == t_bc_type_dirichlet){
                            mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aT];
                            mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aT] * bc_temp[bcid_t].t_dirichlet;
                            mesh.cellcoeff[i][j][k][id_aT] = 0.0f;
                        }else if(bc_temp[bcid_t].bc_t_type == t_bc_type_neumann){
                            mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aT];
                            //通量的负号和cellcoef的负号，负负得正
                            //第二类边界条件一般是给出通量q=-k\frac{\partial u}{\partial x}，因此这里计算的时候也需要消去k 即除以con
                            mesh.cellcoeff[i][j][k][id_bsrc] += mesh.cellcoeff[i][j][k][id_aT] * 
                                                            bc_temp[bcid_t].t_neumann * 2.0f * (z_nodes[k+1] - z_cells[k]) / conductivity;
                            mesh.cellcoeff[i][j][k][id_aT] = 0.0f;
                        }
                    }

                    //底面
                    if(bc_physics_type_b == bc_inlet){
                        //计算的系数a都是为正的，cellcoeff中的除aP的系数都是已经加上负号，因此这里要再次加上负号变成正的
                        mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aB];
                        mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aB] * bc_temp[bcid_b].t_dirichlet;
                        mesh.cellcoeff[i][j][k][id_aB] = 0.0f;
                    } else if(bc_physics_type_b == bc_outlet){
                        mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aB];
                        mesh.cellcoeff[i][j][k][id_aB] = 0.0f;
                    }else if(bc_physics_type_b == bc_wall){
                        if(bc_temp[bcid_b].bc_t_type == t_bc_type_dirichlet){
                            mesh.cellcoeff[i][j][k][id_aP] -= mesh.cellcoeff[i][j][k][id_aB];
                            mesh.cellcoeff[i][j][k][id_bsrc] -= 2.0f * mesh.cellcoeff[i][j][k][id_aB] * bc_temp[bcid_b].t_dirichlet;
                            mesh.cellcoeff[i][j][k][id_aB] = 0.0f;
                        }else if(bc_temp[bcid_b].bc_t_type == t_bc_type_neumann){
                            mesh.cellcoeff[i][j][k][id_aP] += mesh.cellcoeff[i][j][k][id_aB];
                            //通量的负号和cellcoef的负号，负负得正
                            //第二类边界条件一般是给出通量q=-k\frac{\partial u}{\partial x}，因此这里计算的时候也需要消去k 即除以con
                            mesh.cellcoeff[i][j][k][id_bsrc] += mesh.cellcoeff[i][j][k][id_aB] * 
                                                            bc_temp[bcid_b].t_neumann * 2.0f * (z_nodes[k] - z_cells[k]) / conductivity;
                            mesh.cellcoeff[i][j][k][id_aB] = 0.0f;
                        }
                    }
                }
            }
        }
    }

    // Forced cell-centered boundary conditions for 红宝书的例子，例子见书上Fig. 11.7
    for(int k=0; k < cellnum[2]; ++k)
    {
        for(int j=0; j < cellnum[1]; ++j)
        {
            for(int i=0; i < cellnum[0]; ++i)
            {
                if(i == 0){
                    mesh.cellcoeff[i][j][k][id_aP] = 1.0f;
                    mesh.cellcoeff[i][j][k][id_aE] = 0.0f;
                    mesh.cellcoeff[i][j][k][id_aW] = 0.0f;
                    mesh.cellcoeff[i][j][k][id_aN] = 0.0f;
                    mesh.cellcoeff[i][j][k][id_aS] = 0.0f;
                    if(dim == 3){
                        mesh.cellcoeff[i][j][k][id_aT] = 0.0f;
                        mesh.cellcoeff[i][j][k][id_aB] = 0.0f;
                    }
                    mesh.cellcoeff[i][j][k][id_bsrc] = bc_temp[bcid_w].t_dirichlet;
                }else if(i == cellnum[0]-1){
                    mesh.cellcoeff[i][j][k][id_aP] = 1.0f;
                    mesh.cellcoeff[i][j][k][id_aE] = 0.0f;
                    mesh.cellcoeff[i][j][k][id_aW] = 0.0f;
                    mesh.cellcoeff[i][j][k][id_aN] = 0.0f;
                    mesh.cellcoeff[i][j][k][id_aS] = 0.0f;
                    if(dim == 3){
                        mesh.cellcoeff[i][j][k][id_aT] = 0.0f;
                        mesh.cellcoeff[i][j][k][id_aB] = 0.0f;
                    }
                    mesh.cellcoeff[i][j][k][id_bsrc] = bc_temp[bcid_e].t_dirichlet;
                }
                
            }
        }
    }

    // std::cout << "---------- Check cellcoef ----------" << std::endl;
    // std::cout << "cellcoeff[0][0][0][id_aN]:" << mesh.cellcoeff[0][0][0][id_aN] << "\n";
    // std::cout << "cellcoeff[0][0][0][id_aP]:" << mesh.cellcoeff[0][0][0][id_aP] << "\n";
    // std::cout << "cellcoeff[0][0][0][id_bsrc]:" << mesh.cellcoeff[0][0][0][id_bsrc] << "\n";
    // std::cout << "cellcoeff[1][0][0][id_aN]:" << mesh.cellcoeff[1][0][0][id_aN] << "\n";
    // std::cout << "cellcoeff[1][0][0][id_aP]:" << mesh.cellcoeff[1][0][0][id_aP] << "\n";
    // std::cout << "cellcoeff[1][0][0][id_bsrc]:" << mesh.cellcoeff[1][0][0][id_bsrc] << "\n";
    // std::cout << "cellcoeff[2][0][0][id_aN]:" << mesh.cellcoeff[2][0][0][id_aN] << "\n";
    // std::cout << "cellcoeff[2][0][0][id_aP]:" << mesh.cellcoeff[2][0][0][id_aP] << "\n";
    // std::cout << "cellcoeff[2][0][0][id_bsrc]:" << mesh.cellcoeff[2][0][0][id_bsrc] << "\n";
    // std::cout << "cellcoeff[3][0][0][id_aN]:" << mesh.cellcoeff[3][0][0][id_aN] << "\n";
    // std::cout << "cellcoeff[3][0][0][id_aP]:" << mesh.cellcoeff[3][0][0][id_aP] << "\n";
    // std::cout << "cellcoeff[3][0][0][id_bsrc]:" << mesh.cellcoeff[3][0][0][id_bsrc] << "\n";
    // std::cout << "\n";
}
