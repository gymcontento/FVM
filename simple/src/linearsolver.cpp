/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-09 19:28:36
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-02-09 23:15:40
 * FilePath: \simple\src\linearsolver.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_Name_email}, All Rights Reserved. 
 ******************************/
#include "linearsolver.h"
#include "solver_settings.h"
#include <cmath>
#include <fstream>

std::vector<float> LinearSolver::CalScalarNorm2(StructureMesh &mesh, SolverSettings& solversettings, std::string var)
{
    const auto& cellnum = mesh.cellnum_export();
    
    float l2_t, l2_max_t;
    if(var == "temp"){
        l2_t = solversettings.l2_t;
        l2_max_t = solversettings.l2_max_t;
    }

    //计算范数
    float res_t;//残差res_u
    for(int k=0; k < cellnum[2]; ++k)
    {
        for(int j=0; j < cellnum[1]; ++j)
        {
            for(int i=0; i < cellnum[0]; ++i)
            {
                res_t = mesh.t[i][j][k] - mesh.t0[i][j][k];
                l2_t = res_t * res_t;
            }
        }
    }
    std::cout << "l2_t: " << l2_t << "\n";
    int totalcell = cellnum[0] * cellnum[1] * cellnum[2];
    l2_t = std::sqrt(l2_t/float(totalcell));
    std::cout << "l2_t: " << l2_t << "\n";

    l2_max_t = std::max(l2_t,l2_max_t);
    std::vector<float> array_l2{l2_t,l2_max_t};
    
    return array_l2;
}

void LinearSolver::GaussSiedelMethod(int& nlitera_nums,std::vector<std::vector<std::vector<std::vector<float>>>>& cc,
                                     StructureMesh& mesh, PostProcess& postprocess, 
                                     SolverSettings& solversettings, int& iter_num, int& total_iter_num, float& rc, float& torlence)
{
    // 定义与网格相关的变量
    const auto& dim = mesh.dim_export();
    const auto& cellnum = mesh.cellnum_export();
    auto& cellcoef = cc;
    const auto& id_aP = mesh.id_aP;
    const auto& id_aE = mesh.id_aE;
    const auto& id_aW = mesh.id_aW;
    const auto& id_aN = mesh.id_aN;
    const auto& id_aS = mesh.id_aS;
    const auto& id_aT = mesh.id_aT;
    const auto& id_aB = mesh.id_aB;
    const auto& id_bsrc = mesh.id_bsrc;

    double max_norm = -1.e20;

    //求解参数
    const auto& nl_total_steps = solversettings.nlinsol_iter_steps;

    //定义与post对象相关的局部变量
    //残差数据输出频率
    const auto& res_freq = postprocess.ResfreqExport();
      
    float t_E, t_W, t_N, t_S, t_T, t_B;
    double initial_norm;

    for(int i=1; i<iter_num+1; ++i)
    {
        // 保存旧温度值用于计算残差
        auto t_old = mesh.t;

        //初始化范数norm2
        double norm2 = 0.0;

        for(int k=0; k < cellnum[2]; ++k)
        {
            for(int j=0; j < cellnum[1]; ++j)
            {
                for(int i=0; i < cellnum[0]; ++i)
                {
                    if(i == 0){
                        t_W = 0.0f; //防止数组越界，本来对应的ct系数就为0
                    }else{
                        // Gauss-Seidel: 使用已更新的值
                        t_W = mesh.t[i-1][j][k];
                    }
                    if(i == cellnum[0]-1){
                        t_E = 0.0f;
                    }else{
                        // Gauss-Seidel: 使用旧值（尚未更新）
                        t_E = t_old[i+1][j][k];
                    }

                    if(j == 0){
                        t_S = 0.0f; //防止数组越界，本来对应的ct系数就为0
                    }else{
                        // Gauss-Seidel: 使用已更新的值
                        t_S = mesh.t[i][j-1][k];
                    }
                    if(j == cellnum[1]-1){
                        t_N = 0.0f;
                    }else{
                        // Gauss-Seidel: 使用旧值（尚未更新）
                        t_N = t_old[i][j+1][k];
                    }

                    if(dim == 3){
                        if(k == 0){
                            t_T = 0.0f; //防止数组越界，本来对应的ct系数就为0
                        }else{
                            // Gauss-Seidel: 使用已更新的值
                            t_T = mesh.t[i][j][k-1];
                        }
                        if(k == cellnum[2]-1){
                            t_B = 0.0f;
                        }else{
                            // Gauss-Seidel: 使用旧值（尚未更新）
                            t_B = t_old[i][j][k+1];
                        }
                    }

                    float t_new = -(cellcoef[i][j][k][id_aE] * t_E +
                                cellcoef[i][j][k][id_aW] * t_W +
                                cellcoef[i][j][k][id_aN] * t_N +
                                cellcoef[i][j][k][id_aS] * t_S) + cellcoef[i][j][k][id_bsrc]; 
                    if(dim == 3){
                        t_new -=  cellcoef[i][j][k][id_aT] * t_T + cellcoef[i][j][k][id_aB] * t_B;
                    }

                    t_new = t_new / cellcoef[i][j][k][id_aP];

                    float du = rc * (t_new - t_old[i][j][k]);

                    // Gauss-Seidel: 立即更新值
                    mesh.t[i][j][k] = t_old[i][j][k] + du;
                    norm2 += du * du;
                }
            }
        }
        
        //网格数量
        int ncell = cellnum[0] * cellnum[1] * cellnum[2];
        //计算L2范数
        norm2 = double(std::sqrt(norm2 / ncell));

        if(i == 1){
            initial_norm = norm2 + 1.e-20;
        }

        max_norm = std::max(norm2, max_norm) + 1.e-20f;
         
        double rel_norm = norm2 / max_norm;

        if(rel_norm < torlence || i == iter_num){
            total_iter_num += i;
            if(nlitera_nums % res_freq || nlitera_nums == 1 || nlitera_nums == nl_total_steps){
                // 输出到控制台
                std::cout << "it_nl, it, tot_it, norm2, init, max, rel " 
                        << nlitera_nums << " " 
                        << i << " " 
                        << total_iter_num << " " 
                        << norm2 << " " 
                        << initial_norm << " " 
                        << max_norm << " " 
                        << rel_norm << std::endl;
                // 追加写入到文件
                std::ofstream linsol_fid(postprocess.linsol_fname, std::ios::app);  // std::ios::app 表示追加模式
                linsol_fid << nlitera_nums << " " 
                        << i << " " 
                        << total_iter_num << " " 
                        << norm2 << " " 
                        << initial_norm << " " 
                        << max_norm << " " 
                        << rel_norm << std::endl;
                linsol_fid.close();
            }
            break;
        } 
    }
    return;
}

void LinearSolver::JacobiIteraMethod(int& nlitera_nums,std::vector<std::vector<std::vector<std::vector<float>>>>& cc,
                                     StructureMesh& mesh, PostProcess& postprocess, 
                                     SolverSettings& solversettings, int& iter_num, int& total_iter_num, float& rc, float& torlence)
{
    // 定义与网格相关的变量
    const auto& dim = mesh.dim_export();
    const auto& cellnum = mesh.cellnum_export();
    auto& cellcoef = cc;
    const auto& id_aP = mesh.id_aP;
    const auto& id_aE = mesh.id_aE;
    const auto& id_aW = mesh.id_aW;
    const auto& id_aN = mesh.id_aN;
    const auto& id_aS = mesh.id_aS;
    const auto& id_aT = mesh.id_aT;
    const auto& id_aB = mesh.id_aB;
    const auto& id_bsrc = mesh.id_bsrc;

    double max_norm = -1.e20;

    //求解参数
    const auto& nl_total_steps = solversettings.nlinsol_iter_steps;

    //定义与post对象相关的局部变量
    //残差数据输出频率
    const auto& res_freq = postprocess.ResfreqExport();
      
    float t_E, t_W, t_N, t_S, t_T, t_B;
    double initial_norm;

    for(int i=1; i<iter_num+1; ++i)
    {
        //将当前版本复制到旧版本
        auto t0 = mesh.t;

        //初始化范数norm2
        double norm2 = 0.0;

        for(int k=0; k < cellnum[2]; ++k)
        {
            for(int j=0; j < cellnum[1]; ++j)
            {
                for(int i=0; i < cellnum[0]; ++i)
                {
                    if(i == 0){
                        t_W = 0.0f; //防止数组越界，本来对应的ct系数就为0
                    }else{
                        t_W = t0[i-1][j][k];
                    }
                    if(i == cellnum[0]-1){
                        t_E = 0.0f;
                    }else{
                        t_E = t0[i+1][j][k];
                    }

                    if(j == 0){
                        t_S = 0.0f; //防止数组越界，本来对应的ct系数就为0
                    }else{
                        t_S = t0[i][j-1][k];
                    }
                    if(j == cellnum[1]-1){
                        t_N = 0.0f;
                    }else{
                        t_N = t0[i][j+1][k];
                    }

                    if(dim == 3){
                        if(k == 0){
                            t_T = 0.0f; //防止数组越界，本来对应的ct系数就为0
                        }else{
                            t_T = t0[i][j][k-1];
                        }
                        if(i == cellnum[2]-1){
                            t_B = 0.0f;
                        }else{
                            t_B = t0[i][j][k+1];
                        }
                    }

                    float t_new = -(cellcoef[i][j][k][id_aE] * t_E +
                                cellcoef[i][j][k][id_aW] * t_W +
                                cellcoef[i][j][k][id_aN] * t_N +
                                cellcoef[i][j][k][id_aS] * t_S) + cellcoef[i][j][k][id_bsrc]; 
                    if(dim == 3){
                        t_new -=  cellcoef[i][j][k][id_aT] * t_T + cellcoef[i][j][k][id_aB] * t_B;
                    }

                    t_new = t_new / cellcoef[i][j][k][id_aP];

                    float du = rc * (t_new - t0[i][j][k]);

                    mesh.t[i][j][k] = t0[i][j][k] + du;
                    norm2 += du * du;
                }
            }
        }
        
        //网格数量
        int ncell = cellnum[0] * cellnum[1] * cellnum[2];
        //计算L2范数
        norm2 = double(std::sqrt(norm2 / ncell));

        if(i == 1){
            initial_norm = norm2 + 1.e-20;
        }

        max_norm = std::max(norm2, max_norm) + 1.e-20f;
         
        double rel_norm = norm2 / max_norm;

        if(rel_norm < torlence || i == iter_num){
            total_iter_num += i;
            if(nlitera_nums % res_freq || nlitera_nums == 1 || nlitera_nums == nl_total_steps){
                // 输出到控制台
                std::cout << "it_nl, it, tot_it, norm2, init, max, rel " 
                        << nlitera_nums << " " 
                        << i << " " 
                        << total_iter_num << " " 
                        << norm2 << " " 
                        << initial_norm << " " 
                        << max_norm << " " 
                        << rel_norm << std::endl;
                // 追加写入到文件
                std::ofstream linsol_fid(postprocess.linsol_fname, std::ios::app);  // std::ios::app 表示追加模式
                linsol_fid << nlitera_nums << " " 
                        << i << " " 
                        << total_iter_num << " " 
                        << norm2 << " " 
                        << initial_norm << " " 
                        << max_norm << " " 
                        << rel_norm << std::endl;
                linsol_fid.close();
            }
            break;
        } 
    }
    return;
}
