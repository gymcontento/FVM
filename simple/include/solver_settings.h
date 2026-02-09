/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-20 09:32:42
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-02-02 01:12:16
 * FilePath: \simple\include\solver_settings.h
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#ifndef SOLVER_SETTINGS_H
#define SOLVER_SETTINGS_H
#include <string>
#include <vector>

class SolverSettings
{
public:
    //方程类型
    int eqn_type{0};
    int eqn_type_conduct{0};
    int eqn_type_flow{1};
    int eqn_type_conduct_flow{2};

    //对流项离散格式类型
    int conv_scheme{0};
    int conv_scheme_upwind{0};
    int conv_scheme_cd{1};
    int conv_scheme_powerlaw{2};
    int conv_scheme_sou{3};

    //非线性迭代次数
    int nlinsol_iter_steps{1};
    //非线性迭代的收敛误差判定值
    float nlinsol_residual{1.e-6f};
    float nlinsol_residual_mass{1.e-6f};
    float nlinsol_residual_temp{1.e-6f};

    //非线性迭代残差的l2范数
    float l2_curr{0.0};
    float l2_max{-1.e20f};
    float l2_max_u{-1.e20f};
    float l2_max_v{-1.e20f};
    float l2_max_w{-1.e20f};
    float l2_max_p{-1.e20f};
    float l2_max_pp{-1.e20f};
    float l2_max_t{-1.e20f};
    float l2_u{0.0f};
    float l2_v{0.0f};
    float l2_w{0.0f};
    float l2_p{0.0f};
    float l2_pp{0.0f};
    float l2_t{0.0f};

    //时间步长
    float dt{1.0};

    //回路控制变量，用来判断是否可以跳出循环
    bool stop_iter = false;

    //线性求解迭代次数
    // int lin_iter_steps{10};

    //松弛因子
    float relax_coef_u{0.0f};
    float relax_coef_v{0.0f};
    float relax_coef_w{0.0f};
    float relax_coef_p{0.0f};
    float relax_coef_t{0.75f};
    //线性求解器迭代的收敛误差判定值
    float lin_residual_u{1.e-2f};
    float lin_residual_v{1.e-2f};
    float lin_residual_w{1.e-2f};
    float lin_residual_p{1.e-2f};
    float lin_residual_t{1.e-2f};
    //线性迭代的总次数
    int lin_iters_nums_u{100};
    int lin_iters_nums_v{100};
    int lin_iters_nums_w{100};
    int lin_iters_nums_p{100};
    int lin_iters_nums_t{10};
    int total_lin_iters{0};

    void SetEqnType(std::string& eqntype);
    void SetConvScheme(std::string& conv_scheme);
    void SetIterNums(int& steps);
    void SetDeltaT(float& timestep);
    void SetSolverParam(int& nlin_steps, float& nlin_res, std::vector<int>& velocity_iter_nums,
        std::vector<float>& velocity_relax, std::vector<float>& velocity_tol,
        int& temp_iter_num, float& temp_relax, float& temp_tol,
        int& press_iter_num, float& press_relax, float& press_tol,
        float& mass_tol, float& time_tol);
    void CheckSolverSettings();
};

#endif