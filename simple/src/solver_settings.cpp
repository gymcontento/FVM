/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-04 22:46:57
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-02-02 01:11:59
 * FilePath: \simple\src\solver_settings.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "solver_settings.h"
#include <iostream>

void SolverSettings::SetEqnType(std::string& eqntype)
{
    if(eqntype == "conduction"){
        eqn_type = 0;
    }else if(eqntype == "flow"){
        eqn_type = 1;
    }else if(eqntype == "conduction_flow"){
        eqn_type = 2;
    }
}

void SolverSettings::SetIterNums(int &steps)
{
    nlinsol_iter_steps = steps;
}

void SolverSettings::SetDeltaT(float& timestep)
{
    dt = timestep;
}


void SolverSettings::SetSolverParam(int& nlin_steps, float& nlin_res, std::vector<int>& velocity_iter_nums,
        std::vector<float>& velocity_relax, std::vector<float>& velocity_tol,
        int& temp_iter_num, float& temp_relax, float& temp_tol,
        int& press_iter_num, float& press_relax, float& press_tol,
        float& mass_tol, float& time_tol)
{
    nlinsol_iter_steps = nlin_steps;
    nlinsol_residual = nlin_res;
    
    lin_iters_nums_u = velocity_iter_nums[0];
    lin_iters_nums_v = velocity_iter_nums[1];
    lin_iters_nums_w = velocity_iter_nums[2];
    lin_iters_nums_t = temp_iter_num;
    lin_iters_nums_p = press_iter_num;

    relax_coef_u = velocity_relax[0];
    relax_coef_v = velocity_relax[1];
    relax_coef_w = velocity_relax[2];
    relax_coef_p = press_relax;
    relax_coef_t = temp_relax;

    lin_residual_u = velocity_tol[0];
    lin_residual_v = velocity_tol[1];
    lin_residual_w = velocity_tol[2];
    lin_iters_nums_p = press_tol;
    lin_iters_nums_t = temp_tol;

    nlinsol_residual_mass = mass_tol;
    nlinsol_residual_temp = temp_tol;
}

void SolverSettings::CheckSolverSettings()
{
    std::cout << "---------- Solver Settings ----------" << std::endl;
    std::cout << "Conv Scheme: " << conv_scheme << std::endl;
    std::cout << "Nonlinear Steps: " << nlinsol_iter_steps << std::endl;
    // std::cout << "Linear Iteration Nums: " << lin_iter_steps << std::endl;
    // std::cout << "Relax Coefficient: " << relax_coef << std::endl;
    // std::cout << "linear residual: " << linsol_residual << std::endl;
    std::cout << "Nonlinear Residual: " << nlinsol_residual << std::endl;
    std::cout << std::endl;
}

void SolverSettings::SetConvScheme(std::string &str_conv_scheme)
{
    if(str_conv_scheme == "uw" || str_conv_scheme == "UW"){
        conv_scheme = conv_scheme_upwind;
    }else if(str_conv_scheme == "cd" || str_conv_scheme == "CD"){
        conv_scheme = conv_scheme_cd;
    }else if(str_conv_scheme == "pl" || str_conv_scheme == "PL"){
        conv_scheme = conv_scheme_powerlaw;
    }else if(str_conv_scheme == "sou" || str_conv_scheme == "SOU"){
        conv_scheme = conv_scheme_sou;
    }
    std::cout << "---------- Check Conv Scheme ----------" << "\n";
    std::cout << "Conv Scheme: " << conv_scheme << "\n";
}

