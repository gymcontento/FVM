/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-04 22:46:57
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-10 10:48:41
 * FilePath: \FVM\heat_conduction\src\solver_settings.cpp
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

void SolverSettings::SetTimeStep(float& timestep)
{
    dt = timestep;
}

void SolverSettings::SetThermalSolverParam(int& t_iter_nums_s, float& relax_coef_s, 
        float& linsol_residual_s, float& nlinsol_residual_s)
{
    t_iter_nums = t_iter_nums_s;
    relax_coef = relax_coef_s;
    linsol_residual = linsol_residual_s;
    nlinsol_residual = nlinsol_residual_s;
}

void SolverSettings::CheckSolverSettings()
{
    std::cout << "---------- Solver Settings ----------" << std::endl;
    std::cout << "Nonlinear Steps: " << nlinsol_iter_steps << std::endl;
    std::cout << "Linear Iteration Nums: " << t_iter_nums << std::endl;
    std::cout << "Relax Coefficient: " << relax_coef << std::endl;
    std::cout << "linear residual: " << linsol_residual << std::endl;
    std::cout << "Nonlinear Residual: " << nlinsol_residual << std::endl;
    std::cout << std::endl;
}
