/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2025-12-27 21:37:44
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2025-12-28 12:03:38
 * FilePath: \FVM\2d_possion\main.cpp
 * Description: 
 * 
 * Copyright (c) 2025 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "2d_possion.h"

int main()
{
    Possion2d po2d;
    po2d.Init();
    po2d.AssembleLinearSystem();
    
    // Test LU decomposition solver
    po2d.SolvingLinearSystem_LU();
    po2d.OutputResults();
    
    std::cout << "\n" << std::string(50, '=') << "\n" << std::endl;
    
    // Test Gauss-Seidel solver for comparison
    po2d.SolvingLinearSystem_GS();
    po2d.OutputResults();
    
    // Output Gauss-Seidel results to file for Python visualization
    po2d.OutputToFile("gauss_seidel_results.csv");
    
    return 0;
}
