/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2025-12-28 12:20:15
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2025-12-29 00:08:23
 * FilePath: \FVM\1d_convection_diffusion\main.cpp
 * Description: 
 * 
 * Copyright (c) 2025 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "1d_cd.h"
#include <iostream>
#include <vector>
#include <string>

int main(){
    std::vector<std::string> param_files = {
        "case1_params.json",
        "case2_params.json", 
        "case3_params.json",
        "case4_params.json"
    };
    
    std::vector<std::string> output_files = {
        "case1_results.txt",
        "case2_results.txt",
        "case3_results.txt", 
        "case4_results.txt"
    };
    
    std::vector<std::string> case_names = {
        "Case 1: 5 nodes, v=0.1",
        "Case 2: 5 nodes, v=2.5",
        "Case 3: 20 nodes, v=0.1", 
        "Case 4: 20 nodes, v=2.5"
    };
    
    for(int i = 0; i < param_files.size(); ++i) {
        std::cout << "\n" << std::string(80, '=') << std::endl;
        std::cout << case_names[i] << std::endl;
        std::cout << std::string(80, '=') << std::endl;
        
        Convection1d cd1d;
        cd1d.Init(param_files[i]);
        cd1d.AssembleLinearSystem();
        
        // Use Gaussian Elimination for accurate results
        cd1d.SolvingLinearSystem_GE();
        cd1d.OutputResults();
        
        // Save results to file for Python visualization
        cd1d.OutputToFile(output_files[i]);
        
        std::cout << "Results saved to: " << output_files[i] << std::endl;
    }
    
    return 0;
}
