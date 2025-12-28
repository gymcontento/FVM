/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2025-12-28 12:18:51
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2025-12-29 00:08:04
 * FilePath: \FVM\1d_convection_diffusion\1d_cd.h
 * Description: 
 * 
 * Copyright (c) 2025 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include <string>
#include <vector>

class Convection1d
{
public:
    void Init(const std::string& param_file = "parameters.json");
    void AssembleLinearSystem();
    void SolvingLinearSystem_GE();
    void SolvingLinearSystem_GS();
    void OutputResults();
    void OutputToFile(const std::string& filename);
private:
    int x_num;
    double density;
    double velocity;
    double thick;
    double k;
    std::vector<double> range{0.0,0.0};
    std::vector<double> b_val{0.0,0.0};
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> T; 
};
