/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2025-12-27 21:35:53
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2025-12-28 00:30:25
 * FilePath: \FVM\2d_possion\2d_possion.h
 * Description: 
 * 
 * Copyright (c) 2025 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include <iostream>
#include <vector>

class Possion2d{
public:
    void Init();
    void AssembleLinearSystem();
    void SolvingLinearSystem_LU();
    void SolvingLinearSystem_GS();
    void OutputResults();
private:
    int x_size;
    int y_size;
    double k;        // Thermal conductivity
    double thick;
    std::vector<double> b_val{0.0, 0.0, 0.0, 0.0};
    std::vector<double> range{0.0, 0.0};
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> T;  // Temperature solution vector
};
