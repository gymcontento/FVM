/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-15 09:18:29
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-15 13:21:02
 * FilePath: \FVM\convection_diffusion\include\kernel.h
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#ifndef KERNEL_H
#define KERNEL_H
#include <vector>

class Kernel
{
public:
    std::vector<float> CalcuAreaVolume(const float& dx, const float& dy, const float& dz);
    float CalCoefA(const int& conv_scheme, const float& area, const float& dx,
                const float& ul, const float& ur,
                const float& right_cond, const float& left_cond,
                const float& rho, float sign_upwind);
private:
    float a_coef;
};

#endif