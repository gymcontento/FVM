/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-20 09:32:42
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-02-10 01:29:28
 * FilePath: \simple\src\kernel.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "kernel.h"
#include <algorithm>
#include <iostream>

std::vector<float> Kernel::CalcuAreaVolume(const float &dx, const float &dy, const float &dz)
{
    float area_x = dy*dz;
    float area_y = dx*dz;
    float area_z = dx*dy;
    float vol    = dx*area_x;

    std::vector<float> temp{area_x, area_y, area_z, vol};
    
    return temp;
}

float Kernel::CoefPecletPow(float pec){
    float ap = 1.0f - 0.1f * std::abs(pec);
    ap = std::max(0.0f, std::pow(ap, 5.0f));
    return ap;
}

float Kernel::CalCoefA(const int& conv_scheme, const float& area, const float &idx, 
                    const float &ul, const float &ur, const float &left_mul,
                    const float &right_mul, const float &rho, float sign_upwind)
{
    //需求解的场变量采用迎风格式，流速采用中心差分格式，未包含时间项
    //扩散项的系数在界面处的值根据调和平均进行计算
    float diffusion_coef = 2.0f * right_mul * left_mul / (right_mul + left_mul + 1e-12f) * idx;
    float convection_coef = 0.5f * rho * (ul + ur);
    if(conv_scheme == 0){
        a_coef = area * (diffusion_coef + std::max(0.0f, sign_upwind * convection_coef));
    }else if(conv_scheme == 1){
        a_coef = area * (diffusion_coef * (1.0f - 0.5f * std::abs(convection_coef / diffusion_coef))
                        + std::max(0.0f, sign_upwind * convection_coef));
    }else if(conv_scheme == 2){
        a_coef = area * (diffusion_coef * CoefPecletPow(convection_coef / diffusion_coef)
                        + std::max(0.0f, sign_upwind * convection_coef));
    }else if(conv_scheme == 3){
        // 简化SOU格式（不需要phi值，但精度略低）
        float base_term = diffusion_coef + std::max(0.0f, sign_upwind * convection_coef);
        
        // 计算贝克来数
        float peclet = std::abs(convection_coef / (diffusion_coef + 1e-12f));
        
        // 简化的二阶修正：基于贝克来数的修正
        float sou_correction = 0.0f;
        
        if(peclet < 2.0f) {
            // 小贝克来数区：接近中心差分
            sou_correction = 0.5f * std::abs(convection_coef) * (1.0f - 0.5f * peclet);
        } else {
            // 大贝克来数区：增强迎风效应
            sou_correction = 0.25f * std::abs(convection_coef);
        }
        
        a_coef = area * (base_term + sou_correction);
    }
    
    
    return a_coef;
}
