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

float Kernel::CalCoefA(const float& area, const float &dx, 
                    const float &ul, const float &ur, const float &right_cond,
                    const float &left_cond, const float &rho, float sign_upwind)
{
    //需求解的场变量采用迎风格式，流速采用中心差分格式，未包含时间项
    //扩散项的系数在界面处的值根据调和平均进行计算
    float diffusion_coef = 2.0f * right_cond * left_cond / ((right_cond + left_cond + 1e-12f) * dx);
    float convection_coef = 0.5f * rho * (ul + ur);
    float a_coef = area * (diffusion_coef + std::max(0.0f, sign_upwind * convection_coef));
    return a_coef;
}
