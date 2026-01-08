#ifndef KERNEL_H
#define KERNEL_H
#include <vector>

class Kernel
{
public:
    std::vector<float> CalcuAreaVolume(const float& dx, const float& dy, const float& dz);
    float CalCoefA(const float& area, const float& dx,
                const float& ul, const float& ur,
                const float& right_cond, const float& left_cond,
                const float& rho, float sign_upwind);
};

#endif