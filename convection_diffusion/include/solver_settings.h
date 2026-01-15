#ifndef SOLVER_SETTINGS_H
#define SOLVER_SETTINGS_H
#include <string>

class SolverSettings
{
public:
    //方程类型
    int eqn_type{0};
    int eqn_type_conduct{0};
    int eqn_type_flow{1};
    int eqn_type_conduct_flow{2};

    //对流项离散格式类型
    int conv_scheme{0};
    int conv_scheme_upwind{0};
    int conv_scheme_cd{1};
    int conv_scheme_powerlaw{2};
    int conv_scheme_sou{3};

    //非线性迭代次数
    int nlinsol_iter_steps{1};
    //非线性迭代的收敛误差判定值
    float nlinsol_residual {1.e-6f};

    //非线性迭代残差的l2范数
    float l2_curr  = 0.0;
    float l2_max   = -1.e20f;
    float l2_max_t = -1.e20f;
    float l2_t     = 0.0f;

    //时间步长
    float dt{1.0};

    //回路控制变量，用来判断是否可以跳出循环
    bool stop_iter = false;

    //温度线性求解迭代次数
    int t_iter_nums{10};

    //松弛因子
    float relax_coef{0.75f};
    //线性求解器迭代的收敛误差判定值
    float linsol_residual{1.e-2f};
    //线性迭代的总次数
    int linsol_iters_nums{0};

    void SetEqnType(std::string& eqntype);
    void SetConvScheme(std::string& conv_scheme);
    void SetIterNums(int& steps);
    void SetTimeStep(float& timestep);
    void SetThermalSolverParam(int& t_iter_nums_s, float& relax_coef_s, 
        float& linsol_residual_s, float& nlinsol_residual_s);
    void CheckSolverSettings();
};

#endif