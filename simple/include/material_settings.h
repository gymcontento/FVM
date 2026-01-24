#ifndef MATERIAL_SETTINGS_H
#define MATERIAL_SETTINGS_H
#include "structure_mesh.h"

class MaterialSettings
{
private:
    //密度和黏度
    std::vector<std::vector<std::vector<float>>> density;
    std::vector<std::vector<std::vector<float>>> viscosity;

    //传热系数和比热系数
    float conductivity{0.0f};
    float specficheat{0.0f};

    //热源
    float heat_source{0.0f};
public:
    void SetInitialDenstVist(StructureMesh& mesh, float& denst, float& vist);
    void SetCondtSpcHt(float& condt, float& spcHt);

    const std::vector<std::vector<std::vector<float>>>& DensityExport() {return density;}
    const std::vector<std::vector<std::vector<float>>>& ViscosityExport() {return viscosity;}

    //数据导出函数
    const float& ConductivityExport() {return conductivity;}
    const float& SpecficheatExport() {return specficheat;}
    const float& HeatSourceExport() {return heat_source;}
    //查看数据函数
    void CheckDensityViscosity(StructureMesh& mesh);
};

#endif