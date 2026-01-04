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

    std::vector<std::vector<std::vector<float>>> DensityExport() {return density;}
    std::vector<std::vector<std::vector<float>>> ViscosityExport() {return viscosity;}

    float ConductivityExport() {return conductivity;}
    float SpecficheatExport() {return specficheat;}

    float HeatSourceExport() {return heat_source;}

    void CheckDensityViscosity(StructureMesh& mesh);
};

#endif