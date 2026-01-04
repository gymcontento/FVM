#include "material_settings.h"

void MaterialSettings::SetInitialDenstVist(StructureMesh& mesh, float& denst, float& vist)
{
    std::vector<int> cellnum = mesh.cellnum_export();
    density = std::vector<std::vector<std::vector<float>>>(cellnum[0], 
        std::vector<std::vector<float>>(cellnum[1], std::vector<float>(cellnum[2], denst)));
    viscosity = std::vector<std::vector<std::vector<float>>>(cellnum[0], 
        std::vector<std::vector<float>>(cellnum[1], std::vector<float>(cellnum[2], vist))); 
}

void MaterialSettings::SetCondtSpcHt(float& condt, float& spcHt)
{
    conductivity = condt;
    specficheat = spcHt;
}

void MaterialSettings::CheckDensityViscosity(StructureMesh& mesh)
{
    std::vector<int> cellnum = mesh.cellnum_export();
    std::cout << "---------Set Initial density---------" << std::endl;
    for(int i=0; i < cellnum[2]; ++i)
    {
        for(int j=0; j < cellnum[1]; ++j)
        {
            for(int k=0; k < cellnum[0]; ++k)
            {
                std::cout << density[k][j][i] << "\t";
            }
            std::cout << std::endl;
        }
    }
    std::cout << std::endl;
    std::cout << "---------Set Initial viscosity---------" << std::endl;
    for(int i=0; i < cellnum[2]; ++i)
    {
        for(int j=0; j < cellnum[1]; ++j)
        {
            for(int k=0; k < cellnum[0]; ++k)
            {
                std::cout << viscosity[k][j][i] << "\t";
            }
            std::cout << std::endl;
        }
    }
}
