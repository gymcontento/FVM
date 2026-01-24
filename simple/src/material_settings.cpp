/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-04 23:07:39
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-11 17:00:41
 * FilePath: \FVM\heat_conduction\src\material_settings.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "material_settings.h"

void MaterialSettings::SetInitialDenstVist(StructureMesh& mesh, float& denst, float& vist)
{
    const std::vector<int>& cellnum = mesh.cellnum_export();
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
    const std::vector<int>& cellnum = mesh.cellnum_export();
    // std::cout << "---------Set Initial density---------" << std::endl;
    // for(int i=0; i < cellnum[2]; ++i)
    // {
    //     for(int j=0; j < cellnum[1]; ++j)
    //     {
    //         for(int k=0; k < cellnum[0]; ++k)
    //         {
    //             std::cout << density[k][j][i] << "\t";
    //         }
    //         std::cout << std::endl;
    //     }
    // }
    // std::cout << std::endl;
    // std::cout << "---------Set Initial viscosity---------" << std::endl;
    // for(int i=0; i < cellnum[2]; ++i)
    // {
    //     for(int j=0; j < cellnum[1]; ++j)
    //     {
    //         for(int k=0; k < cellnum[0]; ++k)
    //         {
    //             std::cout << viscosity[k][j][i] << "\t";
    //         }
    //         std::cout << std::endl;
    //     }
    // }
}
