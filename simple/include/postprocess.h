/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-04 13:21:22
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-02-04 17:41:26
 * FilePath: \simple\include\postprocess.h
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#ifndef POSTPROCESS_H
#define POSTPROCESS_H
#include <string>
#include "structure_mesh.h"
#include "material_settings.h"
#include "solver_settings.h"

class PostProcess
{
private:
    int resfreq, outfreq;
public:
    //输出文件名
    std::string linsol_fname{"lin.res"};
    std::string nonlinsol_fname{"nonlin.res"};
    std::string vtk_fname_temp{"post_temp.vtk"};

    void PostProcess::WriteVTKCollocated_temp(StructureMesh& mesh, std::string filename = "post_temp.vtk");
    void WriteVTKCollocated_vel(StructureMesh& mesh, std::string filename = "post_vel.vtk");
    void PostProcess::WriteVTKCollocated_pre(StructureMesh& mesh, std::string filename = "post_pre.vtk");
    // void WriteVTKCollocated_temp_Pe_L(StructureMesh& mesh, MaterialSettings& material);
    // void WriteVTKCollocated_temp_Pe_L_center(StructureMesh& mesh, MaterialSettings& material, SolverSettings& solversettings);
    
    void SetOuputFreq(const int& res_freq, const int& out_freq);

    const int& ResfreqExport() {return resfreq;}
    const int& OutfreqExport() {return outfreq;}
};

#endif