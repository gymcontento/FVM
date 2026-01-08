/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-04 13:21:22
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-06 23:07:21
 * FilePath: \FVM\heat_conduction\include\postprocess.h
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#ifndef POSTPROCESS_H
#define POSTPROCESS_H
#include <string>
#include "structure_mesh.h"

class PostProcess
{
public:
    void WriteVTKCollocated_temp(
        StructureMesh& mesh, std::string filename = "post_temp.vtk");
    
    void SetOuputFreq(const int res_freq, const int );
};

#endif