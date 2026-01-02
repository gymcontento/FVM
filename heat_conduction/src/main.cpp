/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-03 00:35:43
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-03 01:27:53
 * FilePath: \FVM\heat_conduction\src\main.cpp
 * Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
 ******************************/
#include "structure_mesh.h"
#include <ctime>
#include <vector>
#include <iomanip>

int main(){
    clock_t start = clock();
    int dim = 2;
    std::vector<int> numcell{64,64,1};
    
    float xmin = 0.0;
    float xmax = 0.833;
    float ymin = 0.0;
    float ymax = 0.833;
    float zmin = 0.0;
    float zmax = 1.0;

    StructureMesh case1;
    case1.CreateMesh(dim,numcell);    

    clock_t end = clock();
    double duration = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << 
        std::fixed << std::setprecision(10) 
        << duration << " seconds" << std::endl;
}