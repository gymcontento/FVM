/******************************
 * Author: gymcontento herry996341591@gmail.com
 * Date: 2026-01-03 00:35:43
 * LastEditors: gymcontento herry996341591@gmail.com
 * LastEditTime: 2026-01-03 21:10:58
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
    std::vector<int> numcell{4,3,1};
    std::vector<std::vector<float>> range = {
        {0.0f, 1.0f}, {0.0, 1.0f}, {0.0, 1.0f}};
    float initalT = 100.0;

    StructureMesh case1;
    case1.CreateMesh(dim, numcell);
    case1.CreateCoordinates(range);
    case1.CreateFieldMeshData();
    case1.SetInitialT(initalT);
    case1.CreateCoeffMeshData();

    clock_t end = clock();
    double duration = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Elapsed Time: " << 
        std::fixed << std::setprecision(10) 
        << duration << " seconds" << std::endl;
}