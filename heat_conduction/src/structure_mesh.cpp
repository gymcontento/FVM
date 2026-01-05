#include "structure_mesh.h"

void StructureMesh::CreateMesh(int& dimensions, std::vector<int>& nums)
    {
        // Number of Spatial dimensions
        dim = dimensions;

        // Number of cells in each direction
        cellnum[0] = nums[0];
        cellnum[1] = nums[1];
        cellnum[2] = nums[2];

        // Number of nodes in each direction
        nodenum[0] = nums[0] + 1;
        nodenum[1] = nums[1] + 1;
        nodenum[2] = nums[2] + 1;

        std::cout << "Number of cells in x, y, z = " << cellnum[0] << ", "
            << cellnum[1] << ", " << cellnum[2] << std::endl;
        std::cout << "Number of nodes in x, y, z = " << nodenum[0] << ", "
            << nodenum[1] << ", " << nodenum[2] << std::endl;
    }

void StructureMesh::CreateCoordinates(std::vector<std::vector<float>>& irange)
{
    range = irange;

    if(dim == 2){
        range[2][0] = 0.0f;
        range[2][1] = 1.0f;
    }

    //Mesh coordinates
    x_nodes.resize(nodenum[0], 0.0f);
    y_nodes.resize(nodenum[1], 0.0f);
    z_nodes.resize(nodenum[2], 0.0f);

    //Cell coordinates
    x_cells.resize(cellnum[0], 0.0f);
    y_cells.resize(cellnum[1], 0.0f);
    z_cells.resize(cellnum[2], 0.0f);

    //Mesh Generation
    dx = (range[0][1] - range[0][0]) / float(cellnum[0]);
    dy = (range[1][1] - range[1][0]) / float(cellnum[1]);
    dz = (range[2][1] - range[2][0]) / float(cellnum[2]);
    //Nodes coordinates generation
    for(int i=0; i < nodenum[0]; ++i)
    {
        x_nodes[i] = range[0][0] + i * dx;
    }
    for(int i=0; i < nodenum[1]; ++i)
    {
        y_nodes[i] = range[1][0] + i * dy;
    }
    for(int i=0; i < nodenum[2]; ++i)
    {
        z_nodes[i] = range[2][0] + i * dz;
    }
    //Cell coordinates generation 
    for(int i=0; i < cellnum[0]; ++i)
    {
        x_cells[i] = 0.5f * (x_nodes[i] + x_nodes[i+1]);
    }
    for(int i=0; i < cellnum[1]; ++i)
    {
        y_cells[i] = 0.5f * (y_nodes[i] + y_nodes[i+1]);
    }
    for(int i=0; i < cellnum[2]; ++i)
    {
        z_cells[i] = 0.5f * (z_nodes[i] + z_nodes[i+1]);
    }

    std::cout << "----------Coordinates Info----------" << std::endl;
    std::cout << "dx, dy, dz: " << dx << ", " << dy << ", "
        << dz << std::endl;
    std::cout <<"xmin, xmax: " << range[0][0] << ", " << range[0][1] << std::endl;
    std::cout <<"ymin, ymax: " << range[1][0] << ", " << range[1][1] << std::endl;
    std::cout <<"zmin, zmax: " << range[2][0] << ", " << range[2][1] << std::endl;
    std::cout << "X coordinates: " << std::endl;
    for(int i=0; i < nodenum[0]; ++i)
    {
        std::cout << x_nodes[i] << "\t";
    }
    std::cout << std::endl;
    std::cout << "Y coordinates: " << std::endl;
    for(int i=0; i < nodenum[1]; ++i)
    {
        std::cout << y_nodes[i] << "\t";
    }
    std::cout << std::endl;
    std::cout << "Z coordinates: " << std::endl;
    for(int i=0; i < nodenum[2]; ++i)
    {
        std::cout << z_nodes[i] << "\t";
    }
    std::cout << std::endl;
}

void StructureMesh::CreateFieldMeshData()
{
    t = std::vector<std::vector<std::vector<float>>>(cellnum[0], 
        std::vector<std::vector<float>>(cellnum[1], std::vector<float>(cellnum[2], 0.0f)));
    t0 = std::vector<std::vector<std::vector<float>>>(cellnum[0], 
        std::vector<std::vector<float>>(cellnum[1], std::vector<float>(cellnum[2], 0.0f)));
    
    std::cout << "---------Create Field Data---------" << std::endl;
    for(int i=0; i < cellnum[2]; ++i)
    {
        for(int j=0; j < cellnum[1]; ++j)
        {
            for(int k=0; k < cellnum[0]; ++k)
            {
                std::cout << t[k][j][i] << "\t";
            }
            std::cout << std::endl;
        }
    }
}

void StructureMesh::SetInitialT(float& T)
{
    t = std::vector<std::vector<std::vector<float>>>(cellnum[0], 
        std::vector<std::vector<float>>(cellnum[1], std::vector<float>(cellnum[2], T)));
    std::cout << "---------Set Initial T---------" << std::endl;
    for(int i=0; i < cellnum[2]; ++i)
    {
        for(int j=0; j < cellnum[1]; ++j)
        {
            for(int k=0; k < cellnum[0]; ++k)
            {
                std::cout << t[k][j][i] << "\t";
            }
            std::cout << std::endl;
        }
    }
}

void StructureMesh::CreateCoeffMeshData()
{
    //Set number of coefficients
    if(dim == 2)
    {
        numcoef = 6;
    }else if(dim == 3)
    {
        numcoef = 8;
    }

    // Coefficient storage positions in four dimensional array
    id_ap = 0;
    id_ae = 1;
    id_aw = 2;
    id_an = 3;
    id_as = 4;
    if(dim == 3){
        id_aT = 5;
        id_aB = 6;
    }
    id_bsrc = numcoef - 1;
    
    std::cout << "----------Set Coeff----------" << std::endl;
    std::cout << "Num of Coefficient: " << numcoef << std::endl;

    cellcoeff = std::vector<std::vector<std::vector<std::vector<float>>>>(
        cellnum[0], std::vector<std::vector<std::vector<float>>>(cellnum[1],
        std::vector<std::vector<float>>(cellnum[2], std::vector<float>(numcoef, 0))));
}

void StructureMesh::CreateSimulationData()
{

}
