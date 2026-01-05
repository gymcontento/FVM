#include "boundarysettings.h"
#include <cmath>

void BoundarySettings::SetBcFacesForCells(StructureMesh &mesh)
{
    const int& dim = mesh.dim_export();
    const std::vector<std::vector<float>>& range = mesh.range_export();
    const std::vector<int>& cellnum = mesh.cellnum_export();
    const std::vector<float>& x_nodes = mesh.xnodes_export();
    const std::vector<float>& y_nodes = mesh.ynodes_export();
    const std::vector<float>& z_nodes = mesh.znodes_export();

    //存储单元边界条件id，若为内部点，则设置为0
    //顺序为 e w n s t b
    bcid = std::vector<std::vector<std::vector<std::vector<int>>>>(
        cellnum[0], std::vector<std::vector<std::vector<int>>>(cellnum[1],
        std::vector<std::vector<int>>(cellnum[2], std::vector<int>(2 * dim, 0))));

    //设置检查容差
    float eps = 1.e-12f;
    float x0, y0, z0;

    //判断每个网格的边界是否是在外边界上，并且记录到bcid中
    for(int k=0; k < cellnum[2]; ++k)
    {
        for(int j=0; j < cellnum[1]; ++j)
        {
            for(int i=0; i < cellnum[0]; ++i)
            {
                // 每个网格的西边界
                if (std::abs(x_nodes[i]-range[0][0]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                {   
                    bcid[i][j][k][fid_w] = bcid_xmin;
                }
                // 每个网格的东边界
                if (std::abs(x_nodes[i+1]-range[0][1]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                {   
                    bcid[i][j][k][fid_e] = bcid_xmax;
                }
                // 每个网格的南边界
                if (std::abs(y_nodes[i]-range[1][0]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                {   
                    bcid[i][j][k][fid_s] = bcid_ymin;
                }
                // 每个网格的北边界
                if (std::abs(y_nodes[i+1]-range[1][1]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                {   
                    bcid[i][j][k][fid_n] = bcid_ymax;
                }
                //如果维度为三维，要多加两个边界判断
                    // if dim == 3:
                   // 每个网格的底边界
                if(dim == 3)
                {
                    //每个网格的底边界
                    if (std::abs(z_nodes[i]-range[2][0]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                    {   
                        bcid[i][j][k][fid_b] = bcid_zmax;
                    }
                    //每个网格的顶边界
                    if (std::abs(z_nodes[i+1]-range[2][1]) < eps) // 如果差值的绝对值极小，就判断在西边界上
                    {   
                        bcid[i][j][k][fid_t] = bcid_zmax;
                    }
                }
            }
        }
    }
}
