#include "2d_cd.h"

int main()
{
    Convection_Diffusion2d cd2d;
    cd2d.Init();
    cd2d.AssembleLinearSystem();
}