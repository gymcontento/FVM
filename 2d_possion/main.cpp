#include "2d_possion.h"

int main()
{
    Possion2d po2d;
    po2d.Init();
    po2d.AssembleLinearSystem();
    
    // Test LU decomposition solver
    po2d.SolvingLinearSystem_LU();
    po2d.OutputResults();
    
    std::cout << "\n" << std::string(50, '=') << "\n" << std::endl;
    
    // Test Gauss-Seidel solver for comparison
    po2d.SolvingLinearSystem_GS();
    po2d.OutputResults();
    
    return 0;
}
