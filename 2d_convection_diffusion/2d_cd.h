#include <string>
#include <vector>

class Convection_Diffusion2d
{
public:
    void Init();
    void AssembleLinearSystem();
    void SolvingLinearSystem_LU();
    void SolvingLinearSystem_GS();
    void OutputResults();
    void OutputToFile(const std::string& filename);
private:
    int x_size;
    int y_size;
    double k;        // Thermal conductivity
    double velocity;
    double density;
    std::vector<double> b_val{0.0, 0.0, 0.0, 0.0};
    std::vector<double> range_x{0.0, 0.0};
    std::vector<double> range_y{0.0, 0.0};
    std::vector<std::vector<double>> A;
    std::vector<double> b;
    std::vector<double> T;  // Temperature solution vector
};