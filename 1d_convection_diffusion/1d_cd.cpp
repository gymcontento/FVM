#include "1d_cd.h"
#include <iostream>
#include <fstream>
#include <iomanip>

void Convection1d::Init(const std::string& param_file)
{
    std::ifstream ifs_(param_file);
    if(!ifs_.is_open())
    {
        std::cerr << "Error: Cannot open " << param_file << " file!" << std::endl;
        return;
    }

    std::string line;
    while(std::getline(ifs_,line))
    {
        size_t colonPos = line.find(':');
        if (colonPos != std::string::npos) {
            std::string key = line.substr(0, colonPos);
            std::string value = line.substr(colonPos + 1);

            key.erase(0, key.find_first_not_of(" \t\""));
            key.erase(key.find_last_not_of(" \t\"") + 1);
            value.erase(0, value.find_first_not_of(" \t\""));
            value.erase(value.find_last_not_of(" \t\",") + 1);

            if (key == "x_begin") range[0] = std::stod(value);
            else if (key == "x_end") range[1] = std::stod(value);
            else if (key == "x_num") x_num = std::stoi(value);
            else if (key == "density") density = std::stod(value);
            else if (key == "velocity") velocity = std::stod(value);
            else if (key == "thick") thick = std::stod(value);
            else if (key == "k") k = std::stod(value);     
            else if (key == "W_boundary") b_val[0] = std::stod(value);
            else if (key == "E_boundary") b_val[1] = std::stod(value);  
        }
    }
    ifs_.close();

    A.resize(x_num, std::vector<double>(x_num, 0.0));
    b.resize(x_num, 0.0);
    T.resize(x_num, 0.0);

    std::cout << "---------Init---------" << std::endl;
    std::cout << "Parameters loaded from parameter.json:" << std::endl;
    std::cout << "X range: [" << range[0] << ", " << range[1] << "]" << std::endl;
    std::cout << "Grid dimensions: " << x_num << " nodes" << std::endl;
    std::cout << "Matrix A size: " << x_num << " x " << x_num << std::endl;
    std::cout << "Vector b size: " << x_num << std::endl;
    std::cout << "Boundary conditions: E=" << b_val[0] << ", W=" << b_val[1] << std::endl;
    std::cout << "Thermal conductivity k = " << k << std::endl;
    std::cout << "Velocity = " << velocity << std::endl;
    std::cout << "Density = " << density << std::endl;
    std::cout << "Thickness = " << thick << std::endl;
}

void Convection1d::AssembleLinearSystem()
{
    double dx = (range[1] - range[0]) / x_num;

    double coe_f = density * velocity / 2.0;
    double coe_d = k / dx;

    for(int i=0; i < x_num; ++i)
    {
        if(i == 0)
        {
            A[i][i] = coe_f + 3.0 * coe_d;
            A[i][i+1] = -(coe_d - coe_f);
            b[i] = (2.0 * coe_f + 2.0 * coe_d) * b_val[0];
        }
        else if(i == x_num-1)
        {
            A[i][i] = -coe_f + 3.0 * coe_d;
            A[i][i-1] = -(coe_f + coe_d);
            b[i] = (2.0 * coe_f - 2.0 * coe_d) * b_val[1];
        }else
        {
            A[i][i] = 2.0 * coe_d;
            A[i][i-1] = -(coe_d + coe_f);
            A[i][i+1] = -(coe_d - coe_f);
        }
    }

    std::cout << "----------Assemble----------" << std::endl;
    std::cout << "Element Number: " << x_num << std::endl;
    std::cout << "Element size: " << dx << std::endl;

    std::cout << "Matrix A: " << std::endl; 
    for(int i=0; i < x_num; ++i)
    {
        for(int j=0; j < x_num; ++j)
        {
            std::cout << A[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << "b: " << std::endl;
    for(int i=0; i < x_num; ++i)
    {
        std::cout << b[i] << "\t";
    }
    std::cout << std::endl; 
}

void Convection1d::SolvingLinearSystem_GE()
{
    std::cout << "----------Gaussian Elimination----------" << std::endl;
    
    // Create copies of A and b for elimination
    std::vector<std::vector<double>> A_temp = A;
    std::vector<double> b_temp = b;
    std::vector<double> x(x_num, 0.0);
    
    // Forward elimination
    for(int k = 0; k < x_num - 1; ++k)
    {
        // Find pivot row
        int max_row = k;
        double max_val = std::abs(A_temp[k][k]);
        
        for(int i = k + 1; i < x_num; ++i)
        {
            if(std::abs(A_temp[i][k]) > max_val)
            {
                max_val = std::abs(A_temp[i][k]);
                max_row = i;
            }
        }
        
        // Check for singularity
        if(max_val < 1e-10)
        {
            std::cerr << "Error: Matrix is singular or nearly singular!" << std::endl;
            return;
        }
        
        // Swap rows if needed
        if(max_row != k)
        {
            std::swap(A_temp[k], A_temp[max_row]);
            std::swap(b_temp[k], b_temp[max_row]);
        }
        
        // Eliminate column k
        for(int i = k + 1; i < x_num; ++i)
        {
            double factor = A_temp[i][k] / A_temp[k][k];
            
            // Update row i
            for(int j = k; j < x_num; ++j)
            {
                A_temp[i][j] -= factor * A_temp[k][j];
            }
            b_temp[i] -= factor * b_temp[k];
        }
    }
    
    // Back substitution
    for(int i = x_num - 1; i >= 0; --i)
    {
        double sum = 0.0;
        for(int j = i + 1; j < x_num; ++j)
        {
            sum += A_temp[i][j] * x[j];
        }
        x[i] = (b_temp[i] - sum) / A_temp[i][i];
    }
    
    // Copy solution to T vector
    T = x;
    
    std::cout << "Gaussian elimination completed successfully." << std::endl;
    
    // Calculate and display residual for verification
    double max_residual = 0.0;
    for(int i = 0; i < x_num; ++i)
    {
        double residual = b[i];
        for(int j = 0; j < x_num; ++j)
        {
            residual -= A[i][j] * T[j];
        }
        max_residual = std::max(max_residual, std::abs(residual));
    }
    std::cout << "Maximum residual: " << max_residual << std::endl;
}

void Convection1d::SolvingLinearSystem_GS()
{
    std::cout << "----------SolvingLinearSystem----------" << std::endl;
    int max_iterations = 1000;
    double tolerance = 0.00001;
    
    // Initialize T vector with zeros or initial guess
    for(int i = 0; i < x_num; ++i)
    {
        T[i] = 0.0;
    }
    
    double max_error = 0.0;
    for(int iter = 0; iter < max_iterations; ++iter)
    {
        max_error = 0.0;
        
        for(int i = 0; i < x_num; ++i)
        {
            double T_old = T[i];
            double sum = 0.0;
            
            // Gauss-Seidel: use updated values for j < i, old values for j > i
            for(int j = 0; j < x_num; ++j)
            {
                if(j != i)
                {
                    sum += A[i][j] * T[j];
                }
            }
            
            // Update T[i] using Gauss-Seidel formula
            T[i] = (b[i] - sum) / A[i][i];
            
            // Calculate error
            double error = std::abs(T[i] - T_old);
            max_error = std::max(max_error, error);
        }

        if(max_error < tolerance)
        {
            std::cout << "Gauss-Seidel converged in " << iter + 1 << " iterations." << std::endl;
            std::cout << "Final residual: " << max_error << std::endl;
            return;
        }

        // Print progress every 100 iterations
        if((iter + 1) % 100 == 0) {
            std::cout << "Iteration " << iter + 1 << ", max error: " << max_error << std::endl;
        }
    }
    
    std::cout << "Gauss-Seidel did not converge after " << max_iterations << " iterations." << std::endl;
    std::cout << "Final residual: " << max_error << std::endl;  
}

void Convection1d::OutputToFile(const std::string& filename)
{
    std::ofstream ofs(filename);
    if(!ofs.is_open())
    {
        std::cerr << "Error: Cannot open " << filename << " for writing!" << std::endl;
        return;
    }
    
    // Write parameters
    ofs << "# Parameters" << std::endl;
    ofs << "x_num=" << x_num << std::endl;
    ofs << "velocity=" << velocity << std::endl;
    ofs << "density=" << density << std::endl;
    ofs << "k=" << k << std::endl;
    ofs << "thick=" << thick << std::endl;
    ofs << "x_begin=" << range[0] << std::endl;
    ofs << "x_end=" << range[1] << std::endl;
    ofs << "W_boundary=" << b_val[0] << std::endl;
    ofs << "E_boundary=" << b_val[1] << std::endl;
    ofs << std::endl;
    
    // Write numerical results
    ofs << "# Numerical Results" << std::endl;
    double dx = (range[1] - range[0]) / x_num;
    for(int i = 0; i < x_num; ++i)
    {
        double x = range[0] + (i + 0.5) * dx;  // Cell center
        ofs << x << "\t" << T[i] << std::endl;
    }
    
    ofs.close();
}

void Convection1d::OutputResults()
{
    std::cout << "---------Solution Results---------" << std::endl;
    std::cout << "Temperature distribution:" << std::endl;
    
    // Output temperature field in grid format
    for(int i = 0; i < x_num; ++i) {

        std::cout << "T[" << i << "] = " << T[i] << "\t";
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
    
    // Calculate and output some statistics
    double T_min = T[0];
    double T_max = T[0];
    double T_avg = 0.0;
    
    for(int i = 0; i < x_num; ++i) {
            T_min = std::min(T_min, T[i]);
            T_max = std::max(T_max, T[i]);
            T_avg += T[i];
    }
    T_avg /= x_num;
    
    std::cout << "Temperature statistics:" << std::endl;
    std::cout << "Minimum temperature: " << T_min << std::endl;
    std::cout << "Maximum temperature: " << T_max << std::endl;
    std::cout << "Average temperature: " << T_avg << std::endl;
    std::cout << std::endl;
    
    // Output corner temperatures for verification
    std::cout << "Corner temperatures:" << std::endl;
    std::cout << "Bottom-left (SW): " << T[0] << std::endl;
    std::cout << "Bottom-right (SE): " << T[x_num - 1] << std::endl;
    
    // Verify solution by computing residual
    std::cout << "\nSolution verification:" << std::endl;
    double max_residual = 0.0;
    for(int i = 0; i < x_num; ++i) {
        double residual = b[i];
        for(int j = 0; j < x_num; ++j) {
            residual -= A[i][j] * T[j];
        }
        max_residual = std::max(max_residual, std::abs(residual));
        std::cout << "Residual[" << i << "] = " << residual << std::endl;
    }
    std::cout << "Max residual: " << max_residual << std::endl;
}
