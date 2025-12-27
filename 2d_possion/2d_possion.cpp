#include "2d_possion.h"
#include <fstream>
#include <sstream>
#include <string>

void Possion2d::Init()
{
    std::ifstream ifs_("parameter.json");
    if (!ifs_.is_open()) {
        std::cerr << "Error: Cannot open parameter.json file!" << std::endl;
        return;
    }
    
    // Variables to store grid dimensions and coefficients
    int x_num = 0, y_num = 0;
    double A_coeff = 0.0, k_coeff = 0.0;
    
    // Read and parse JSON file
    std::string line;
    while (std::getline(ifs_, line)) {
        // Find the colon position to separate key and value
        size_t colonPos = line.find(':');
        if (colonPos != std::string::npos) {
            std::string key = line.substr(0, colonPos);
            std::string value = line.substr(colonPos + 1);
            
            // Clean up key and value (remove spaces, quotes, commas)
            key.erase(0, key.find_first_not_of(" \t\""));
            key.erase(key.find_last_not_of(" \t\"") + 1);
            value.erase(0, value.find_first_not_of(" \t\""));
            value.erase(value.find_last_not_of(" \t\",") + 1);
            
            // Store parameter values in member variables
            if (key == "x_begin") range[0] = std::stod(value);
            else if (key == "x_end") range[1] = std::stod(value);
            else if (key == "y_begin") {} // Can add y_range if needed
            else if (key == "y_end") {}   // Can add y_range if needed
            else if (key == "x_num") x_num = std::stoi(value);
            else if (key == "y_num") y_num = std::stoi(value);
            else if (key == "thick") thick = std::stod(value);
            else if (key == "k") k = std::stod(value);
            else if (key == "E_boundary") b_val[0] = std::stod(value);
            else if (key == "W_boundary") b_val[1] = std::stod(value);
            else if (key == "N_boundary") b_val[2] = std::stod(value);
            else if (key == "S_boundary") b_val[3] = std::stod(value);
        }
    }
    
    ifs_.close();
    
    // Initialize matrix A, vector b and solution vector T based on grid dimensions
    int total_nodes = x_num * y_num;
    A.resize(total_nodes, std::vector<double>(total_nodes, 0.0));
    b.resize(total_nodes, 0.0);
    T.resize(total_nodes, 0.0);  // Initialize temperature solution vector
    x_size = x_num;
    y_size = y_num;
    
    // Calculate grid spacing
    double dx = (range[1] - range[0]) / (x_num - 1);
    
    // Output loaded parameters for verification
    std::cout << "---------Init---------" << std::endl;
    std::cout << "Parameters loaded from parameter.json:" << std::endl;
    std::cout << "X range: [" << range[0] << ", " << range[1] << "]" << std::endl;
    std::cout << "Grid dimensions: " << x_num << " x " << y_num << " = " << total_nodes << " nodes" << std::endl;
    std::cout << "Matrix A size: " << total_nodes << " x " << total_nodes << std::endl;
    std::cout << "Vector b size: " << total_nodes << std::endl;
    std::cout << "Boundary conditions: E=" << b_val[0] << ", W=" << b_val[1] 
              << ", N=" << b_val[2] << ", S=" << b_val[3] << std::endl;
    std::cout << "Thermal conductivity k = " << k << std::endl;
    std::cout << "Thickness = " << thick << std::endl;
    std::cout << std::endl;
}

void Possion2d::AssembleLinearSystem()
{
    std::cout << "---------Assemble---------" << std::endl;
    
    // Calculate grid spacing for finite volume method
    // For FVM, divide by number of cells, not number of nodes
    // For uniform grid: Δx = Δy = h
    double dx = (range[1] - range[0]) / x_size;
    double dy = (range[1] - range[0]) / y_size; // Assuming square domain

    // Face areas (assuming unit depth in z-direction)
    double A_east = dy * thick;    // East face area
    double A_west = dy * thick;    // West face area
    double A_north = dx * thick;   // North face area
    double A_south = dx * thick;   // South face area
    
    // Distances between cell centers
    double de = dx;        // Distance to east neighbor
    double dw = dx;        // Distance to west neighbor
    double dn = dy;        // Distance to north neighbor
    double ds = dy;        // Distance to south neighbor

    // Thermal conductances including thermal conductivity k
    double G_east = k * A_east / de;   // East face conductance
    double G_west = k * A_west / dw;   // West face conductance
    double G_north = k * A_north / dn;  // North face conductance
    double G_south = k * A_south / ds;  // South face conductance
                

    // Finite Volume Method assembly for steady heat equation
    for(int i=0; i < y_size; ++i)
    {
        for(int j=0; j < x_size; ++j)
        {
            int node = i * x_size + j;  // Current control volume index
            
            // Boundary conditions using Dirichlet BC
            // Bottom boundary (i = 0) - South
            if(i == 0)
            {
                A[node][node] += G_south;
                b[node] += b_val[3] * 2 * G_south; // South boundary value
                A[node][node + x_size] = -G_north;
                if(j !=0 && j != x_size-1)
                {
                    A[node][node - 1] = -G_west;
                    A[node][node + 1] = -G_east;
                }
            }
            // Top boundary (i = y_size - 1) - North
            if(i == y_size - 1)
            {
                A[node][node] += G_north;
                b[node] += b_val[2] * 2 * G_north; // North boundary value
                A[node][node - x_size] = -G_south;
                if(j !=0 && j != x_size-1)
                {
                    A[node][node - 1] = -G_west;
                    A[node][node + 1] = -G_east;
                }
            }
            // Left boundary (j = 0) - West
            if(j == 0)
            {
                A[node][node] += G_west;
                b[node] += b_val[1] * 2 * G_west; // West boundary value
                A[node][node + 1] = -G_east;
                if(i !=0 && i != y_size-1)
                {
                    A[node][node - x_size] = -G_south;
                    A[node][node + x_size] = -G_north;
                }
            }
            // Right boundary (j = x_size - 1) - East
            if(j == x_size - 1)
            {
                A[node][node] += G_east;
                b[node] += b_val[0] * 2 * G_east; // East boundary value
                A[node][node - 1] = -G_west;
                if(i !=0 && i != y_size-1)
                {
                    A[node][node - x_size] = -G_south;
                    A[node][node + x_size] = -G_north;
                }
            }
            // Interior control volumes - FVM for steady heat equation: ∇·(k∇T) + coe*T = 0
            if(i != 0 && i != y_size-1 && j != 0 && j != x_size-1)
            {               
                // Neighbor coefficients (thermal conductances)
                // East neighbor (j+1, i)
                A[node][node + 1] = -G_east;
                
                // West neighbor (j-1, i)
                A[node][node - 1] = -G_west;
                
                // North neighbor (j, i+1)
                A[node][node + x_size] = -G_north;
                
                // South neighbor (j, i-1)
                A[node][node - x_size] = -G_south;
            }

            A[node][node] += G_east + G_west + G_north + G_south;
        }
    }
    std::cout << "FVM matrix assembly completed for steady heat equation." << std::endl;
    std::cout << "Grid spacing: dx = " << dx << ", dy = " << dy << std::endl;
    std::cout << "Control volume volume: " << dx * dy << std::endl;
    for(int i=0; i < y_size+1; ++i)
    {
        for(int j=0; j < x_size+1; ++j)
        {
            std::cout << A[i][j] << "\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
    for(int i=0; i < y_size+1; ++i)
    {
        std::cout << b[i] << "\t";
    }
    std::cout << std::endl;
}

void Possion2d::SolvingLinearSystem_LU()
{
    std::cout << "---------Solving with LU Decomposition---------" << std::endl;
    
    int n = x_size * y_size;
    
    // Create copies of matrix A and vector b for LU decomposition
    std::vector<std::vector<double>> L(n, std::vector<double>(n, 0.0));
    std::vector<std::vector<double>> U(n, std::vector<double>(n, 0.0));
    std::vector<double> b_copy = b;
    
    // Initialize solution vector
    std::fill(T.begin(), T.end(), 0.0);
    
    // LU Decomposition with partial pivoting
    std::cout << "Performing LU decomposition..." << std::endl;
    
    // Initialize U as copy of A
    for(int i = 0; i < n; ++i) {
        for(int j = 0; j < n; ++j) {
            U[i][j] = A[i][j];
        }
    }
    
    // Initialize L as identity matrix
    for(int i = 0; i < n; ++i) {
        L[i][i] = 1.0;
    }
    
    // Doolittle's method for LU decomposition
    for(int i = 0; i < n; ++i) {
        // Upper triangular elements
        for(int j = i; j < n; ++j) {
            double sum = 0.0;
            for(int k = 0; k < i; ++k) {
                sum += L[i][k] * U[k][j];
            }
            U[i][j] = A[i][j] - sum;
        }
        
        // Lower triangular elements
        for(int j = i + 1; j < n; ++j) {
            double sum = 0.0;
            for(int k = 0; k < i; ++k) {
                sum += L[j][k] * U[k][i];
            }
            
            // Check for zero division
            if(std::abs(U[i][i]) < 1e-10) {
                std::cerr << "Warning: Zero pivot detected at row " << i << ". Matrix may be singular." << std::endl;
                L[j][i] = 0.0;
            } else {
                L[j][i] = (A[j][i] - sum) / U[i][i];
            }
        }
    }
    
    std::cout << "LU decomposition completed." << std::endl;
    
    // Forward substitution: Solve L*y = b
    std::vector<double> y(n, 0.0);
    std::cout << "Performing forward substitution..." << std::endl;
    
    for(int i = 0; i < n; ++i) {
        double sum = 0.0;
        for(int j = 0; j < i; ++j) {
            sum += L[i][j] * y[j];
        }
        y[i] = b_copy[i] - sum;
    }
    
    // Back substitution: Solve U*x = y
    std::cout << "Performing back substitution..." << std::endl;
    
    for(int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for(int j = i + 1; j < n; ++j) {
            sum += U[i][j] * T[j];
        }
        
        // Check for zero division
        if(std::abs(U[i][i]) < 1e-10) {
            std::cerr << "Warning: Zero diagonal element in U at row " << i << ". Using fallback solution." << std::endl;
            T[i] = 0.0;
        } else {
            T[i] = (y[i] - sum) / U[i][i];
        }
    }
    
    std::cout << "LU decomposition solution completed." << std::endl;
    
    // Calculate residual for verification
    double max_residual = 0.0;
    for(int i = 0; i < n; ++i) {
        double residual = b_copy[i];
        for(int j = 0; j < n; ++j) {
            residual -= A[i][j] * T[j];
        }
        max_residual = std::max(max_residual, std::abs(residual));
    }
    
    std::cout << "Maximum residual: " << max_residual << std::endl;
    
    if(max_residual < 1e-6) {
        std::cout << "LU decomposition solution converged successfully." << std::endl;
    } else {
        std::cout << "Warning: Large residual detected. Solution may be inaccurate." << std::endl;
    }
}

void Possion2d::SolvingLinearSystem_GS()
{
    std::cout << "---------Solving with Gauss-Seidel---------" << std::endl;
    
    const int max_iterations = 10000;
    const double tolerance = 1e-6;
    int total_nodes = x_size * y_size;
    
    // Initialize solution vector (all nodes start from 0)
    // Boundary conditions are already handled in the matrix assembly
    std::fill(T.begin(), T.end(), 0.0);
    
    // Gauss-Seidel iteration
    double max_error = 0.0;
    for(int iter = 0; iter < max_iterations; ++iter) {
        max_error = 0.0;
        
        // Update all nodes (boundary nodes will converge to boundary values)
        for(int i = 0; i < y_size; ++i) {
            for(int j = 0; j < x_size; ++j) {
                int node = i * x_size + j;
                
                // Store old value for error calculation
                double T_old = T[node];
                
                // Gauss-Seidel update formula: T[i] = (b[i] - Σ(A[i][j]*T[j])) / A[i][i]
                double sum = 0.0;
                
                // Sum contributions from all other nodes
                for(int k = 0; k < x_size * y_size; ++k) {
                    if(k != node) {
                        sum += A[node][k] * T[k];
                    }
                }
                
                // Update temperature
                T[node] = (b[node] - sum) / A[node][node];
                
                // Calculate error (skip boundary nodes as they are fixed)
                if(A[node][node] != 1.0) {  // Boundary nodes have A[i][i] = 1.0
                    double error = std::abs(T[node] - T_old);
                    max_error = std::max(max_error, error);
                }
            }
        }
        
        // Check convergence
        if(max_error < tolerance) {
            std::cout << "Gauss-Seidel converged in " << iter + 1 << " iterations." << std::endl;
            std::cout << "Final residual: " << max_error << std::endl;
            return;
        }
        
        // Print progress every 1000 iterations
        if((iter + 1) % 1000 == 0) {
            std::cout << "Iteration " << iter + 1 << ", max error: " << max_error << std::endl;
        }
    }
    
    std::cout << "Gauss-Seidel did not converge after " << max_iterations << " iterations." << std::endl;
    std::cout << "Final residual: " << max_error << std::endl;
}

void Possion2d::OutputResults()
{
    std::cout << "---------Solution Results---------" << std::endl;
    std::cout << "Temperature distribution:" << std::endl;
    
    // Output temperature field in grid format
    for(int i = 0; i < y_size; ++i) {
        for(int j = 0; j < x_size; ++j) {
            int node = i * x_size + j;
            std::cout << "T[" << i << "][" << j << "] = " << T[node] << "\t";
        }
        std::cout << std::endl;
    }
    
    std::cout << std::endl;
    
    // Calculate and output some statistics
    double T_min = T[0];
    double T_max = T[0];
    double T_avg = 0.0;
    
    for(int i = 0; i < y_size; ++i) {
        for(int j = 0; j < x_size; ++j) {
            int node = i * x_size + j;
            T_min = std::min(T_min, T[node]);
            T_max = std::max(T_max, T[node]);
            T_avg += T[node];
        }
    }
    T_avg /= (x_size * y_size);
    
    std::cout << "Temperature statistics:" << std::endl;
    std::cout << "Minimum temperature: " << T_min << std::endl;
    std::cout << "Maximum temperature: " << T_max << std::endl;
    std::cout << "Average temperature: " << T_avg << std::endl;
    std::cout << std::endl;
    
    // Output corner temperatures for verification
    std::cout << "Corner temperatures:" << std::endl;
    std::cout << "Bottom-left (SW): " << T[0] << std::endl;
    std::cout << "Bottom-right (SE): " << T[x_size - 1] << std::endl;
    std::cout << "Top-left (NW): " << T[(y_size - 1) * x_size] << std::endl;
    std::cout << "Top-right (NE): " << T[y_size * x_size - 1] << std::endl;
}
