#include "2d_cd.h"
#include <iostream>
#include <fstream>

void Convection_Diffusion2d::Init()
{
    std::ifstream ifs_("parameters.json");
    if(!ifs_.is_open())
    {
        std::cerr << "Error: Cannot open parameter.json file!" << std::endl;
        return; 
    }

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
            if (key == "x_begin") range_x[0] = std::stod(value);
            else if (key == "x_end") range_x[1] = std::stod(value);
            else if (key == "y_begin") range_y[0] = std::stod(value);// Can add y_range if needed
            else if (key == "y_end") range_y[1 ] = std::stod(value);   // Can add y_range if needed
            else if (key == "x_num") x_num = std::stoi(value);
            else if (key == "y_num") y_num = std::stoi(value);
            else if (key == "k") k = std::stod(value);
            else if (key == "velocity") velocity = std::stod(value);
            else if (key == "density") density = std::stod(value);
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
    // Output loaded parameters for verification
    std::cout << "---------Init---------" << std::endl;
    std::cout << "Parameters loaded from parameter.json:" << std::endl;
    std::cout << "X range: [" << range_x[0] << ", " << range_x[1] << "]" << std::endl;
    std::cout << "Y range: [" << range_y[0] << ", " << range_y[1] << "]" << std::endl;
    std::cout << "Grid dimensions: " << x_num << " x " << y_num << " = " << total_nodes << " nodes" << std::endl;
    std::cout << "Matrix A size: " << total_nodes << " x " << total_nodes << std::endl;
    std::cout << "Vector b size: " << total_nodes << std::endl;
    std::cout << "Boundary conditions: E=" << b_val[0] << ", W=" << b_val[1] 
              << ", N=" << b_val[2] << ", S=" << b_val[3] << std::endl;
    std::cout << "Thermal conductivity k = " << k << std::endl;
    std::cout << "Velocity: " << velocity << std::endl;
    std::cout << "Density = " << density << std::endl;
    std::cout << std::endl;
}

void Convection_Diffusion2d::AssembleLinearSystem()
{
    double dx = (range_x[1] - range_x[0]) / x_size;
    double dy = (range_y[1] - range_y[0]) / y_size;

    double coe_fe = density * velocity; 
    double coe_fw = density * velocity;
    double coe_fs = density * velocity;
    double coe_fn = density * velocity;

    double coe_de = k / dx;
    double coe_dw = k / dx;
    double coe_dn = k / dy;
    double coe_ds = k / dy;

    for(int i=0; i < y_size; ++i)
    {
        for(int j=0; j < x_size; ++j)
        {
            int node = i * x_size + j;
            //South Boundary
            if(i == 0)
            {
                A[node][node] += coe_fs / 2.0 + coe_ds;
                b[node] += b_val[3] * (coe_fs + 2.0 * coe_ds);
                A[node][node + x_size] = -(coe_dn - coe_fn / 2.0);
                if(j !=0 && j != x_size - 1)
                {
                    A[node][node - 1] = -(coe_dw - coe_fw / 2.0);
                    A[node][node + 1] = -(coe_de - coe_fe / 2.0); 
                } 
            }
            // Top boundary (i = y_size - 1) - North
            if(i == y_size - 1)
            {
                A[node][node] += -coe_fn / 2.0 + coe_dn;
                b[node] += b_val[2] * (-coe_fn + 2.0 * coe_dn); // North boundary value
                A[node][node - x_size] = -(coe_ds - coe_fs / 2.0);
                if(j !=0 && j != x_size-1)
                {
                    A[node][node - 1] = -(coe_dw - coe_fw / 2.0);
                    A[node][node + 1] = -(coe_de - coe_fe / 2.0);
                }
            }
            // Left boundary (j = 0) - West
            if(j == 0)
            {
                A[node][node] += coe_fw / 2.0 + coe_dw;
                b[node] += b_val[1] * (coe_fw + 2.0 * coe_dw); // West boundary value
                A[node][node + 1] = -(coe_dw - coe_fw / 2.0);
                if(i !=0 && i != y_size-1)
                {
                    A[node][node - x_size] = -(coe_ds - coe_fs / 2.0);
                    A[node][node + x_size] = -(coe_dn - coe_fn / 2.0);
                }
            }
            // Right boundary (j = x_size - 1) - East
            if(j == x_size - 1)
            {
                A[node][node] += -coe_fe / 2.0 + coe_de;
                b[node] += b_val[0] * (-coe_fe + 2.0 * coe_de); // East boundary value
                A[node][node - 1] = -(coe_de - coe_fe / 2.0);
                if(i !=0 && i != y_size-1)
                {
                    A[node][node - x_size] = -(coe_ds - coe_fs / 2.0);
                    A[node][node + x_size] = -(coe_dn - coe_fn / 2.0);
                }
            }
            // Interior control volumes - FVM for steady heat equation: ∇·(k∇T) + coe*T = 0
            if(i != 0 && i != y_size-1 && j != 0 && j != x_size-1)
            {               
                // Neighbor coefficients (thermal conductances)
                // East neighbor (j+1, i)
                A[node][node + 1] = -(coe_de - coe_fe / 2.0);
                
                // West neighbor (j-1, i)
                A[node][node - 1] = -(coe_dw - coe_fw / 2.0);
                
                // North neighbor (j, i+1)
                A[node][node + x_size] = -(coe_dn - coe_fn / 2.0);
                
                // South neighbor (j, i-1)
                A[node][node - x_size] = -(coe_ds - coe_fs / 2.0);
            }

            A[node][node] += coe_de + coe_de + coe_de + coe_de
                 + (coe_fe - coe_fw + coe_fn - coe_fs) / 2.0;
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
