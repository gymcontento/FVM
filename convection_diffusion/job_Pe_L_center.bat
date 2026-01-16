@echo off
setlocal enabledelayedexpansion
REM 
REM @Author: gymcontento herry996341591@gmail.com
REM @Date: 2026-01-16 10:39:19
REM @LastEditors: gymcontento herry996341591@gmail.com
REM @LastEditTime: 2026-01-16 10:54:30
REM @FilePath: \convection_diffusion\job_Pe_L_center.bat
REM @Description: Run C++ simulation with varying Pe_L parameters
REM 
REM Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 

REM Remove old results
del *.res post* *.dat 2>nul

REM Build the project
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
cd ..

REM Loop for negative velocity (uf = -1.0)
set initialUF=-1.0
for %%c in (0.1 0.1111111111111111 0.125 0.14285714285714285 0.16666666666666666 0.2 0.25 0.3333333333333333 0.5 1.0) do (
    set conductivity=%%c
    echo Running simulation with uf = !initialUF!, conductivity = !conductivity!
    
    REM Modify main.cpp parameters using PowerShell
    copy src\main.cpp src\main_ori.cpp >nul
    powershell -Command "(Get-Content src\main.cpp) -replace 'float initialUF = 1\.0f; //uf', 'float initialUF = !initialUF!f; //uf' | Set-Content src\main.cpp"
    powershell -Command "(Get-Content src\main.cpp) -replace 'float conductivity\{1000000\.0f\};', 'float conductivity{!conductivity!f};' | Set-Content src\main.cpp"
    
    REM Rebuild with new parameters
    cd build
    cmake --build . --config Release
    cd ..
    
    REM Run the simulation
    .\build\Release\demo.exe
    
    REM Restore original main.cpp
    move /Y src\main_ori.cpp src\main.cpp >nul
)

REM Loop for positive velocity (uf = 1.0)
set uf=1.0
for %%c in (1.0 0.5 0.3333333333333333 0.25 0.2 0.16666666666666666 0.14285714285714285 0.125 0.1111111111111111 0.1) do (
    set con=%%c
    echo Running simulation with uf = !uf!, conductivity = !con!
    
    REM Modify main.cpp parameters using PowerShell
    copy src\main.cpp src\main_ori.cpp >nul
    powershell -Command "(Get-Content src\main.cpp) -replace 'float initialUF = 1\.0f; //uf', 'float initialUF = !uf!f; //uf' | Set-Content src\main.cpp"
    powershell -Command "(Get-Content src\main.cpp) -replace 'float conductivity\{1000000\.0f\};', 'float conductivity{!con!f};' | Set-Content src\main.cpp"
    
    REM Rebuild with new parameters
    cd build
    cmake --build . --config Release
    cd ..
    
    REM Run the simulation
    .\build\Release\demo.exe
    
    REM Restore original main.cpp
    move /Y src\main_ori.cpp src\main.cpp >nul
)

echo All simulations completed!
