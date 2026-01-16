@echo off
setlocal enabledelayedexpansion
REM 
REM @Author: gymcontento herry996341591@gmail.com
REM @Date: 2026-01-16 12:47:54
REM @LastEditors: gymcontento herry996341591@gmail.com
REM @LastEditTime: 2026-01-16 12:50:30
REM @FilePath: \convection_diffusion\job_Pe_L.bat
REM @Description: Run C++ simulation with varying conductivity (Pe_L)
REM 
REM Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 

REM Remove old results
del *.res post* temp_x_*.dat analytical_temp_x_*.dat 2>nul

REM Build the project initially
cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
cmake --build . --config Release
cd ..

REM Loop through different conductivity values
for %%c in (100000.0 1.0 0.5 0.25 0.1 0.01) do (
    set conductivity=%%c
    
    REM Calculate Pe_L using PowerShell
    for /f %%p in ('powershell -Command "[math]::Round(1.0 / !conductivity!, 10)"') do set Pe_L=%%p
    
    echo Running simulation with conductivity = !conductivity!, Pe_L = !Pe_L!
    
    REM Modify main.cpp parameters using PowerShell
    copy src\main.cpp src\main_ori.cpp >nul
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

echo All simulations completed!
