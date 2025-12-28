# 高斯塞德尔求解结果可视化

## 项目概述

本项目实现了2D泊松方程的高斯塞德尔求解器，并提供Python可视化脚本来展示求解结果。

## 文件说明

### C++源文件
- `2d_possion.h` - 主要类定义
- `2d_possion.cpp` - 有限体积法实现和高斯塞德尔求解器
- `main.cpp` - 主程序入口
- `parameter.json` - 参数配置文件

### Python可视化文件
- `visualize_gauss_seidel.py` - 可视化脚本
- `gauss_seidel_results.csv` - 求解结果数据文件

### 生成的图形文件
- `gauss_seidel_contour_3d.png` - 等高线图和3D表面图
- `gauss_seidel_heatmap.png` - 热力图

## 使用方法

### 1. 编译并运行C++程序
```bash
cd 2d_possion
mkdir -p build && cd build
cmake ..
cmake --build .
cd ..
./build/Debug/2d_possion.exe
```

### 2. 运行Python可视化脚本
```bash
cd 2d_possion
python visualize_gauss_seidel.py
```

## 参数配置

在 `parameter.json` 中可以配置以下参数：
- `x_begin`, `x_end` - X方向域范围
- `y_begin`, `y_end` - Y方向域范围  
- `x_num`, `y_num` - 网格点数量
- `k` - 热导率
- `thick` - 厚度
- `E_boundary`, `W_boundary`, `N_boundary`, `S_boundary` - 东、西、北、南边界条件

## 求解结果

### 当前配置结果
- **网格大小**: 4×4 = 16个节点
- **域范围**: [0, 4] × [0, 4]
- **边界条件**: 
  - 东边界: 200°C
  - 西边界: 100°C  
  - 北边界: 150°C
  - 南边界: 250°C

### 温度统计
- **最小温度**: 123.845°C (位置: 0.00, 2.67)
- **最大温度**: 226.155°C (位置: 2.67, 0.00)
- **平均温度**: 175.000°C
- **标准差**: 30.946°C

### 收敛信息
- **高斯塞德尔迭代次数**: 32次
- **最终残差**: 7.064e-07
- **收敛容差**: 1e-6

## 可视化特性

Python脚本提供三种可视化方式：

1. **等高线图** - 显示温度等值线和数值标注
2. **3D表面图** - 立体展示温度分布
3. **热力图** - 使用颜色映射显示温度梯度

所有图形都使用 `coolwarm` 颜色映射，蓝色表示低温，红色表示高温。

## 依赖项

### C++
- CMake
- C++11或更高版本的编译器

### Python
- numpy
- pandas
- matplotlib

安装Python依赖：
```bash
pip install numpy pandas matplotlib
```

## 技术细节

### 数值方法
- **离散化**: 有限体积法 (FVM)
- **求解器**: 高斯塞德尔迭代法
- **边界条件**: Dirichlet边界条件
- **收敛准则**: 最大残差 < 1e-6

### 可视化技术
- **数据格式**: CSV格式，包含坐标和温度值
- **网格重建**: 基于数据点创建规则网格
- **颜色映射**: matplotlib的coolwarm映射
- **图形输出**: 高分辨率PNG格式 (300 DPI)

## 注意事项

1. 运行C++程序前确保 `parameter.json` 文件存在
2. 运行Python脚本前确保已生成 `gauss_seidel_results.csv` 文件
3. 图形文件会自动保存到当前目录
4. Python脚本会显示交互式图形窗口
