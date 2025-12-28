#!/usr/bin/env python3
"""
高斯塞德尔求解结果可视化脚本
Visualize Gauss-Seidel solver results for 2D Poisson equation
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.ticker as ticker

def load_data(filename):
    """
    从CSV文件加载温度数据
    Load temperature data from CSV file
    """
    # 读取数据，跳过以#开头的注释行
    data = pd.read_csv(filename, comment='#', header=None, 
                      names=['x', 'y', 'temperature'])
    
    return data

def create_grid(data):
    """
    创建规则的网格用于可视化
    Create regular grid for visualization
    """
    # 获取唯一的x和y坐标
    x_unique = np.sort(data['x'].unique())
    y_unique = np.sort(data['y'].unique())
    
    # 创建网格
    X, Y = np.meshgrid(x_unique, y_unique)
    
    # 将温度数据重塑为网格形式
    temp_grid = np.zeros_like(X)
    for i, y in enumerate(y_unique):
        for j, x in enumerate(x_unique):
            mask = (data['x'] == x) & (data['y'] == y)
            temp_grid[i, j] = data[mask]['temperature'].values[0]
    
    return X, Y, temp_grid, x_unique, y_unique

def plot_contour(X, Y, temp_grid, x_unique, y_unique):
    """
    绘制等高线图
    Plot contour map
    """
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))
    
    # 使用coolwarm颜色映射
    cmap = 'coolwarm'
    
    # 等高线图
    contour = ax1.contourf(X, Y, temp_grid, levels=20, cmap=cmap)
    ax1.contour(X, Y, temp_grid, levels=20, colors='black', alpha=0.3, linewidths=0.5)
    
    # 添加颜色条
    cbar1 = plt.colorbar(contour, ax=ax1)
    cbar1.set_label('温度 (°C)', fontsize=12)
    
    # 设置标题和标签
    ax1.set_title('高斯塞德尔求解结果 - 等高线图', fontsize=14, fontweight='bold')
    ax1.set_xlabel('X 坐标', fontsize=12)
    ax1.set_ylabel('Y 坐标', fontsize=12)
    
    # 添加网格点
    ax1.plot(X, Y, 'ko', markersize=3, alpha=0.5)
    
    # 在每个网格点显示温度值
    for i, y in enumerate(y_unique):
        for j, x in enumerate(x_unique):
            ax1.text(x, y, f'{temp_grid[i,j]:.1f}', 
                    ha='center', va='center', fontsize=8, 
                    bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.7))
    
    # 3D表面图
    ax2 = fig.add_subplot(122, projection='3d')
    surf = ax2.plot_surface(X, Y, temp_grid, cmap=cmap, alpha=0.9)
    
    # 设置3D图标题和标签
    ax2.set_title('高斯塞德尔求解结果 - 3D表面图', fontsize=14, fontweight='bold')
    ax2.set_xlabel('X 坐标', fontsize=12)
    ax2.set_ylabel('Y 坐标', fontsize=12)
    ax2.set_zlabel('温度 (°C)', fontsize=12)
    
    # 添加颜色条
    cbar2 = plt.colorbar(surf, ax=ax2, shrink=0.5, aspect=5)
    cbar2.set_label('温度 (°C)', fontsize=12)
    
    plt.tight_layout()
    return fig

def plot_heatmap(X, Y, temp_grid):
    """
    绘制热力图
    Plot heatmap
    """
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # 使用coolwarm颜色映射
    cmap = 'coolwarm'
    
    # 热力图
    im = ax.imshow(temp_grid, extent=[X.min(), X.max(), Y.min(), Y.max()], 
                   origin='lower', cmap=cmap, aspect='auto')
    
    # 添加颜色条
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('温度 (°C)', fontsize=12)
    
    # 设置标题和标签
    ax.set_title('高斯塞德尔求解结果 - 热力图', fontsize=14, fontweight='bold')
    ax.set_xlabel('X 坐标', fontsize=12)
    ax.set_ylabel('Y 坐标', fontsize=12)
    
    # 添加网格
    ax.grid(True, alpha=0.3)
    
    # 在每个单元格中心显示温度值
    for i in range(temp_grid.shape[0]):
        for j in range(temp_grid.shape[1]):
            x_pos = X[i, j]
            y_pos = Y[i, j]
            ax.text(x_pos, y_pos, f'{temp_grid[i, j]:.1f}', 
                   ha='center', va='center', fontsize=10, fontweight='bold',
                   bbox=dict(boxstyle='round,pad=0.3', facecolor='white', alpha=0.8))
    
    plt.tight_layout()
    return fig

def print_statistics(data, temp_grid, X, Y):
    """
    打印统计信息
    Print statistics
    """
    print("\n" + "="*50)
    print("高斯塞德尔求解结果统计信息")
    print("="*50)
    print(f"最小温度: {temp_grid.min():.3f} °C")
    print(f"最大温度: {temp_grid.max():.3f} °C")
    print(f"平均温度: {temp_grid.mean():.3f} °C")
    print(f"标准差: {temp_grid.std():.3f} °C")
    
    # 找到最高温和最低温的位置
    min_idx = np.unravel_index(np.argmin(temp_grid), temp_grid.shape)
    max_idx = np.unravel_index(np.argmax(temp_grid), temp_grid.shape)
    
    print(f"\n最低温度位置: ({X[min_idx]:.2f}, {Y[min_idx]:.2f})")
    print(f"最高温度位置: ({X[max_idx]:.2f}, {Y[max_idx]:.2f})")
    
    # 边界条件验证
    print(f"\n边界条件验证:")
    print(f"东边界 (x=4): 温度范围 [{temp_grid[:, -1].min():.1f}, {temp_grid[:, -1].max():.1f}] °C")
    print(f"西边界 (x=0): 温度范围 [{temp_grid[:, 0].min():.1f}, {temp_grid[:, 0].max():.1f}] °C")
    print(f"北边界 (y=4): 温度范围 [{temp_grid[-1, :].min():.1f}, {temp_grid[-1, :].max():.1f}] °C")
    print(f"南边界 (y=0): 温度范围 [{temp_grid[0, :].min():.1f}, {temp_grid[0, :].max():.1f}] °C")

def main():
    """
    主函数
    Main function
    """
    # 设置中文字体（如果支持）
    plt.rcParams['font.sans-serif'] = ['SimHei', 'DejaVu Sans']
    plt.rcParams['axes.unicode_minus'] = False
    
    # 加载数据
    filename = 'gauss_seidel_results.csv'
    print(f"正在加载数据文件: {filename}")
    
    try:
        data = load_data(filename)
        print(f"成功加载 {len(data)} 个数据点")
    except FileNotFoundError:
        print(f"错误: 找不到文件 {filename}")
        print("请先运行C++程序生成数据文件")
        return
    
    # 创建网格
    X, Y, temp_grid, x_unique, y_unique = create_grid(data)
    
    # 打印统计信息
    print_statistics(data, temp_grid, X, Y)
    
    # 绘制图形
    print("\n正在生成可视化图形...")
    
    # 等高线图和3D表面图
    fig1 = plot_contour(X, Y, temp_grid, x_unique, y_unique)
    fig1.savefig('gauss_seidel_contour_3d.png', dpi=300, bbox_inches='tight')
    print("保存图形: gauss_seidel_contour_3d.png")
    
    # 热力图
    fig2 = plot_heatmap(X, Y, temp_grid)
    fig2.savefig('gauss_seidel_heatmap.png', dpi=300, bbox_inches='tight')
    print("保存图形: gauss_seidel_heatmap.png")
    
    # 显示图形
    plt.show()
    
    print("\n可视化完成！")

if __name__ == "__main__":
    main()
