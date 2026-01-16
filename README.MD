<!--
 * @Author: gymcontento herry996341591@gmail.com
 * @Date: 2025-12-28 01:09:24
 * @LastEditors: gymcontento herry996341591@gmail.com
 * @LastEditTime: 2026-01-17 00:55:40
 * @FilePath: \FVM\README.MD
 * @Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
-->
# Description
文件heat_conduction中的代码基于曾导SJTU的从0到1学CFD中的有限体积法代码框架，完成了二维对流扩散方程（传热方程）无压力项的有限体积法求解。

## 2d_possion： 二维泊松方程有限体积法求解
更新可视化部分

## 1d_convection_diffusion: 一维对流扩散方程有限体积法求解
更新可视化部分

## heat_conduction: 二维对流扩散方程（传热方程）无压力项有限体积法求解
### 26.1.2
添加了热传导的网格序号功能，暂时不考虑程序设计的单一性、模块化等问题，先实现功能为主。

### 26.1.3
#### 数据结构：添加了节点坐标、单元中心点坐标、单元尺寸、温度场、有限体积法系数矩阵等数据结构。
#### 函数：并添加了创建坐标函数、温度场函数、初始化温度函数、有限体积法系数矩阵函数。

### 26.1.4
#### 类：添加了求解器类、材料类
#### 函数：提供了初始化求解器参数和材料参数的函数

### 26.1.5
#### 类：添加了边界条件类
#### 函数：提供了识别单元边界条件的函数
同时修改了类之间成员变量的调用方式（从赋值改为引用）但未运行通过测试

### 26.1.8
#### 类：添加了组装线性系统类 
#### 函数：提供了组装线性系统矩阵和向量的函数（未处理边界条件）
检查了boundary类bcid矩阵的正确性，找到错误，并进行了修改

### 26.1.11
#### 类：添加了求解类-雅可比迭代法
#### 函数：提供了雅可比迭代法求解线性系统的函数
#### 测试：完成了热传导问题的完整测试，验证了程序的正确性
对比了https://www.theseus-fe.com/resources/validations/basic-heat-transfer 2D conduction in rectangular plate的解析解，在0.0833处和0.5处的结果一致
#### 添加了.plot文件，程序会生成数据，可以使用gnuplot打开plot文件对比数值解和解析解

### 26.1.15
#### 添加了Convection_Diffusion文件，完成了一维对流扩散方程的有限体积法求解
#### 加入了对流项的迎风格式和中心差分格式
#### 还需对程序在不同Peclect数下的结果进行测试

### 26.1.16
#### 在Convection_Diffusion文件Solver中添加了高斯-赛德尔迭代法
#### 添加了plot文件，方便使用gnuplot进行可视化
#### 添加了不同Peclet数下程序的测试脚本（Windows）
#### 验证了不同Peclet数下程序的结果，见pe_uw_result.pdf，结果和解析解一致，结果体现迎风格式的数值耗散特性，相比于解析解，更加早的出现弯曲