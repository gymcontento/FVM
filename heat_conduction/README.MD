<!--
 * @Author: gymcontento herry996341591@gmail.com
 * @Date: 2026-01-03 01:35:49
 * @LastEditors: gymcontento herry996341591@gmail.com
 * @LastEditTime: 2026-01-04 00:01:08
 * @FilePath: \FVM\heat_conduction\README.MD
 * @Description: 
 * 
 * Copyright (c) 2026 by ${git_name_email}, All Rights Reserved. 
-->
# Description
This directory contains code and documentation related to heat conduction simulations. The files included provide information on how to set up, run, and analyze heat conduction models using various numerical methods. Users can find example scripts, input files, and output data to help them understand the principles of heat conduction and apply them to their own projects.
## 26.1.2
添加了热传导的网格序号功能，暂时不考虑程序设计的单一性、模块化等问题，先实现功能为主。

## 26.1.3
数据结构：添加了节点坐标、单元中心点坐标、单元尺寸、温度场、有限体积法系数矩阵等数据结构。
函数：并添加了创建坐标函数、温度场函数、初始化温度函数、有限体积法系数矩阵函数。

## 26.1.4
类：添加了求解器类、材料类
函数：提供了初始化求解器参数和材料参数的函数

## 26.1.5
类：添加了边界条件类
函数：提供了识别单元边界条件的函数
同时修改了类之间成员变量的调用方式（从赋值改为引用）
但未运行通过测试

## 26.1.8
类：添加了组装线性系统类 
函数：提供了组装线性系统矩阵和向量的函数（未处理边界条件）
检查了boundary类bcid矩阵的正确性，找到错误，并进行了修改