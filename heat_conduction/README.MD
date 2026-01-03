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