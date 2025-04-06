# AutoMapper 项目说明文档

## 项目概述

AutoMapper 是一套工具包，专为自动化创建 LAMMPS（Large-scale Atomic/Molecular Massively Parallel Simulator）中 bond/react 功能所需的文件而设计。该工具主要处理分子动力学模拟中的化学反应映射问题，能够处理创建原子并支持新的类型 ID 字符串。

## 项目结构

项目由多个 Python 模块组成，每个模块负责特定的功能：

- `AutoMapper.py`：主程序入口，提供命令行接口
- `LammpsUnifiedCleaner.py`：统一不同文件中的类型定义
- `LammpsToMolecule.py`：将 LAMMPS 输入文件转换为分子文件格式
- `MapProcessor.py`：创建反应前后的分子文件和映射文件
- `LammpsSearchFuncs.py`：提供搜索 LAMMPS 数据的函数
- `LammpsTreatmentFuncs.py`：处理 LAMMPS 数据的通用函数
- `AtomObjectBuilder.py`：构建原子对象
- `QueueFuncs.py`：队列处理函数
- `PathSearch.py`：路径搜索算法

## 主要功能

AutoMapper 提供三个主要工具：

1. **clean**：统一两个或多个文件之间的类型（如原子、键、角度等），并移除未使用的系数
2. **molecule**：将 LAMMPS 输入文件转换为 LAMMPS 分子文件格式
3. **map**：创建反应前后的分子文件和映射文件，满足使用 `bond/react` 的全部要求

## 工作流程

### 1. 清理和统一数据（clean）

`LammpsUnifiedCleaner.py` 中的 `file_unifier` 函数负责统一不同文件中的类型定义：

1. 加载并清理数据文件
2. 提取头部信息
3. 初始化 `Data` 类对象
4. 合并所有类型（原子、键、角度等）
5. 更新各部分的类型编号
6. 更新头部信息
7. 保存清理后的数据文件
8. 处理设置文件中的系数
9. 保存清理后的系数文件

### 2. 转换为分子文件（molecule）

`LammpsToMolecule.py` 负责将 LAMMPS 数据文件转换为分子文件格式：

1. 加载清理后的数据文件
2. 提取相关原子和它们的连接信息
3. 重新编号原子
4. 生成分子文件格式的输出

### 3. 创建映射文件（map）

`MapProcessor.py` 负责创建反应前后的映射关系：

1. 加载反应前后的数据
2. 构建原子对象，包括邻居关系
3. 使用广度优先搜索（BFS）算法寻找路径
4. 确定需要映射的原子
5. 处理需要删除或创建的原子
6. 生成映射文件

## 核心类和函数

### Data 类 (`LammpsUnifiedCleaner.py`)

管理 LAMMPS 数据文件的内容，提供方法获取和修改各种类型：

- `__init__`：初始化数据结构
- `get_atom_types`、`get_bond_types` 等：获取各种类型
- `change_mass_types`、`change_atom_types` 等：更改类型编号
- `change_header`：更新头部信息

### 队列处理 (`QueueFuncs.py`)

提供队列操作功能，用于广度优先搜索：

- `Queue` 类：基本队列实现
- `run_queue`：运行队列处理映射关系

### 路径搜索 (`MapProcessor.py`)

- `bfs`：广度优先搜索算法，用于寻找原子间的路径
- `get_byproducts`：确定副产物
- `keep_all_neighbours`：保留所有邻居原子

### 数据处理函数 (`LammpsTreatmentFuncs.py`)

- `clean_data`：清理数据，移除空行和注释
- `refine_data`：根据原子 ID 集合精炼数据
- `save_text_file`：保存处理后的数据到文件

### 映射处理 (`MapProcessor.py`)

- `output_map`：生成映射文件输出
- `get_byproducts`：识别副产物

## 使用示例

### 清理和统一数据

```bash
AutoMapper.py . clean pre-reaction.data post-reaction.data --coeff_file system.in.settings

## 映射处理过程详解

`MapProcessor.py` 中的映射处理流程分为以下几个关键步骤：

1. **初始分子创建**：
   - 使用 `lammps_to_molecule` 函数分别创建反应前和反应后的分子文件
   - 处理删除原子列表，确保前后原子数量匹配

2. **初始映射创建**：
   - 调用 `map_from_path` 函数生成初始原子映射关系
   - 使用广度优先搜索(BFS)算法寻找原子间的连接路径

3. **部分结构处理**：
   - 检测成键原子是否在环结构中(`is_cyclic`)
   - 判断是否为开环反应(`is_ring_opening`)
   - 保留成键原子周围4个键范围内的原子(`keep_all_neighbours`)

4. **边缘原子处理**：
   - 查找初始边缘原子(`find_edge_atoms`)
   - 验证边缘原子是否需要扩展(`verify_edge_atoms`)
   - 必要时扩展边缘区域(`extend_edge_atoms`)

5. **副产物处理**：
   - 识别并处理副产物原子(`get_byproducts`)

6. **最终输出**：
   - 生成映射文件(`output_map`)
   - 处理部分结构的重新编号

## 边缘原子判断机制

边缘原子的判断和处理是映射过程中的关键环节，具体逻辑如下：

1. **边缘原子定义**：
   - 边缘原子是指位于部分结构边界上的原子
   - 氢原子(H)不会被识别为边缘原子
   - 判断标准：如果一个原子的任一直接邻居(1st neighbour)不在部分原子集中，则该原子被标记为边缘原子

2. **边缘验证流程**：
   - 检查边缘原子及其邻居的原子类型是否在反应前后发生变化
   - 根据类型变化距离边缘的远近，决定需要扩展的范围：
     - 边缘原子本身类型变化：扩展3层邻居
     - 第一邻居类型变化：扩展2层邻居
     - 第二邻居类型变化：扩展1层邻居

3. **边缘扩展处理**：
   - 将需要扩展的邻居原子加入部分原子集
   - 更新映射关系以包含新加入的原子
   - 确保扩展后的边缘区域足够大，避免类型变化区域过于靠近边界

4. **特殊处理**：
   - 对于开环反应，会保留整个环结构及其邻居
   - 副产物原子会被自动包含在部分结构中
   - 删除原子和创建原子会被特殊处理，不受边缘判断影响
