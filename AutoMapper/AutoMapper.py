#!/usr/bin/env python3
##############################################################################
# 开发者: Matthew Bone
# 最后更新: 03/08/2021
# 更新者: Matthew Bone
#
# 联系方式:
# Bristol Composites Institute (BCI)
# Department of Aerospace Engineering - University of Bristol
# Queen's Building - University Walk
# Bristol, BS8 1TR
# U.K.
# Email - matthew.bone@bristol.ac.uk
#
# 文件描述:
# 此包装器用于从命令行调用所有 AutoMapper 系列工具。
# Linux 用户可以考虑将此仓库的位置添加到他们的 PATH 中，
# 这样就可以从任何位置调用此包装器。
##############################################################################

import sys
# 检查 Python 版本
pythonVersion = sys.version_info
if pythonVersion[0] == 2: # Python 2 检查 - 由于 shebang 可能不需要
    print('此脚本不兼容 Python 2。请升级到 Python 3.6+')
    sys.exit()

if pythonVersion[0] == 3 and pythonVersion[1] < 6: # 关于早期 Python 3 版本的说明
    print('注意：此代码使用 Python 3.6 添加的插入顺序字典。AutoMapper 可能在早期 Python 3 版本中工作，但结果可能有所不同。')

# 检查是否安装了 natsort
try:
    import natsort
except ModuleNotFoundError:
    print('在 Python 模块中未找到 natsort 包。请继续之前安装 natsort。')
    sys.exit()

import argparse
from AutoMapper.LammpsUnifiedCleaner import file_unifier
from AutoMapper.LammpsToMolecule import lammps_to_molecule
from AutoMapper.MapProcessor import map_processor

# 初始化解析器
parser = argparse.ArgumentParser(description='运行预处理工具用于 LAMMPS 模拟，使用 fix bond/react')

# 命令行参数列表
parser.add_argument('directory', metavar='directory', type=str, nargs=1, help='文件目录，可在 bash 中使用 . 或 $PWD 找到')
parser.add_argument('tool', metavar='tool', type=str, nargs=1, choices=['clean', 'molecule', 'map'], help='要使用的工具名称。可选工具: clean, molecule, map')
parser.add_argument('data_files', metavar='data_files', nargs='+', help='要操作的文件名。如果工具是 "map" 则必须按反应前-反应后顺序提供两个文件。如果是 "clean" 可以是任意顺序的文件列表')
parser.add_argument('--coeff_file', metavar='coeff_file', nargs=1, help='"clean" 工具参数: 要清理的系数文件')
parser.add_argument('--save_name', metavar='save_name', nargs='+', help='"molecule" 和 "map" 工具参数: 新文件的文件名')
parser.add_argument('--ba', metavar='bonding_atoms', nargs='+', help='"map" 工具参数: 将参与创建新键的原子ID，用空格分隔。映射时分子文件之间的原子顺序必须相同')
parser.add_argument('--ebt', metavar='elements_by_type', nargs='+', help='"map" 工具参数: 元素符号序列，顺序与数据文件中指定的类型相同，用空格分隔')
parser.add_argument('--da', metavar='delete_atoms', nargs='+', help='"map" 工具可选参数: 键形成后将被删除的原子ID，用空格分隔')
parser.add_argument('--debug', action='store_true', help='"map" 工具可选参数: 打印调试语句，包含路径搜索和映射处理器的信息')
parser.add_argument('--ca', metavar='create_atoms', nargs='+', help='"map" 工具可选参数: 键形成后将被创建的原子ID，用空格分隔')

# 从解析器获取参数
args = parser.parse_args()
# 从列表中取出必需参数 - 其他参数稍后处理
tool = args.tool[0]
directory = args.directory[0]

# 如果某些工具缺少参数则抛出错误
if tool == 'clean' and (args.coeff_file is None):
    parser.error('"clean" 工具需要 --coeff_file 参数')

if tool == 'molecule' and args.save_name is None:
    parser.error('"molecule" 工具需要 --save_name 参数')

if tool == 'molecule' and len(args.data_files) > 1:
    parser.error('molecule 工具只能接受 1 个 data_file 作为输入')

if tool == 'map' and (len(args.data_files) != 2 or len(args.save_name) != 2):
    parser.error('map 工具需要 2 个 data_files 和 2 个 save_names (用于反应前和反应后分子文件)')

if tool == 'map' and (len(args.ba) < 4 or args.ebt is None):
    parser.error('map 工具需要 --ba (成键原子) 参数指定至少 4 个原子ID，以及 --ebt (按类型的元素) 参数')

# 统一数据文件清理
if tool == "clean":  
    print(f'数据文件列表: {args.data_files}')
    file_unifier(directory, args.coeff_file[0], args.data_files)

# 生成分子数据文件
elif tool == "molecule":
    lammps_to_molecule(directory, args.data_files[0], args.save_name[0])

# 组合分子和映射创建代码
elif tool == 'map':
    map_processor(directory, args.data_files[0], args.data_files[1], args.save_name[0], args.save_name[1], args.ba[:2], args.ba[2:], args.da, args.ebt, args.ca, args.debug)

# 打印消息显示 AutoMapper 已完成
print('AutoMapper 任务完成')