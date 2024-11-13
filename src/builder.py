# builder.py
# 该模块用于生成体系与反应模板文件

from src.filewriter import *
from src.execute import *
from src.readdata import get_map_type_str

def construct_system(path_dict, box_len = 60):
    """构建分子系统。

    Args:
        path_dict: 包含反应物、溶剂及产物的文件名、SMILES 和路径信息的字典。
            - names: 各种分子的文件名
            - num: 各种分子的数量
            - paths: 各种路径信息
        box_len: 模拟盒子的长度。

    """
    # 从 path_dict 中读取所需信息
    data_PATH = path_dict['paths']['data']
    lt_PATH = path_dict['paths']['lt']
    sys_PATH = path_dict['paths']['sys']
    
    # 动态生成 name_list 和 num_list
    name_list = [path_dict['names']['r1'], path_dict['names']['r2'], path_dict['names']['sol']]
    num_list = [path_dict['num']['r1'], path_dict['num']['r2'], path_dict['num']['sol']]
    
    # 执行各个步骤
    exec_ltemplify(data_PATH, lt_PATH, name_list)
    write_sys_lt(lt_PATH, name_list, num_list)
    write_packmol_inp(sys_PATH, path_dict['paths']['mol'], name_list, num_list, box_len)
    exec_packmol(sys_PATH)
    
    

def construct_reaction_map(path_dict):
    """构建反应映射模板。

    Args:
        path_dict: 包含反应物、溶剂及产物的文件名、SMILES 和路径信息的字典。
            - names: 各种分子的文件名
            - paths: 各种路径信息

    """
    # 从 path_dict 中读取所需信息
    lt_PATH = path_dict['paths']['lt']
    sys_PATH = path_dict['paths']['sys']
    data_PATH = path_dict['paths']['data']
    map_PATH = path_dict['paths']['map']
    
    
    
    reaction_index_dicts_list = path_dict['reaction_index_dicts_list']

    # 获取 pre_list 和 post_list 分子名列表
    pre_list = [path_dict['names']['r1'], path_dict['names']['r2']]
    post_lists = [[item, path_dict['names']['byp']] for item in path_dict['names']['p']]
    
    # 构建两个分子为体系的 Lammps DATA 文件
    for index, post_list in enumerate(post_lists):
        # 构建产物 LT 文件
        exec_ltemplify(data_PATH, lt_PATH, post_list)
        write_template_lt(lt_PATH, pre_list, post_list, index)
    
    exec_moltemplate(lt_PATH, sys_PATH, data_PATH, loop = len(post_lists))
    
    # 清理整理力场信息
    exec_AutoMapper_clean(data_PATH, map_PATH)
    
    # 读取清理后的 data 获取元素种类字符串
    atom_ele_type_str = get_map_type_str(data_PATH + 'cleanedsystem.data')
    exec_AutoMapper(map_PATH, reaction_index_dicts_list, atom_ele_type_str)

