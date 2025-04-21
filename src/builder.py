# builder.py
# 该模块用于生成体系与反应模板文件

from src.filewriter import *
from src.execute import *
from src.readdata import get_map_type_str
import logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def construct_system(path_dict, box_len = 60):
    """构建分子系统。

    Args:
        path_dict: 包含反应物、溶剂及产物的文件名、SMILES 和路径信息的字典。
            - file_name: 按顺序存储的分子文件名列表 [r1, r2, sol, p1, p2, ..., byp, p1, p2, ...]
            - type: 各种分子的类型字典
            - num: 各种分子的数量字典
            - paths: 各种路径信息
        box_len: 模拟盒子的长度。

    """
    # 从 path_dict 中读取所需信息
    data_PATH = path_dict['paths']['data']
    lt_PATH = path_dict['paths']['lt']
    sys_PATH = path_dict['paths']['sys']
    
    # 获取反应物和溶剂的文件名和数量
    # 按照file_name列表的顺序，前三个元素分别是r1, r2, sol
    name_list = path_dict['file_name'][:3]  # 获取r1, r2, sol
    num_list = [path_dict['num'][name] for name in name_list]
    
    # 执行各个步骤
    write_sys_lt(lt_PATH, name_list, num_list, box_len)
    write_packmol_inp(sys_PATH, path_dict['paths']['mol'], name_list, num_list, box_len)
    exec_packmol(sys_PATH)
    exec_moltemplate_system(lt_PATH, sys_PATH, data_PATH)

def prepare_reaction_templates(path_dict):
    """准备反应映射模板。
    
    构建反应前后的分子列表和反应映射模板，并将其添加到path_dict中。
    
    Args:
        path_dict: 包含反应物、溶剂及产物的文件名、SMILES 和路径信息的字典。
            - file_name: 按顺序存储的分子文件名列表 [r1, r2, sol, p1, p2, ..., byp, p1, p2, ...]
            - type: 各种分子的类型字典
            - paths: 各种路径信息
            - p_smiles_num: 产物分子的数量
            
    Returns:
        无返回值，但会在path_dict中添加以下键：
        - pre_lists: 反应前分子列表
        - post_lists: 反应后分子列表
        - map_templates: 反应映射模板字典
    """
    # 获取 pre_list 和 post_list 分子名列表    
    p_count = path_dict['p_smiles_num']
    # 产物在file_name中的位置是从索引3开始，到3+p_count结束
    p_names = path_dict['file_name'][3:3+p_count]

    # 根据file_name列表的顺序，产物后一个元素是byp
    byp_name = path_dict['file_name'][3+p_count]

    # 根据产物和二次产物进行循环
    # 输出首次反应的反应物和产物信息
    # 输出首次反应信息
    logging.info("首次反应信息:")
    pre_lists = []
    post_lists = []
    
    for p_name in p_names:
        name_list = p_name.split('_')
        # 解析命名获取各反应物和产物
        r1_name = name_list[0] + '_' + name_list[1]  # r1_1
        r2_name = name_list[2] + '_' + name_list[3]  # r2_1
        logging.info(f"Map：反应物 [{r1_name}, {r2_name}], 产物 [{p_name}, {byp_name}]")
        pre_lists.append([r1_name, r2_name])
        post_lists.append([p_name, byp_name])
    
    # 二次反应产物
    if 'p_smiles_num_2nd' in path_dict:
        p_count_2nd = path_dict['p_smiles_num_2nd']
        second_reaction_p_names = path_dict['file_name'][4+p_count:4+p_count+p_count_2nd]
        
        # 输出二次反应的反应物和产物信息
        logging.info("二次反应信息:")
        for p_name in second_reaction_p_names:
            # 例如r1_1_r2_1_r2_1_2nd_1
            # 前r1_1是r1, r2_1是r2, r1_1_r2_1是一次次产物作为反应物1, r2_1是二次产物作为反应物2
            name_list = p_name.split('_')
            # 解析命名获取各反应物和产物
            # r1_1
            r1_name = name_list[0] + '_' + name_list[1]
            # r2_1
            r2_name = name_list[2] + '_' + name_list[3]
            # r1_1_r2_1
            first_product = r1_name + '_' + r2_name + '_' + name_list[4]
            logging.info(f"Map：二级反应物 [{first_product}, {r2_name}], 产物 [{p_name}, {byp_name}]")
            pre_lists.append([first_product, r2_name])
            post_lists.append([p_name, byp_name])
    
    # 添加反应映射模板字典
    map_templates = {}
    for i in range(len(pre_lists)):
        # 构建映射模板：[反应物1，反应物2，产物，副产物]
        map_templates[f'map_{i+1}'] = [pre_lists[i][0], pre_lists[i][1],
                                    post_lists[i][0], post_lists[i][1]]
    # 将映射模板添加到path_dict
    path_dict['map_templates'] = map_templates
    logging.info("反应映射模板准备完成")
    return path_dict

def construct_reaction_map(path_dict):
    """构建反应映射模板。

    Args:
        path_dict: 包含反应物、溶剂及产物的文件名、SMILES 和路径信息的字典。
            - file_name: 按顺序存储的分子文件名列表 [r1, r2, sol, p1, p2, ..., byp, p1, p2, ...]
            - type: 各种分子的类型字典
            - paths: 各种路径信息
            - p_smiles_num: 产物分子的数量

    """
    # 从 path_dict 中读取所需信息
    lt_PATH = path_dict['paths']['lt']
    data_PATH = path_dict['paths']['data']
    map_PATH = path_dict['paths']['map']
    reaction_index_dicts_list = path_dict['reaction_index_dicts_list']
    
    # 获取映射模板字典
    map_templates = path_dict['map_templates']
    
    for map_name in map_templates.keys():
        # 构建反应前后体系 LT 文件
        pre_list = map_templates[map_name][:2]
        pre_map_name = f"pre_{map_name}"
        write_react_lt(lt_PATH, pre_map_name, pre_list)
        post_list = map_templates[map_name][2:]
        post_map_name = f"post_{map_name}"
        write_react_lt(lt_PATH, post_map_name, post_list)
        exec_moltemplate_reaction(lt_PATH, data_PATH, map_name)
    
    # 清理整理力场信息
    exec_AutoMapper_clean(data_PATH, map_PATH, map_templates.keys())
    
    # 读取清理后的 data 获取元素种类字符串
    atom_ele_type_str = get_map_type_str(data_PATH + 'cleanedsystem.data')
    exec_AutoMapper(map_PATH, map_templates.keys(), reaction_index_dicts_list, atom_ele_type_str)
    path_dict['elem_type_str'] = atom_ele_type_str
    return path_dict