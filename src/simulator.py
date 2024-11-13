# simulator.py
# 该模块为模拟模块，包括初始化路径、创建写入文件路径、执行模拟等

import os
import logging
import ast
from pysimm.system import System
from pysimm.lmps import Simulation

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def simulation_init(db, id, tmp_path):
    """初始化模拟，包括从数据库读取反应物和溶剂的 SMILES，创建文件路径。

    Args:
        db: 数据库对象，用于获取反应物和溶剂的 SMILES。
        id: 反应的唯一标识符。
        data_path: 数据文件的输出路径。
        tmp_path: 临时文件路径，用于存放模拟结果。

    Returns:
        path_dict: 包含反应物、溶剂及产物的文件名和 SMILES 的字典，包括路径信息。
    """
    
    # 从数据库中读取反应物和溶剂的 SMILES 和其他信息
    try:
        reactant1_smiles = db.get_value_by_id(id, 'reactant1_smiles')
        reactant1_id = db.get_value_by_id(id, 'reactant1_key')
        r1_num = db.get_value_by_id(id, 'r1_num')

        reactant2_smiles = db.get_value_by_id(id, 'reactant2_smiles')
        reactant2_id = db.get_value_by_id(id, 'reactant2_key')
        r2_num = db.get_value_by_id(id, 'r2_num')

        solvent_smiles = db.get_value_by_id(id, 'solvent_smiles')
        solvent_id = db.get_value_by_id(id, 'solvent_key')
        sol_num = db.get_value_by_id(id, 'sol_num')

        product_smiles = db.get_value_by_id(id, 'product_smiles')
        byproduct_smiles = db.get_value_by_id(id, 'byproduct_smiles')
        
        reaction_index_dicts_str = db.get_value_by_id(id, 'reaction_index_dicts')
        reaction_index_dicts_list = ast.literal_eval(reaction_index_dicts_str)

    except KeyError as e:
        logging.error(f"未找到 ID {id} 的数据: {e}")
        raise

    product_smiles_list = product_smiles.split(';')

    # 定义文件名变量
    r1_file_name = f'r1_{reactant1_id}'
    r2_file_name = f'r2_{reactant2_id}'
    sol_file_name = f'sol_{solvent_id}'
    p_file_name_list = [f'{r1_file_name}_{r2_file_name}_{i}' for i, _ in enumerate(product_smiles_list)]
    byp_file_name = 'byp'

    # 返回字典
    path_dict = {
        'names': {
            'r1': r1_file_name,
            'r2': r2_file_name,
            'sol': sol_file_name,
            'p': p_file_name_list,
            'byp': byp_file_name
        },
        'smiles': {
            'r1': reactant1_smiles,
            'r2': reactant2_smiles,
            'sol': solvent_smiles,
            'p': product_smiles_list,
            'byp': byproduct_smiles
        },
        'num': {
            'r1': r1_num,
            'r2': r2_num,
            'sol': sol_num
        },
        'reaction_index_dicts_list': reaction_index_dicts_list
    }

    # 创建反应系统文件路径
    base_PATH = os.getcwd()
    reaction_path = os.path.join(base_PATH, tmp_path, str(id))
    # 确保路径为绝对路径
    reaction_path = os.path.abspath(reaction_path)
    os.makedirs(reaction_path, exist_ok=True)

    # 创建各个子目录路径
    paths = {
        'lt': os.path.join(reaction_path, 'lt') + '/',
        'mol': os.path.join(reaction_path, 'mol') + '/',
        'sys': os.path.join(reaction_path, 'sys') + '/',
        'map': os.path.join(reaction_path, 'map') + '/',
        'data': os.path.join(reaction_path, 'data') + '/',
        'run': os.path.join(reaction_path, 'run') + '/'
    }

    # 创建所有子目录
    for path in paths.values():
        os.makedirs(path, exist_ok=True)
        logging.info(f"创建目录: {path}")

    # 将路径信息添加到 path_dict 中
    path_dict['paths'] = paths
    return path_dict

# 继承 pysimm 包 System 类
class SystemModule(System):
    def __init__(self):
        super().__init__()




# 继承 pysimm 包 Simulation 类
class SimulationModule(Simulation):
    def __init__(self):
        super().__init__()
