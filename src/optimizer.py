# optimizer.py
# 该模块用于优化小分子结构并输出结果

import os
import logging
from src.molecular import MolecularModule
from pysimm import system, forcefield, lmps

logging.basicConfig(level=logging.INFO)
def is_contains_file(directory, filename):
    """检查指定目录中是否包含给定文件名（忽略后缀）。
    
    Args:
        directory: 要检查的目录路径。
        filename: 要查找的文件名。
    
    Returns:
        如果目录中存在该文件名，则返回 True，否则返回 False。
    """
    # 获取不带后缀的文件名
    name_without_extension = os.path.splitext(filename)[0]
    # 列出目录中的所有文件和目录
    for item in os.listdir(directory):
        current_name_without_extension = os.path.splitext(item)[0]
        # 如果找到指定的文件名，返回 True
        if current_name_without_extension == name_without_extension:
            return True
    return False

def assign_mol_sys_ff(mol_path, file_name, data_out):
    """应用 GAFF2 力场并优化分子结构。
    
    Args:
        mol_path: 输入的分子文件路径。
        file_name: 输入的分子文件名（不带后缀）。
        data_out: 输出数据的路径。
    """
    pwd = os.getcwd()
    mol_file = f'{file_name}.mol'
    data_file = f'{data_out}{file_name}.data'
    try:
        os.chdir(mol_path)
        mol_system = system.read_mol(mol_file)
        mol_system.apply_forcefield(f=forcefield.Gaff2(), charges="gasteiger")
        lmps.quick_min(mol_system, np=1, min_style='cg', name='min_cg')
        mol_system.write_lammps(data_file)
        logging.info("优化分子结构成功")
    except Exception as e:
        logging.error(f"优化分子结构时出错: {e}")
        raise
    finally:
        os.chdir(pwd)

def assign_mol_sys_ff_ions(mol_path, file_name, data_out, ff_dict):
    """为离子分子分配力场并优化结构。
    
    Args:
        mol_path: 输入的离子分子文件路径。
        file_name: 输入的离子分子文件名（不带后缀）。
        data_out: 输出数据的路径。
        ff_dict: 力场字典，用于离子类型的分配。
    """
    pwd = os.getcwd()
    mol_file = f'{file_name}.mol'
    data_file = f'{data_out}{file_name}.data'
    try:
        os.chdir(mol_path)
        mol_system_ions = system.read_mol(mol_file)
        for i in mol_system_ions.particles:
            for key in ff_dict:
                if i.elem == key:
                    i.type_name = ff_dict.get(key)
                    i.set(bonds=None)
        mol_system_ions.bonds.count = 0
        mol_system_ions.apply_forcefield(f=forcefield.Gaff2())
        mol_system_ions.write_lammps(data_file)
        logging.info("指认离子力场成功")
    except Exception as e:
        logging.error(f"指认离子力场时出错: {e}")
        raise
    finally:
        os.chdir(pwd)
    
def generate_molecular_structure(smiles, file_name, mol_path, g_type):
    """生成分子结构、图片并保存到指定路径。
    
    Args:
        smiles: 分子的 SMILES 字符串。
        file_name: 输出文件名（不带后缀）。
        mol_path: 分子文件的输出路径。
        g_type: 分子类型:
                反应物: r1(氨基), r2(酰氯)
                溶剂: sol
                产物: p
                副产物: byp
    """
    logging.info(f'检测到新的{g_type}分子')
    mol = MolecularModule(smiles, g_type)
    mol.rkmol_print(mol.mol, mol.smiles, f'{mol_path}{file_name}')

def get_mol_info(path_dict, g_type):
    """根据分子类型获取文件名和 SMILES 字符串。
    
    Args:
        path_dict: 包含分子名称和 SMILES 的字典。
        g_type: 分子的类型。
    
    Returns:
        file_name: 分子文件名。
        smiles: 分子的 SMILES 字符串。
    """
    file_name = path_dict['names'][g_type]
    smiles = path_dict['smiles'][g_type]
    return file_name, smiles

def add_mol_data(path_dict, g_type, mol_path, data_out, mol_ff = {}):
    """根据分子类型添加分子数据并优化结构。
    
    Args:
        path_dict: 包含分子名称和 SMILES 的字典。
        g_type: 分子的类型。
        mol_path: 分子文件的输出路径。
        data_out: 输出数据的路径。
        mol_ff: 离子的 GAFF2 力场字典。
    """
    file_name, smiles = get_mol_info(path_dict, g_type)
    if not (g_type == 'p' or g_type == 'byp'):
        if not is_contains_file(mol_path, file_name):
            generate_molecular_structure(smiles, file_name, mol_path, g_type)
            assign_mol_sys_ff(mol_path, file_name, data_out)
    elif g_type == 'p':
        for name, smile in zip(file_name, smiles):
            if not is_contains_file(mol_path, name):
                generate_molecular_structure(smile, name, mol_path, g_type)
                assign_mol_sys_ff(mol_path, name, data_out)
    else:
        if not is_contains_file(mol_path, file_name):
            generate_molecular_structure(smiles, file_name, mol_path, g_type)
            assign_mol_sys_ff_ions(mol_path, file_name, data_out, mol_ff)
    logging.info(f"已添加{g_type}分子 DATA 文件至相应目录")

def optimize_structure(path_dict, mol_ff = {'Cl': 'Cl', 'H': 'hx'}):
    """根据反应物的 SMILES 优化分子结构。
    
    Args:
        path_dict: 包含反应物、溶剂及产物的 SMILES 和文件名的字典。
        data_path: 数据输出路径。
        mol_ff: 离子的力场字典（默认为 {'Cl': 'Cl', 'H': 'hx'}）。
    """
    mol_path = path_dict['paths']['mol']
    data_out = path_dict['paths']['data']

    # 输出 SMILES 结构到 data/mol 文件夹
    add_mol_data(path_dict, 'r1', mol_path, data_out)
    add_mol_data(path_dict, 'r2', mol_path, data_out)
    add_mol_data(path_dict, 'sol', mol_path, data_out)
    add_mol_data(path_dict, 'p', mol_path, data_out)
    add_mol_data(path_dict, 'byp', mol_path, data_out, mol_ff)
    