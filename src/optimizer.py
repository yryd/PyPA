# optimizer.py
# 该模块用于优化小分子结构并输出结果，以及初始化分子属性

import os
import logging
from src.molecular import MolecularModule
from src.execute import exec_ltemplify
from pysimm import system, forcefield, lmps

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

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

def apply_forcefield_and_optimize(mol_path, file_name, data_out, ff_dict=None):
    """分子应用GAFF2力场并优化分子结构，支持普通分子和含力场指认的离子化合物。
    
    Args:
        mol_path: 输入的分子文件路径。
        file_name: 输入的分子文件名（不带后缀）。
        data_out: 输出数据的路径。
        ff_dict: 离子的力场字典。
    """
    pwd = os.getcwd()
    mol_file = f'{file_name}.mol'
    data_file = f'{data_out}{file_name}.data'
    
    try:
        os.chdir(mol_path)
        mol_system = system.read_mol(mol_file)
        
        if ff_dict:
            # 处理离子化合物
            for particle in mol_system.particles:
                for key in ff_dict:
                    if particle.elem == key:
                        particle.type_name = ff_dict.get(key)
                        particle.set(bonds=None)
            mol_system.bonds.count = 0
            mol_system.apply_forcefield(f=forcefield.Gaff2())
            logging.info("指认离子力场成功")
        else:
            # 处理普通分子
            mol_system.apply_forcefield(f=forcefield.Gaff2(), charges="gasteiger")
            lmps.quick_min(mol_system, np=1, min_style='cg', name='min_cg')
            logging.info("优化分子结构成功")
            
        mol_system.write_lammps(data_file)
    except Exception as e:
        logging.error(f"{'指认离子力场' if ff_dict else '优化分子结构'}时出错: {e}")
        raise
    finally:
        os.chdir(pwd)

def create_molecule_file(file_name, g_type, obj_mol, mol_path, data_out, mol_ff=None):
    """处理单类分子：生成结构文件并应用力场。
    
    Args:
        file_name: 分子文件名。
        g_type: 分子类型。
        obj_mol: 分子对象。
        mol_path: 分子文件输出路径。
        data_out: 数据输出路径。
        mol_ff: 离子力场字典，默认为None。
    """
    # 处理所有类型分子
    if file_name and not is_contains_file(mol_path, file_name):
        logging.info(f'检测到新的{g_type}分子')
        obj_mol.rkmol_print(obj_mol.mol, obj_mol.smiles, f'{mol_path}{file_name}')
        apply_forcefield_and_optimize(mol_path, file_name, data_out, mol_ff)
        logging.info(f"已添加{g_type}分子 DATA 文件至相应目录")

def init_mol_prop(path_dict, file_names, mol_ff = None):
    """
    初始化分子属性，通过SMILES实例化分子并计算其属性，优化并生成结构
    
    Args:
        path_dict: 包含基本目录结构、反应信息的字典
        file_names: 分子文件名列表
        mol_ff: 离子力场字典，默认为None
        
    Returns:
        path_dict: 包含基本目录结构、反应信息、分子属性的补充字典
        mols: 分子对象列表
    """
    logging.info('实例化反应物分子...')

    # 通过 SMILES 实例化反应物 MolecularModule 类
    mols = []
    for name in file_names:
        mol = MolecularModule(path_dict['smiles'][name], path_dict['type'][name])
        # 仅对反应物和溶剂分子计算属性
        if path_dict['type'][name] in ['r1', 'r2', 'sol']:
            # 显式调用计算分子属性
            mol.cal_mol_prop()
        # 优化生成分子结构
        create_molecule_file(name, path_dict['type'][name], mol, path_dict['paths']['mol'], path_dict['paths']['data'], mol_ff)
        exec_ltemplify(path_dict['paths']['data'], path_dict['paths']['lt'], name)
        mols.append(mol)
    logging.info('实例化反应物分子完成')
    
    # 将分子性质添加到 path_dict 中
    # 构造各个性质的临时字典
    molecular_properties = [
        'group_smiles', 'rings', 'count_group', 'group_min_distances',
        'group_max_distances', 'ar_sp3_balance', 'aromatic_rings', 'molfinger'
    ]
    
    # 使用字典推导式构建属性字典
    property_dicts = {
        prop: {name: getattr(mol, prop) for name, mol in zip(file_names, mols)}
        for prop in molecular_properties
    }
    
    # 将属性字典更新到path_dict中
    for prop in property_dicts:
        if prop in path_dict:
            path_dict[prop].update(property_dicts[prop])
        else:
            path_dict[prop] = property_dicts[prop]
    
    return path_dict, mols


if __name__ == "__main__":
    # 测试配置 - 验证分子初始化、力场应用和优化流程
    file_names = ['r1_1', 'r2_1', 'sol_1']
    mol_ff = None  # 离子力场字典，测试使用None表示普通分子
    
    path_dict = {
        'file_name': file_names,
        'type': {
            'r1_1': 'r1',
            'r2_1': 'r2',
            'sol_1': 'sol'
        },
        'smiles': {
            'r1_1': "C",
            'r2_1': "C", 
            'sol_1': "C"
        },
        'num': {
            'r1_1': 100,
            'r2_1': 100,
            'sol_1': 100
        },
        'paths': {
            'mol': './mol/',
            'data': './data/'
        }
    }
    
    # 初始化分子属性并验证流程
    path_dict, mols = init_mol_prop(path_dict, file_names, mol_ff)
    logging.info("测试完成 - 分子初始化流程验证通过")
    
