# simulator.py
# 该模块为模拟模块，包括初始化路径、创建写入文件路径、执行模拟等

import os
import logging
import ast
import time
from src.filewriter import write_file, combin_files
from pysimm.system import System, read_lammps
from pysimm.lmps import Simulation
from src.readdata import read_sys_info

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

        product_smiles = "[H]c1c(C(=O)Cl)c([H])c(C(=O)N2C([H])([H])C([H])([H])N([H])C([H])([H])C2([H])[H])c([H])c1C(=O)Cl;[H]c1c(C(Cl)=O)c([H])c(C(N2C([H])([H])C([H])([H])N(Cc3c(c(C(Cl)=O)c(c(c3[H])C(Cl)=O)[H])[H])C([H])([H])C2([H])[H])=O)c([H])c1C(Cl)=O"
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
        'run': os.path.join(reaction_path, 'run') + '/',
        'logs': os.path.join(reaction_path, 'logs') + '/',
        'result': os.path.join(reaction_path, 'run') + '/result/'
    }

    # 创建所有子目录
    for path in paths.values():
        os.makedirs(path, exist_ok=True)
        logging.info(f"创建目录: {path}")

    # 将路径信息添加到 path_dict 中
    path_dict['paths'] = paths
    return path_dict

def path_init(reactant1_smiles, reactant2_smiles, solvent_smiles, 
                          reactant1_ratio, reactant2_ratio, solvent_ratio,
                          num, tmp_path='tmp/', id = "test", reactant1_key = "1",
                          reactant2_key = "1", solvent_key = "1"):
    """
    初始化模拟，直接接收反应物和溶剂的SMILES字符串及其比例作为输入，创建文件路径，存储反应结构与条件。

    Args:
        reactant1_smiles: 反应物1的SMILES字符串
        reactant2_smiles: 反应物2的SMILES字符串
        solvent_smiles: 溶剂的SMILES字符串
        reactant1_ratio: 反应物1的比例
        reactant2_ratio: 反应物2的比例
        solvent_ratio: 溶剂的比例
        num: 反应物的分子数量基数
        tmp_path: 临时文件路径，用于存放模拟结果

    Returns:
        path_dict: 包含基本目录、反应条件、反应物结构的字典
    """
    logging.info('初始化模拟路径...')    
    # 创建反应系统文件路径
    base_PATH = os.getcwd()
    # 确保路径为绝对路径
    base_PATH = os.path.abspath(base_PATH)
    reaction_path = os.path.join(base_PATH, tmp_path, id)

    os.makedirs(reaction_path, exist_ok=True)

    # 创建各个子目录路径
    paths = {
        'lt': os.path.join(reaction_path, 'lt') + '/',
        'mol': os.path.join(reaction_path, 'mol') + '/',
        'sys': os.path.join(reaction_path, 'sys') + '/',
        'map': os.path.join(reaction_path, 'map') + '/',
        'data': os.path.join(reaction_path, 'data') + '/',
        'run': os.path.join(reaction_path, 'run') + '/',
        'logs': os.path.join(reaction_path, 'logs') + '/',
        'result': os.path.join(reaction_path, 'run') + '/result/',
        'xlink': os.path.join(reaction_path, 'run') + '/xlink/',
        'insert': os.path.join(reaction_path, 'run') + '/insert/',
        'MSD': os.path.join(reaction_path, 'run') + '/MSD/'
    }

    # 创建所有子目录
    for path in paths.values():
        os.makedirs(path, exist_ok=True)
        logging.info(f"创建目录: {path}")

    # 计算各组分的分子数量
    r1_num = reactant1_ratio * num
    r2_num = reactant2_ratio * num
    sol_num = solvent_ratio * num
    
    # 创建反应物基本信息字典
    # 使用文件名作为键
    r1_name = f'r1_{reactant1_key}'
    r2_name = f'r2_{reactant2_key}'
    sol_name = f'sol_{solvent_key}'
    # file_name 按顺序排放r1, r2, sol, p, byp 其中p可能多个
    path_dict = {
        'file_name': [r1_name, r2_name, sol_name],
        'type': {
            r1_name: 'r1',
            r2_name: 'r2',
            sol_name: 'sol'
        },
        'smiles': {
            r1_name: reactant1_smiles,
            r2_name: reactant2_smiles,
            sol_name: solvent_smiles
        },
        'num': {
            r1_name: r1_num,
            r2_name: r2_num,
            sol_name: sol_num
        }
    }


    # 将路径信息添加到 path_dict 中
    path_dict['paths'] = paths
    path_dict['paths']['template'] = base_PATH + '/data/template/'
    logging.info('初始化模拟路径完成')
    
    return path_dict


def simulation_file_collect(path_dict):
    data_PATH = path_dict['paths']['data']
    map_PATH = path_dict['paths']['map']
    xlink_PATH = path_dict['paths']['xlink']
    str_list = []
    str_list.append(f'cp {data_PATH}cleanedsystem.data {xlink_PATH}system.data\n')
    str_list.append(f'cp {data_PATH}cleanedsystem.in.settings {xlink_PATH}system.in.settings\n')
    # str_list.append(f'cp {data_PATH}system.in.init {xlink_PATH}system.in.init\n')
    for map_name in path_dict['map_templates'].keys():
        str_list.append(f'cp {map_PATH}mol.pre_{map_name} {xlink_PATH}mol.pre_{map_name}\n')
        str_list.append(f'cp {map_PATH}mol.post_{map_name} {xlink_PATH}mol.post_{map_name}\n')
        str_list.append(f'cp {map_PATH}txt.{map_name} {xlink_PATH}txt.{map_name}\n')
    
    write_file(xlink_PATH, 'getfile.sh', str_list)
    
    try:
        os.system(f'. {xlink_PATH}getfile.sh > /dev/null 2>&1')
        # 将 system.data 和 system.in.settings 的内容合并到 sys_init.lmps
        combin_files(xlink_PATH)
    except Exception as e:
        logging.error(f"校验生成输入文件完整性出错: {e}")
        raise

    sys_info = read_sys_info(xlink_PATH + 'sys_init.lmps')
    path_dict['sys_info'] = sys_info
    logging.info('模拟文件收集完成')
    return path_dict


# 继承 pysimm 包 System 类
class SystemModule(System):
    def __init__(self, path_dict, params):
        super().__init__()
        self.run_PATH = path_dict['paths']['run']
        self.data_all = self.run_PATH + 'sys_init.lmps'
        self.updata_lmps(params)
    
    def updata_lmps(self, params):
        pysimm_system = read_lammps(self.data_all, **params)
        # 将 pysimm_system 的属性更新到当前 SystemModule 实例
        self.__dict__.update(pysimm_system.__dict__)

# 继承 pysimm 包 Simulation 类
class SimulationModule(Simulation):
    def __init__(self, path_dict, system):
        super().__init__(system, log= f'{path_dict["paths"]["logs"]}steps.log', custom  = True)
        self.run_PATH = path_dict['paths']['run']
        self.result_PATH = path_dict['paths']['result']
        self.map_PATH = path_dict['paths']['map']
        self.molnum = path_dict['num']
        self.molsnames = path_dict['names']
        
    @staticmethod
    def init_writer():
        init_str = ''
        init_str += '#' * 80 + '\n'
        init_str += "units           real\n"
        init_str += "pair_style      lj/cut 12.0\n"
        init_str += "atom_style      full\n"
        init_str += "bond_style      harmonic\n"
        init_str += "angle_style     harmonic\n"
        init_str += "dihedral_style  fourier\n"
        init_str += "improper_style  cvff\n"
        init_str += "\n"
        init_str += "dimension       3\n"
        init_str += "boundary        p p p\n"
        init_str += "neigh_modify    every 1 delay 0 check yes\n"
        init_str += "neighbor        2.5 bin\n"
        init_str += "read_data       sys_init.lmps\n"
        init_str += '#' * 80 + '\n'
        return init_str

    def input_conditions(self, simulation_params):
        self.input_conditions_start(simulation_params)
        self.input_conditions_custum()
        self.input_conditions_end(simulation_params)

    def input_conditions_start(self, simulation_params):
        # 初始化
        init_str = self.init_writer()
        self.add_custom(init_str)
        
        self.add_custom(f"write_data      {self.result_PATH}system_start.lmp")
        # 能量最小化
        if simulation_params.get('minimization_enabled'):
            self.add_min(name = 'min_sys', **simulation_params['minimization_params'])
        # 速度设置部分
        if simulation_params.get('minimization_enabled'):
            self.add_custom(f"\n\nvelocity        all create {simulation_params['velocity']['temperature']} {simulation_params['velocity']['seed']} rot yes dist {simulation_params['velocity']['distribution']}\n\n")

        # 松弛结构
        if simulation_params.get('thermodynamics_enabled'):
            self.add_md(name = 'relax_md', **simulation_params['thermodynamics_params'])

        # 输出
        self.add_custom(f"\ndump            mydump1 all xtc {simulation_params['output']['dump_frequency']} {self.result_PATH}{simulation_params['output']['dump_file']}\n\n")

    def input_conditions_custum(self):
        react_len = 3.25
        names = self.molsnames
        nums = self.molnum
        r1_name = names['r1']
        r2_name = names['r2']
        sol_name = names['sol']
        p_names_list = names['p']
        byp_name = names['byp']
        r1_num = nums['r1']
        r2_num = nums['r2']
        sol_num = nums['sol']
        post_type = len(p_names_list)
        
        str_tmp0 = ""
        str_tmp0 += f"molecule        pre {self.map_PATH}pre_mol.data\n"
        for i in range(post_type):
            str_tmp0 += f"molecule        post_{i} {self.map_PATH}post_{i}_mol.data\n"
        self.add_custom(str_tmp0)
        
        str_tmp1 = ""
        str_tmp1 += "timestep        1\n"
        str_tmp1 += "fix             md_normal  all npt temp 298.15 298.15 100.0 iso 1.0 1.0 1000.0\n"
        str_tmp1 += "run             2000\n"
        str_tmp1 += "unfix           md_normal\n"
        self.add_custom(str_tmp1)
        
        str_tmp2 = ""
        str_tmp2 += "fix             xlink_fix all bond/react stabilization yes statted_grp .03 &\n"
        for i in range(post_type):
            # 每 100 步反应一次，反应距离为 [0, react_len]
            str_tmp2 += f"                    react rxn1 all 100 0 {react_len} pre post_{i} {self.map_PATH}automap_{i}.data stabilize_steps 100\n"
        str_tmp2 += "fix             nvt_md all nvt temp 300.0 300.0 100.0\n"
        str_tmp2 += "run             100000\n"
        str_tmp2 += "unfix           nvt_md\n"
        str_tmp2 += "fix             npt_md  all npt temp 298.15 298.15 100.0 iso 1.0 1.0 1000.0\n"
        str_tmp2 += "run             100000\n"
        str_tmp2 += "unfix           npt_md\n"
        str_tmp2 += "unfix           xlink_fix\n"
        self.add_custom(str_tmp2)
        
        
    def input_conditions_end(self, simulation_params):
        str_end = ''
        str_end += "\n\n"
        str_end += f"write_restart   {self.result_PATH}{simulation_params['restart']['restart_file']}\n"
        str_end += f"write_data      {self.result_PATH}{simulation_params['restart']['data_file']}\n\n"
        self.add_custom(str_end)

    @staticmethod
    def add_gpu(input_str):
        # 定义需要查找的关键词
        keywords = [' lj/cut ', ' lj/cut/coul/long ', ' nvt ', ' npt ', ' nve ']

        # 遍历关键词并替换
        for keyword in keywords:
            if keyword in input_str:
                # 替换关键词，添加后缀 '/gpu'
                input_str = input_str.replace(keyword, f"{keyword[:-1]}/gpu ")
        return input_str
    
    def get_lammps_in(self, simulation_params, gpu = False):
        
        self.input_conditions(simulation_params)
        write_str = self.input
        if gpu:
            write_str = self.add_gpu(self.input)
        with open(f'{self.run_PATH}in.lammps', 'w', encoding='utf-8') as file:
            file.write(write_str)
