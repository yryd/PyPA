# single.py
# 该模块为单独运行模块，接收反应物和溶剂的SMILES字符串及其比例作为输入，并生成所有必要的目录和文件

import os
import logging
from src.filewriter import write_file
from src.generator import generate_reaction_smile

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def run_single_simulation():
    """
    运行单个反应模拟任务，初始化路径及反应信息作为输入。
    
    Args:
        path_dict: 包含基本目录结构、反应信息的字典
        
    Returns:
        path_dict: 包含基本目录结构、反应信息、产物信息、分子信息、运行结果的补充字典
    """
    
    # 示例SMILES字符串
    r1_smiles = "C1CNCCN1"  # 哌嗪
    r2_smiles = "O=C(Cl)C1=CC(C(=O)Cl)=CC(C(=O)Cl)=C1"  # TMC
    sol_smiles = "CCCCCC"  # 己烷作为溶剂
    
    from src.simulator import path_init
    # 先初始化模拟环境，获取环境变量path_dict
    path_dict = path_init(r1_smiles, r2_smiles, sol_smiles, reactant1_ratio=3, reactant2_ratio=2, solvent_ratio=1, num=50)
    from src.optimizer import init_mol_prop
    # 优化反应物分子结构，更新环境变量
    path_dict, r_mols = init_mol_prop(path_dict, path_dict['file_name'])
    
    # 生成优化产物与副产物分子结构，更新环境变量
    from src.generator import init_product_info
    path_dict, _ = init_product_info(path_dict, r_mols[0], r_mols[1])
    
    # Step 4: 构建分子系统
    from src.builder import construct_system, construct_reaction_map, prepare_reaction_templates
    construct_system(path_dict, box_len = 60)
    logging.info('体系坐标构建完成')
    
    # Step 5: 生成反应模板
    path_dict = prepare_reaction_templates(path_dict)
    path_dict = construct_reaction_map(path_dict)
    logging.info('反应映射模板构建完成')    
    
    # Step 6: 检测输入文件完整性，复制到运行目录，生成初始输入文件：sys_init.lmps
    from src.simulator import simulation_file_collect
    path_dict = simulation_file_collect(path_dict)
    logging.info('拓扑、映射等输入文件准备完成')
    
    # 生成Lammps输入文件
    from src.template import add_all_lammps_params, write_lammps_template
    path_dict = add_all_lammps_params(path_dict)
    write_lammps_template(path_dict)
    logging.info('Lammps输入文件生成完成')
    
    # Step 7: 运行Lammps模拟
    from src.lammpsrun import run_lammps_simulation
    
    # 运行LAMMPS模拟，可以根据需要调整参数
    run_lammps_simulation(path_dict, ntasks=2, use_gpu=True, ngpus=2, input_file=path_dict['params']['xlink']['x_file_name'], output_dir=path_dict['paths']['xlink'])
    run_lammps_simulation(path_dict, ntasks=2, use_gpu=True, ngpus=2, input_file=path_dict['params']['insert']['h2o_file_name'], output_dir=path_dict['paths']['insert'])
    
    
    
    print(path_dict)
    import json
    with open('data.json', 'w') as f:
        json.dump(path_dict, f, indent = 4)
    return path_dict

# 示例用法
if __name__ == "__main__":
    run_single_simulation()