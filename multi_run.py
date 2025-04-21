# multi_run.py
# 该模块用于并行运行多个LAMMPS模拟任务

import os
import logging
import argparse
from src.simulator import path_init, simulation_file_collect, run_lammps_simulation, run_parallel_simulations
from src.optimizer import init_mol_prop
from src.generator import init_product_info
from src.builder import construct_system, construct_reaction_map, prepare_reaction_templates
from src.template import add_all_lammps_params, write_lammps_template

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def prepare_simulation(r1_smiles, r2_smiles, sol_smiles, reactant1_ratio, reactant2_ratio, solvent_ratio, num, task_id):
    """
    准备单个模拟任务
    
    Args:
        r1_smiles: 反应物1的SMILES字符串
        r2_smiles: 反应物2的SMILES字符串
        sol_smiles: 溶剂的SMILES字符串
        reactant1_ratio: 反应物1的比例
        reactant2_ratio: 反应物2的比例
        solvent_ratio: 溶剂的比例
        num: 反应物的分子数量基数
        task_id: 任务ID
        
    Returns:
        path_dict: 包含模拟环境的字典
    """
    # 初始化模拟环境
    path_dict = path_init(r1_smiles, r2_smiles, sol_smiles, 
                          reactant1_ratio, reactant2_ratio, solvent_ratio, 
                          num, tmp_path='tmp/', id=str(task_id))
    
    # 优化反应物分子结构
    path_dict, r_mols = init_mol_prop(path_dict, path_dict['file_name'])
    
    # 生成优化产物与副产物分子结构
    path_dict, _ = init_product_info(path_dict, r_mols[0], r_mols[1])
    
    # 构建分子系统
    construct_system(path_dict, box_len=60)
    logging.info('体系坐标构建完成')
    
    # 生成反应模板
    path_dict = prepare_reaction_templates(path_dict)
    path_dict = construct_reaction_map(path_dict)
    logging.info('反应映射模板构建完成')
    
    # 检测输入文件完整性，复制到运行目录
    path_dict = simulation_file_collect(path_dict)
    logging.info('拓扑、映射等输入文件准备完成')
    
    # 生成Lammps输入文件
    path_dict = add_all_lammps_params(path_dict)
    write_lammps_template(path_dict)
    logging.info('Lammps输入文件生成完成')
    
    return path_dict

def run_multi_simulations(tasks, max_workers=None):
    """
    运行多个模拟任务
    
    Args:
        tasks: 任务列表，每个任务是一个字典，包含模拟参数
        max_workers: 最大并行任务数
        
    Returns:
        results: 包含所有模拟结果的列表
    """
    # 准备所有模拟任务的配置
    simulation_configs = []
    
    for i, task in enumerate(tasks):
        # 准备模拟环境
        path_dict = prepare_simulation(
            task['r1_smiles'],
            task['r2_smiles'],
            task['sol_smiles'],
            task['reactant1_ratio'],
            task['reactant2_ratio'],
            task['solvent_ratio'],
            task['num'],
            task.get('id', f'task_{i}')
        )
        
        # 创建模拟配置
        config = {
            'path_dict': path_dict,
            'ntasks': task.get('ntasks', 4),
            'use_gpu': task.get('use_gpu', True),
            'gpu_options': task.get('gpu_options', "1"),
            'binsize': task.get('binsize', None)
        }
        
        simulation_configs.append(config)
    
    # 并行运行所有模拟任务
    results = run_parallel_simulations(simulation_configs, max_workers)
    
    return results

def main():
    parser = argparse.ArgumentParser(description='并行运行多个LAMMPS模拟任务')
    parser.add_argument('--max_workers', type=int, default=None, help='最大并行任务数')
    args = parser.parse_args()
    
    # 示例任务列表
    tasks = [
        {
            # 任务1：哌嗪与TMC反应，溶剂为己烷
            'r1_smiles': "C1CNCCN1",  # 哌嗪
            'r2_smiles': "O=C(Cl)C1=CC(C(=O)Cl)=CC(C(=O)Cl)=C1",  # TMC
            'sol_smiles': "CCCCCC",  # 己烷
            'reactant1_ratio': 3,
            'reactant2_ratio': 2,
            'solvent_ratio': 1,
            'num': 50,
            'id': 'task_1',
            'ntasks': 4,
            'use_gpu': True
        },
        # 可以添加更多任务...
    ]
    
    # 运行多个模拟任务
    results = run_multi_simulations(tasks, args.max_workers)
    
    # 输出结果
    for i, result in enumerate(results):
        status = result['simulation_result']['status']
        logging.info(f"任务 {i+1} 状态: {status}")
        
        if status == 'success':
            output_files = result['simulation_result']['output_files']
            logging.info(f"输出文件: {output_files}")
        else:
            error = result['simulation_result']['error']
            logging.error(f"错误信息: {error}")

if __name__ == "__main__":
    main()