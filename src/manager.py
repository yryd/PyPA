# manager.py
# 并行运行数据库单行反应模拟任务

import logging
import concurrent.futures
from src.database import DatabaseModule
from src.simulator import simulation_init
from src.optimizer import optimize_structure
from src.builder import construct_system, construct_reaction_map

from src.filewriter import write_lammps_in


from src.template_generator import generate_reaction_templates

from src.analyzer import analyze_diffusion


logging.basicConfig(level=logging.INFO)

def run_simulation(id, tmp_path = 'tmp/'):
    """
    运行单个反应模拟任务。
    
    Args:
        id: 数据库的反应 id 编号。
    """
    logging.info(f'反应id: {id}开始')
    # Step 1: 连接数据库
    db = DatabaseModule()
    logging.info('连接数据库完成')
    # Step 2: 初始化模拟的文件名、路径等，产生并优化小分子结构
    path_dict = simulation_init(db, id, tmp_path)
    logging.info('初始化模拟路径完成')
    
    # Step 3: 优化结构与指认力场
    optimize_structure(path_dict, mol_ff = {'Cl': 'Cl', 'H': 'hx'})
    logging.info('小分子结构优化完成')
    
    # Step 4: 构建分子系统
    construct_system(path_dict)
    logging.info('体系坐标构建完成')
    
    # Step 5: 生成反应模板
    construct_reaction_map(path_dict)
    logging.info('反应映射模板构建完成')
    
    # # Step 6: 生成 Lammps 输入文件
    # write_lammps_in(data_PATH, map_PATH, run_PATH)
    # logging.info('Lammps 输入文件写入完成')
    
    # # Step 6: 并行运行 LAMMPS 模拟
    
    # run_lammps_simulation("in.lammps")

    # # Step 7: 对运行结果分析扩散和孔径
    # analyze_diffusion("simulation_results.data")

    # # Step 8: 将结果存储在数据库中
    # simulation_data = (reactant1, reactant2, solvent, 10, 10, 100, 300, 1.0, 500, 1.23, 0.45)
    # store_simulation_results("simulation.db", simulation_data)
    
    # print(f"Running simulation for {id}")


def run_parallel(db_path='data/database/reaction.db', max_workers = 3):
    """
    使用线程池并行运行多个反应模拟任务。
    
    Args:
        db_path: 数据库文件路径，默认为 'data/database/reaction.db'。
        max_workers: 最大线程数，默认为 3。
    """
    db = DatabaseModule(db_path)
    indexs = db.get_column_list('id')
    with concurrent.futures.ThreadPoolExecutor(max_workers) as executor:
        executor.map(run_simulation, indexs)


if __name__ == "__main__":
    run_simulation(1, tmp_path = '../tmp', data_path = '..data/')