# manager.py
# 并行运行数据库单行反应模拟任务

import logging
import concurrent.futures
from src.database import DatabaseModule
from src.simulator import *
from src.optimizer import optimize_structure
from src.builder import construct_system, construct_reaction_map

from pysimm.system import System, read_lammps
from pysimm.lmps import Simulation

from src.filewriter import write_lammps_in

from pysimm import forcefield


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
    construct_system(path_dict, box_len = 60)
    logging.info('体系坐标构建完成')
    
    # Step 5: 生成反应模板
    construct_reaction_map(path_dict)
    logging.info('反应映射模板构建完成')
    
    # Step 6: 检测输入文件完整性，复制到运行目录
    simulation_file_collect(path_dict)
    logging.info('拓扑、映射文件准备完成')

    # Step 7: 新建模拟，生成 Lammps 输入文件
    sys_params = {
        "quiet": False,  # 是否打印状态信息
        "atom_style": "full",  # 原子样式
        "pair_style": "lj/cut",  # 配对样式
        "bond_style": "harmonic",  # 键样式
        "angle_style": "harmonic",  # 角样式
        "dihedral_style": "fourier",  # 二面角样式
        "improper_style": "cvff",  # 不当样式（共轭平面性）
        "set_types": True,  # 是否将类型对象化
        "name": "System_all"  # 系统名称
    }
    system = SystemModule(path_dict, sys_params)

    simulation_params = {
        'minimization_enabled': True,  # 是否进行最小化
        'minimization_params': {
            'min_style': 'cg',  # 最小化方法
            'etol': 1e-5,  # 能量容忍度
            'ftol': 1e-5,  # 力容忍度
            'maxiter': 10000,  # 最大迭代次数
            'maxeval': 100000,  # 最大评估次数
            'unfix': True  # 使用后解除 fix
        },
        'velocity_enabled': True, # 是否进行速度随机
        'velocity': {
            'temperature': 300.0,  # 初始速度温度（单位：K）
            'seed': 4928459,  # 随机种子
            'distribution': 'gaussian'  # 速度分布类型
        },
        'thermodynamics_enabled': True,  # 是否进行松弛
        'thermodynamics_params': {
            'ensemble': 'nvt',  # 松弛动力学系综类型
            'temperature': {  # 温度相关参数
                'start': 298.15,  # 初始温度（单位：K）
                'stop': 298.15,  # 结束温度（单位：K）
                'damp': 100.0,  # 温度阻尼
            },
            'pressure': {  # 压力相关参数
                'iso': 'iso',  # 压力类型
                'start': 1.0,  # 初始压力（单位：atm）
                'stop': 1.0,  # 结束压力（单位：atm）
                'damp': 1000.0,  # 压力阻尼
            },
            'time_step': 1.0,  # 时间步长（单位：fs）
            'total_time': 100000,  # 总模拟时间（单位：fs）
            'unfix': True  # 使用后解除 fix
        },
        'output': {
            'dump_file': 'traj_npt.xtc',  # 轨迹文件名称
            'dump_frequency': 100,  # 轨迹输出频率
            'thermo_frequency': 100  # 热力学输出频率
        },
        'restart': {
            'restart_file': 'system_end.rst',  # 重启文件名称
            'data_file': 'system_end.lmp'  # 数据文件名称
        }
    }
    simulation = SimulationModule(path_dict, system)
    simulation.get_lammps_in(simulation_params, gpu = False)
    logging.info('Lammps 输入文件写入完成')
    
    # # Step 8: 运行 LAMMPS 模拟
    os.chdir(path_dict['paths']['run'])
    os.system('mpirun -np 2 lmp_mpi -echo screen -in in.lammps')

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