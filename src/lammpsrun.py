# lammpsrun.py
# 该模块为LAMMPS运行模块，包括运行单个LAMMPS模拟和并行运行多个LAMMPS模拟

import os
import logging
import subprocess
import multiprocessing

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_lammps_simulation(path_dict, ntasks=4, use_gpu=True, ngpus=1, binsize=None, input_file=None, output_dir=None):
    """
    运行LAMMPS模拟
    
    Args:
        path_dict: 包含路径信息的字典
        ntasks: MPI进程数量
        use_gpu: 是否使用GPU加速
        gpu_options: GPU选项，默认为"1"
        binsize: 邻居列表bin大小，如果为None则不设置
        input_file: 输入文件的完整路径，如果为None则使用默认路径和文件名
        output_dir: 输出文件的目录路径，如果为None则使用默认路径
        
    Returns:
        path_dict: 更新后的路径字典，包含模拟结果信息
    """
    logging.info('开始运行LAMMPS模拟...')
    
    # 获取运行目录
    run_PATH = output_dir if output_dir is not None else path_dict['paths']['run']
    input_file = input_file if input_file is not None else path_dict['params']['xlink']['x_file_name']
    # 构建LAMMPS命令
    cmd = ['mpirun', '-np', str(ntasks), 'lmp_mpi']
    
    # 添加GPU加速选项
    if use_gpu:
        cmd.extend(['-sf', 'gpu', '-pk', 'gpu', str(ngpus)])
        if binsize is not None:
            cmd.extend(['binsize', str(binsize)])
    
    # 添加输入文件
    cmd.extend(['-in', input_file])

    
    # 切换到运行目录
    current_dir = os.getcwd()
    os.chdir(run_PATH)
    
    try:
        # 运行LAMMPS命令
        logging.info(f"执行命令: {' '.join(cmd)}")
        process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        # 实时输出日志
        while True:
            output = process.stdout.readline()
            if output == '' and process.poll() is not None:
                break
            if output:
                logging.info(output.strip())
        
        # 获取返回码
        return_code = process.poll()
        
        # 检查是否成功
        if return_code == 0:
            logging.info('LAMMPS模拟成功完成')
        else:
            stderr = process.stderr.read()
            logging.error(f'LAMMPS模拟失败，返回码: {return_code}\n错误信息: {stderr}')
    
    except Exception as e:
        logging.error(f"运行LAMMPS时出错: {e}")
    
    finally:
        # 恢复原始工作目录
        os.chdir(current_dir)
    
    return

def run_parallel_simulations(simulation_configs, max_workers=None):
    """
    并行运行多个LAMMPS模拟任务
    
    Args:
        simulation_configs: 包含多个模拟配置的列表，每个配置是一个字典，包含path_dict和其他参数
        max_workers: 最大并行任务数，默认为CPU核心数
        
    Returns:
        results: 包含所有模拟结果的列表
    """
    logging.info(f'开始并行运行{len(simulation_configs)}个LAMMPS模拟任务...')
    
    # 如果未指定最大并行任务数，则使用CPU核心数
    if max_workers is None:
        max_workers = multiprocessing.cpu_count()
    
    # 限制最大并行任务数不超过配置数量
    max_workers = min(max_workers, len(simulation_configs))
    
    # 创建进程池
    pool = multiprocessing.Pool(processes=max_workers)
    
    # 提交任务
    results = []
    for config in simulation_configs:
        path_dict = config.get('path_dict')
        ntasks = config.get('ntasks', 4)
        use_gpu = config.get('use_gpu', True)
        gpu_options = config.get('gpu_options', "1")
        binsize = config.get('binsize', None)
        
        # 异步提交任务
        result = pool.apply_async(run_lammps_simulation, 
                                 args=(path_dict, ntasks, use_gpu, gpu_options, binsize))
        results.append(result)
    
    # 关闭进程池，不再接受新任务
    pool.close()
    
    # 等待所有任务完成
    pool.join()
    
    # 获取所有任务的结果
    final_results = [result.get() for result in results]
    
    logging.info('所有LAMMPS模拟任务已完成')
    return final_results