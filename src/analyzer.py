# analyzer.py
# 该模块用于分析分子动力学模拟结果，包括扩散系数计算和孔径分析

import os
import logging
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress
import MDAnalysis as mda
from MDAnalysis.analysis import msd
import csv
import sys

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

# 设置递归深度限制，用于孔径计算中的深度优先搜索
sys.setrecursionlimit(1000000)


class PoreSizeAnalyzer:
    """孔径分析器，用于计算多孔材料的孔径分布和平均孔径"""
    
    def __init__(self, arr):
        """初始化孔径分析器
        
        Args:
            arr: 布尔类型的三维数组，表示材料的结构
        """
        self.white = True
        self.row_num = arr.shape[0]
        self.col_num = arr.shape[0]
        self.walked_set = set()
        self.roming_set = set()
        self.dfs_num = 0
        self.array = arr.astype(bool)
        self.list = []
        self.xy_list = []

    def dfs(self, x, y, rgb):
        """深度优先搜索算法，用于寻找连通区域
        
        Args:
            x: x坐标
            y: y坐标
            rgb: 目标值
        """
        self.roming_set.add(tuple([x, y]))
        if tuple([x,y]) in self.walked_set:  # 重复遍历检查
            return
        if rgb != self.array[x][y]:  # 目标值检查
            return

        self.walked_set.add(tuple([x, y]))
        if (x == self.col_num - 1):
            self.dfs(0, y, rgb)  # x
        else:
            self.dfs(x + 1, y, rgb)  # x
        if (y == self.row_num - 1):
            self.dfs(x, 0, rgb)  # y
        else:
            self.dfs(x, y + 1, rgb)  # y
        if (x == - self.col_num):
            self.dfs(self.col_num - 1, y, rgb)  # -x
        else:
            self.dfs(x - 1, y, rgb)  # -x
        if (y == -self.row_num):
            self.dfs(x, self.row_num - 1, rgb)  # -y
        else:
            self.dfs(x, y - 1, rgb)  # -y
        return

    def walk(self):
        """遍历整个数组，寻找所有连通区域
        
        Returns:
            tuple: (list, xy_list) 孔径大小列表和对应的起始坐标
        """
        for y in range(self.col_num):
            for x in range(self.row_num):
                rgb = self.array[x][y]
                if tuple([x, y]) in self.roming_set:
                    continue
                if rgb != self.white:
                    self.dfs(x, y, rgb)
                    num = len(self.walked_set)
                    self.list.append(num)
                    self.xy_list.append([x,y])
                    self.walked_set.clear()
        self.roming_set.clear()
        return (self.list, self.xy_list)


def read_xyz_file(file_path):
    """读取XYZ文件并转换为CSV格式
    
    Args:
        file_path: XYZ文件的路径
        
    Returns:
        str: 生成的CSV文件路径
    """
    logging.info(f'读取XYZ文件: {file_path}')
    data = []
    with open(file_path, encoding='utf-8') as file:
        for line in file:
            items = line.split()
            # 去除数字行与空行
            if len(items) > 2:
                name = items[0]
                X = float(items[1])
                Y = float(items[2])
                Z = float(items[3])
                data.append([name, X, Y, Z])
    
    # 写入CSV文件
    output_path = file_path[:-4] + '_xyz.csv'
    headers = ['Name', 'X', 'Y', 'Z']
    rows = []
    for i in data:
        name = i[0]
        X = i[1]
        Y = i[2]
        Z = i[3]
        row = [name, X, Y, Z]
        rows.append(row)
    
    with open(output_path, 'w', newline='') as f:
        f_csv = csv.writer(f)
        f_csv.writerow(headers)
        f_csv.writerows(rows)
    
    logging.info(f'XYZ文件已转换为CSV: {output_path}')
    return output_path


def read_csv_data(file_path):
    """读取CSV文件数据
    
    Args:
        file_path: CSV文件路径
        
    Returns:
        numpy.ndarray: 读取的数据
    """
    with open(file_path, encoding='utf-8') as f:
        # 读取文件，支持str，跳过首行，读取索引为0,1,2,3的列
        data = np.loadtxt(f, str, delimiter=",", skiprows=1, usecols=(0,1,2,3))
        return data


def process_csv_data(data, atom_types):
    """处理CSV数据，筛选指定原子类型
    
    Args:
        data: CSV数据
        atom_types: 要筛选的原子类型列表
        
    Returns:
        numpy.ndarray: 处理后的数据
    """
    index = data[..., 0]
    content = data[..., 1:]
    all_list = []
    for name in atom_types:
        row = np.where(index == name)
        all_list = all_list + row[0].tolist()
    use_line = np.array(all_list)
    content = content[use_line, ...]
    # 注意以字符串读取进来的要转换回数据
    rev_data = content.astype(np.float64)
    return rev_data


def create_3d_histogram(data, cell_length, resolution):
    """创建三维频率直方图
    
    Args:
        data: 坐标数据
        cell_length: 晶格长度
        resolution: 分辨率
        
    Returns:
        tuple: (density, density_bool) 密度数组和布尔数组
    """
    len_int = int(cell_length)
    num_int = int(cell_length / resolution) + 1
    gridx = np.linspace(0, len_int, num_int)
    gridy = np.linspace(0, len_int, num_int)
    gridz = np.linspace(0, len_int, num_int)
    density, edges = np.histogramdd(data, bins=[gridx, gridy, gridz])
    density_bool = density.astype(bool)
    return (density, density_bool)


def analyze_plane(plane_data):
    """分析平面数据
    
    Args:
        plane_data: 平面数据
        
    Returns:
        tuple: (zero_in_vol, var3) 空隙率和方差
    """
    index_list = np.nonzero(plane_data)
    none_zero = index_list[0].size
    length = plane_data.shape[0]
    zero = length * length - none_zero
    zero_in_vol = zero / (length * length)
    
    one_line = []
    for deny in range(length):
        for denx in range(length):
            if plane_data[denx][deny] != 0:
                one_line.append(plane_data[denx][deny])
    var3 = np.var(np.asarray(one_line))
    return (zero_in_vol, var3)


def calculate_pore_size(data):
    """计算孔径
    
    Args:
        data: 数据数组
        
    Returns:
        float: 平均孔径
    """
    pore_size_obj = PoreSizeAnalyzer(data)
    (size_list, xy_list) = pore_size_obj.walk()
    myset = set(size_list)
    size_all = 0
    count = 0
    for item in myset:
        size_all += item * size_list.count(item)
        count += size_list.count(item)
    average = size_all / count if count > 0 else 0
    return average


def calculate_average_pore_size(xyz_file_path, atom_types, cell_length, resolution):
    """计算平均孔径
    
    Args:
        xyz_file_path: XYZ文件路径
        atom_types: 要分析的原子类型列表
        cell_length: 晶格长度
        resolution: 分辨率
        
    Returns:
        dict: 包含各方向平均孔径和孔隙率的字典
    """
    logging.info(f'开始计算平均孔径: {xyz_file_path}')
    
    # 如果输入是xyz文件，先转换为csv
    if xyz_file_path.endswith('.xyz'):
        csv_file_path = read_xyz_file(xyz_file_path)
    else:
        csv_file_path = xyz_file_path
    
    # 读取CSV数据
    csv_data = read_csv_data(csv_file_path)
    
    # 处理数据
    processed_data = process_csv_data(csv_data, atom_types)
    
    # 创建3D直方图
    density, density_bool = create_3d_histogram(processed_data, cell_length, resolution)
    data_int = density.astype(int)
    length = data_int.shape[0]
    
    # 分析各个方向的孔径
    avarage_x_gol = 0
    avarage_y_gol = 0
    avarage_z_gol = 0
    zero_in_vol = 0
    
    for i in range(length):
        # X方向
        yz_data = data_int[i, :, :]
        zero_in_vol_x, _ = analyze_plane(yz_data)
        avarage_x = calculate_pore_size(yz_data)
        avarage_x_gol += avarage_x
        
        # Y方向
        xz_data = data_int[:, i, :]
        zero_in_vol_y, _ = analyze_plane(xz_data)
        avarage_y = calculate_pore_size(xz_data)
        avarage_y_gol += avarage_y
        
        # Z方向
        xy_data = data_int[:, :, i]
        zero_in_vol_z, _ = analyze_plane(xy_data)
        avarage_z = calculate_pore_size(xy_data)
        avarage_z_gol += avarage_z
        
        zero_in_vol += (zero_in_vol_x + zero_in_vol_y + zero_in_vol_z) / 3
    
    # 计算平均值
    avarage_all = (avarage_x_gol + avarage_y_gol + avarage_z_gol) / (3 * length)
    avarage_all_x = avarage_x_gol / length
    avarage_all_y = avarage_y_gol / length
    avarage_all_z = avarage_z_gol / length
    zero_in_vol_ava = zero_in_vol / length
    
    # 结果字典
    results = {
        'average_pore_size': avarage_all * resolution,  # 转换为实际尺寸
        'average_pore_size_x': avarage_all_x * resolution,
        'average_pore_size_y': avarage_all_y * resolution,
        'average_pore_size_z': avarage_all_z * resolution,
        'porosity': zero_in_vol_ava
    }
    
    logging.info(f'平均孔径计算完成: {results["average_pore_size"]:.3f} nm')
    return results


def calculate_diffusion_coefficient(data_file, trajectory_file, selection='all', start_frame=0, end_frame=1000, timestep=1):
    """计算扩散系数
    
    Args:
        data_file: 数据文件路径
        trajectory_file: 轨迹文件路径
        selection: 原子选择字符串，默认为'all'
        start_frame: 起始帧，默认为0
        end_frame: 结束帧，默认为1000
        timestep: 时间步长，默认为1
        
    Returns:
        dict: 包含扩散系数和误差的字典
    """
    logging.info(f'开始计算扩散系数: {trajectory_file}')
    
    # 加载动力学帧
    u = mda.Universe(data_file, trajectory_file, format="LAMMPS")
    
    # 执行MSD分析
    msd_analysis = msd.EinsteinMSD(u, select=selection, msd_type='xyz', fft=True)
    msd_analysis.run(start=start_frame, stop=end_frame)
    
    # 获取MSD时间序列
    msd_values = msd_analysis.results.timeseries
    
    # 获取帧数并设置时间步长
    nframes = msd_analysis.n_frames
    start_index = int(start_frame / timestep)
    end_index = min(int(end_frame / timestep), nframes)

    lagtimes = np.arange(nframes) * timestep  # 创建时间轴

    # 执行线性回归
    linear_model = linregress(lagtimes[start_index:end_index], msd_values[start_index:end_index])
    slope = linear_model.slope
    error = linear_model.stderr
    
    # 计算扩散系数
    D = slope * 1 / (2 * msd_analysis.dim_fac)
    
    # 结果字典
    results = {
        'diffusion_coefficient': D,
        'error': error,
        'msd_values': msd_values,
        'lagtimes': lagtimes,
        'linear_model': linear_model
    }
    
    logging.info(f'扩散系数计算完成: {D:.6f}')
    return results


def plot_msd_curve(results, output_path=None):
    """绘制MSD曲线
    
    Args:
        results: 扩散系数计算结果
        output_path: 输出路径，默认为None
        
    Returns:
        None
    """
    lagtimes = results['lagtimes']
    msd_values = results['msd_values']
    linear_model = results['linear_model']
    start_index = 1  # 从第二个点开始拟合
    end_index = len(lagtimes)
    
    plt.figure(figsize=(10, 6))
    plt.plot(lagtimes, msd_values, label='MSD')
    plt.plot(lagtimes[start_index:end_index], 
             linear_model.intercept + linear_model.slope * lagtimes[start_index:end_index], 
             label='Linear Fit', linestyle='--')
    plt.xlabel(r'$\tau$ (ps)', fontsize=12)
    plt.ylabel('MSD (Å²)', fontsize=12)
    plt.title('Mean Square Displacement', fontsize=14)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.7)
    
    if output_path:
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        logging.info(f'MSD曲线已保存至: {output_path}')
    else:
        plt.show()
    
    plt.close()


def save_results_to_csv(results, output_path):
    """将结果保存为CSV文件
    
    Args:
        results: 结果字典
        output_path: 输出路径
        
    Returns:
        None
    """
    df = pd.DataFrame([results])
    df.to_csv(output_path, index=False)
    logging.info(f'结果已保存至: {output_path}')


def analyze_h2o_diffusion(data_file, trajectory_file, output_folder='分析结果'):
    """分析H2O的扩散系数
    
    Args:
        data_file: 数据文件路径
        trajectory_file: 轨迹文件路径
        output_folder: 输出文件夹，默认为'分析结果'
        
    Returns:
        dict: 扩散系数计算结果
    """
    # 创建输出文件夹
    os.makedirs(output_folder, exist_ok=True)
    
    # 计算扩散系数
    results = calculate_diffusion_coefficient(data_file, trajectory_file)
    
    # 绘制MSD曲线
    plot_path = os.path.join(output_folder, 'H2O_msd_plot.png')
    plot_msd_curve(results, plot_path)
    
    # 保存MSD数据
    msd_data = pd.DataFrame({
        't': results['lagtimes'], 
        'MSD': results['msd_values']
    })
    msd_csv_path = os.path.join(output_folder, 'H2O_msd.csv')
    msd_data.to_csv(msd_csv_path, index=False)
    
    # 保存扩散系数结果
    results_dict = {
        'dataset': 'H2O',
        'D': results['diffusion_coefficient'],
        'error': results['error']
    }
    results_csv_path = os.path.join(output_folder, 'H2O_扩散系数.csv')
    save_results_to_csv(results_dict, results_csv_path)
    
    return results


def analyze_pa_pore_size(xyz_file_path, atom_types=['C', 'N', 'O'], cell_length=100, resolution=5, output_folder='分析结果'):
    """分析PA的平均孔径
    
    Args:
        xyz_file_path: XYZ文件路径
        atom_types: 要分析的原子类型列表，默认为['C', 'N', 'O']
        cell_length: 晶格长度，默认为100
        resolution: 分辨率，默认为5
        output_folder: 输出文件夹，默认为'分析结果'
        
    Returns:
        dict: 孔径计算结果
    """
    # 创建输出文件夹
    os.makedirs(output_folder, exist_ok=True)
    
    # 计算平均孔径
    results = calculate_average_pore_size(xyz_file_path, atom_types, cell_length, resolution)
    
    # 保存结果
    results_csv_path = os.path.join(output_folder, 'PA_孔径分析.csv')
    save_results_to_csv(results, results_csv_path)
    
    # 生成结果文本
    result_text = f"""体平均孔径: {results['average_pore_size']:.3f} nm
"""
    result_text += f"""X方向平均孔径: {results['average_pore_size_x']:.3f} nm
"""
    result_text += f"""Y方向平均孔径: {results['average_pore_size_y']:.3f} nm
"""
    result_text += f"""Z方向平均孔径: {results['average_pore_size_z']:.3f} nm
"""
    result_text += f"""体平均孔隙率: {results['porosity']:.3f}"""
    
    # 保存结果文本
    result_txt_path = os.path.join(output_folder, 'PA_孔径分析.txt')
    with open(result_txt_path, 'w', encoding='utf-8') as f:
        f.write(result_text)
    
    logging.info(f'PA孔径分析结果已保存至: {output_folder}')
    return results


if __name__ == "__main__":
    # 设置输出文件夹
    output_folder = '分析结果'
    os.makedirs(output_folder, exist_ok=True)
    
    # 从template.py中获取文件路径
    from src.template import generate_simulation_params
    params = generate_simulation_params()
    h2o_dcd_path = params['H2O_dump']
    pa_xyz_path = params['final_PA_xyz']
    
    # 分析H2O的扩散系数
    data_file = 'path/to/data/file.data'  # 需要替换为实际的数据文件路径
    h2o_results = analyze_h2o_diffusion(data_file, h2o_dcd_path, output_folder)
    
    # 分析PA的平均孔径
    pa_results = analyze_pa_pore_size(pa_xyz_path, output_folder=output_folder)
    
    logging.info('分析完成')