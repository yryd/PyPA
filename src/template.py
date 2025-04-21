# template.py
# 为lammps生成模板文件

import logging
import jinja2
import os

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
def generate_xlink_lammps_params(path_dict):
    """
    根据输入字典生成LAMMPS模板参数
    
    Args:
        path_dict (dict): 输入配置字典
        
    Returns:
        dict: 包含所有LAMMPS模板参数的字典
    """
    xlink_PATH = path_dict['paths']['xlink']

    
    # 根据map_templates动态生成分子模板和反应参数
    molecule_templates = {}
    reactions_list = []
    
    # 遍历map_templates生成分子模板和反应参数
    for i, map_key in enumerate(path_dict['map_templates'].keys(), 1):
        # 添加分子模板
        pre_mol_name = f'pre{i}'
        post_mol_name = f'post{i}'
        
        molecule_templates[pre_mol_name] = xlink_PATH + f'mol.pre_{map_key}'
        molecule_templates[post_mol_name] = xlink_PATH + f'mol.post_{map_key}'
        
        # 添加反应参数
        reactions_list.append({
            'name': f'rxn{i}',
            'nevery': 100,
            'mincutoff': 0.0,
            'maxcutoff': 5.0,
            'pre_mol': pre_mol_name,
            'post_mol': post_mol_name,
            'map_file': xlink_PATH + f'txt.{map_key}',
            'stabilize_steps': 60
        })
    
    # 获取元素列表
    elements = path_dict['elem_type_str']
    
    # 获取分子数量范围
    molecule_numbers = {}
    r1_count = 0
    r2_count = 0
    
    # 从path_dict中获取反应物数量
    if 'num' in path_dict:
        for mol_name, count in path_dict['num'].items():
            if mol_name.startswith('r1'):
                r1_count += count
            elif mol_name.startswith('r2'):
                r2_count += count
    
    if r1_count > 0:
        molecule_numbers['N_single'] = [1, r1_count]
    else:
        molecule_numbers['N_single'] = [1, 150]
        
    if r2_count > 0:
        molecule_numbers['Cl_single'] = [r1_count + 1, r1_count + r2_count]
    else:
        molecule_numbers['Cl_single'] = [151, 250]
    
    # 定义模板参数
    params = {
        # 基本文件设置
        'x_file_name': xlink_PATH + 'in.xlink.lammps',
        'input_file': xlink_PATH + 'sys_init.lmps',
        'log_file': xlink_PATH + 'xlink_simulation.log',
        
        # 分子模板设置
        'molecule_templates': molecule_templates,
        
        # 分子数量范围
        'molecule_numbers': molecule_numbers,
        
        # 元素符号列表
        'elements': elements,
        
        # 模拟条件
        'temperature': 298.15,
        'pressure': 1.0,
        'timestep': 1.0,
        # 运行时间设置
        'run_time': {
            'equilibration': 50000,
            'production': 100000
        },
        
        # 输出文件设置
        'output_files': {
            'balance_log': xlink_PATH + 'balance_xlink.log',
            'equilibration': xlink_PATH + 'system_after_bal_npt.xyz',
            'dump_xyz': xlink_PATH + 'dump_xlink.xyz',
            'dump_custom': xlink_PATH + 'dump_xlink.custom',
            'react_log': xlink_PATH + 'steps_react.log',
            'final_data': xlink_PATH + 'system_after_xlink.lmp',
            'final_xyz': xlink_PATH + 'system_after_xlink_nosol.xyz',
            'final_nosol_data': xlink_PATH + 'system_after_xlink_nosol.lmp',
            'final_restart': xlink_PATH + 'system_after_xlink_nosol.rst'
        },
        
        # 反应参数设置
        'reactions': {
            'noact_grp': 0.03,
            'reactions': reactions_list
        }
    }
    
    return params

def generate_insertH2O_lammps_params(path_dict):
    """
    根据输入字典生成LAMMPS Na2SO4水溶液模板参数
    
    Args:
        path_dict (dict): 输入配置字典
        
    Returns:
        dict: 包含所有LAMMPS模板参数的字典
    """
    insert_PATH = path_dict['paths']['insert']
    xlink_PATH = path_dict['paths']['xlink']
    tplt_PATH = path_dict['paths']['template']
    
    # 获取元素字符串
    elements = path_dict['elem_type_str']
    
    # 定义模板参数
    params = {
        # 基本文件设置
        'h2o_file_name': insert_PATH + 'in.insertH2O.lammps',
        'input_file': xlink_PATH + 'system_after_xlink_nosol.lmp',
        'log_file': insert_PATH + 'h2o_simulation.log',
        
        # 分子数据文件
        'Na_data_file': tplt_PATH + 'Na.data',
        'SO4_data_file': tplt_PATH + 'SO4.data',
        'H2O_data_file': tplt_PATH + 'OPC3.data',
        
        # 分子模板文件
        'Na_mol_file': insert_PATH + 'Na.txt',
        'SO4_mol_file': insert_PATH + 'SO4.txt',
        'H2O_mol_file': insert_PATH + 'OPC3.txt',
        
        # 插入分子数量
        'Na_count': 20,
        'SO4_count': 10,
        'H2O_count': 600,
        
        # 分子重叠参数
        'Na_overlap': 2,
        'SO4_overlap': 2,
        'H2O_overlap': 2,
        
        # 元素符号列表
        'elements': elements + ' Na S O O H',
        
        # 模拟条件
        'temperature': 298.15,
        'pressure': 1.0,
        'timestep': 1.0,
        
        # 运行时间设置
        'equilibration_steps': 20000,
        'production_steps': 500000,
        
        # 输出文件设置
        'balance_log': insert_PATH + 'balance_H2O.log',
        'dump_custom': insert_PATH + 'Na2SO4_balance.custom',
        'H2O_dump': insert_PATH + 'H2O_dcd.dcd',
        'Na2SO4_dump': insert_PATH + 'Na2SO4_dcd.dcd',
        'final_data': insert_PATH + 'Na2SO4_balance.lmp',
        'final_xyz': insert_PATH + 'Na2SO4_balance.xyz',
        'final_PA_xyz': insert_PATH + 'PA.xyz'
    }
    
    return params

def generate_mol_lammps_params(path_dict):
    """
    根据输入字典生成OPC3水分子和SO4分子模板参数
    
    Args:
        path_dict (dict): 输入配置字典
        
    Returns:
        tuple: 包含OPC3水分子和SO4分子模板参数的字典元组
    """
    insert_PATH = path_dict['paths']['insert']
    angle_types = path_dict['sys_info']['angle_types']
    bond_types = path_dict['sys_info']['bond_types']
    
    # 定义OPC3水分子模板参数
    opc3_params = {
        'opc3_file_name': insert_PATH + 'OPC3.txt',
        'H2O_angle_types': angle_types + 2,
        'H2O_bond_types': bond_types + 2
    }
    
    # 定义SO4分子模板参数
    so4_params = {
        'so4_file_name': insert_PATH + 'SO4.txt',
        'SO4_angle_types': angle_types + 1,
        'SO4_bond_types': bond_types + 1
    }
    
    return opc3_params, so4_params

def add_all_lammps_params(path_dict):
    """
    生成所有LAMMPS模拟所需的参数并添加到path_dict中
    
    Args:
        path_dict (dict): 输入配置字典
        
    Returns:
        dict: 更新后的配置字典
    """
    
    # 生成各种参数
    xlink_params = generate_xlink_lammps_params(path_dict)
    insert_params = generate_insertH2O_lammps_params(path_dict)
    opc3_params, so4_params = generate_mol_lammps_params(path_dict)
    
    # 将参数添加到path_dict中
    path_dict['params'] = {
        'xlink': xlink_params,
        'insert': insert_params,
        'opc3': opc3_params,
        'so4': so4_params
    }
    
    return path_dict

def write_lammps_template(path_dict):
    """
    根据path_dict中的参数生成LAMMPS模板文件

    Args:
        path_dict (dict): 包含所有参数的字典
    """
    # 获取模板参数
    xlink_params = path_dict['params']['xlink']
    insert_params = path_dict['params']['insert']
    opc3_params = path_dict['params']['opc3']
    so4_params = path_dict['params']['so4']
    tplt_PATH = path_dict['paths']['template']
    insert_PATH = path_dict['paths']['insert']

    # 创建Jinja2环境
    env = jinja2.Environment(loader=jinja2.FileSystemLoader(tplt_PATH), trim_blocks=True, lstrip_blocks=True)
    
    # 加载并渲染xlink模板
    template_xlink = env.get_template('in.xlink.template.lammps')
    rendered_xlink = template_xlink.render(**xlink_params)
    # 写入输出文件
    with open(xlink_params['x_file_name'], 'w') as f:
        f.write(rendered_xlink)
    logging.info(f"已生成LAMMPS输入文件: {xlink_params['x_file_name']}")
    
    # 加载并渲染OPC3水分子模板
    template_opc3 = env.get_template('OPC3.template.txt')
    rendered_opc3 = template_opc3.render(**opc3_params)
    # 写入输出文件
    with open(opc3_params['opc3_file_name'], 'w') as f:
        f.write(rendered_opc3)
    logging.info(f"已生成分子模板文件: {opc3_params['opc3_file_name']}")
    
    # 加载并渲染SO4分子模板
    template_so4 = env.get_template('SO4.template.txt')
    rendered_so4 = template_so4.render(**so4_params)
    # 写入输出文件
    with open(so4_params['so4_file_name'], 'w') as f:
        f.write(rendered_so4)
    logging.info(f"已生成分子模板文件: {so4_params['so4_file_name']}")

    # 复制Na分子文件
    os.system(f'cp {tplt_PATH}Na.txt {insert_PATH}')

    # 加载并渲染insertH2O模板
    template_insert = env.get_template('in.insertH2O.template.lammps')
    rendered_insert = template_insert.render(**insert_params)
    # 写入输出文件
    with open(insert_params['h2o_file_name'], 'w') as f:
        f.write(rendered_insert)
    logging.info(f"已生成LAMMPS输入文件: {insert_params['h2o_file_name']}")