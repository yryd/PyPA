# execute.py
from AutoMapper.call_automapper import run_automapper_clean, run_automapper_map

import logging
import os
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def exec_ltemplify(data_PATH, lt_PATH, name):
    """根据分子 DATA 文件生成分子 LT 分子模板文件并指认 GAFF2 力场。
    
    Args:
        data_PATH: 分子 DATA 文件的路径。
        lt_PATH: 输出 LAMMPS 模板文件的路径。
        name: 分子名称。
    """
    data_str = '.data'
    lt_str = '.lt'
    str_list = []
    
    str_list.append(f'ltemplify.py -name "{name} inherits GAFF2" -molid "1"  -ignore-coeffs -ignore-angles -ignore-bond-types -ignore-masses {data_PATH}{name}{data_str} > {lt_PATH}{name}{lt_str}\n')
    # 将命令写入到脚本文件
    with open(lt_PATH + 'ltemplify_my.sh', 'w') as file:
        file.writelines(str_list)

    try:
        os.system(f'. {lt_PATH}ltemplify_my.sh > /dev/null 2>&1')
        logging.info(f"构建{name} lt文件完成")
    except Exception as e:
        logging.error(f"构建lt文件出错: {e}")
        raise
    # os.remove(f"{lt_PATH}{name}{data_str}")

def exec_packmol(sys_PATH):
    """运行 Packmol 以建立分子系统模型。
    
    Args:
        sys_PATH: Packmol 输入文件的路径。
    """
    with open(sys_PATH + 'run_packmol.sh', 'w') as file:
        file.write(f'packmol < {sys_PATH}system.inp')
    try:
        logging.info("等待packmol运行...")
        os.system(f'. {sys_PATH}run_packmol.sh > /dev/null 2>&1')
    except Exception as e:
        logging.error(f"packmol运行失败: {e}")
        raise
    return

def exec_moltemplate_system(lt_PATH, sys_PATH, data_PATH):
    """使用 moltemplate 通过 LT 文件，构建系统的 LAMMPS DATA 文件。
    
    Args:
        lt_PATH: LAMMPS 模板文件的路径。
        sys_PATH: 系统文件的路径。
        data_PATH: 输出数据文件的路径。
    """
    # 使用 str_list 列表构建命令内容
    str_list = []
    str_list.append(f'moltemplate.sh -pdb {sys_PATH}system.pdb -atomstyle "full" {lt_PATH}system.lt\n')

    # 将命令写入到脚本文件
    with open(data_PATH + 'run_moltemplate_system.sh', 'w') as file:
        file.writelines(str_list)

    tmp_PATH = os.getcwd()
    os.chdir(data_PATH)

    # 执行脚本并丢弃输出信息
    try:
        logging.info("等待构建系统 DATA 文件...")
        os.system(f'. {data_PATH}run_moltemplate_system.sh > /dev/null 2>&1')
    except Exception as e:
        logging.error(f"构建系统模板出错: {e}")
        raise

    os.chdir(tmp_PATH)
    os.system(f'rm -rf {data_PATH}output_ttree')
    return

def exec_moltemplate_reaction(lt_PATH, data_PATH, map_name):
    """使用 moltemplate 通过 LT 文件，构建反应前后的 LAMMPS DATA 文件。
    
    Args:
        lt_PATH: LAMMPS 模板文件的路径。
        data_PATH: 输出数据文件的路径。
        map_name: 反应模板名称。
    """
    # 使用 str_list 列表构建命令内容
    str_list = []
    str_list.append(f'moltemplate.sh -atomstyle "full" {lt_PATH}pre_{map_name}.lt\n')
    str_list.append(f'moltemplate.sh -atomstyle "full" {lt_PATH}post_{map_name}.lt\n')

    # 将命令写入到脚本文件
    with open(data_PATH + 'run_moltemplate_reaction.sh', 'w') as file:
        file.writelines(str_list)

    tmp_PATH = os.getcwd()
    os.chdir(data_PATH)

    # 执行脚本并丢弃输出信息
    try:
        logging.info("等待构建反应模板 DATA 文件...")
        os.system(f'. {data_PATH}run_moltemplate_reaction.sh > /dev/null 2>&1')
        logging.info(f"构建反应模板 {map_name} DATA 文件完成...")
    except Exception as e:
        logging.error(f"构建反应模板DATA出错: {e}")
        raise

    # 删除不需要的文件
    os.remove(f'pre_{map_name}.in')
    os.remove(f'pre_{map_name}.in.init')
    os.remove(f'pre_{map_name}.in.settings')

    os.remove(f'post_{map_name}.in')
    os.remove(f'post_{map_name}.in.init')
    os.remove(f'post_{map_name}.in.settings')

    os.chdir(tmp_PATH)
    return

def exec_AutoMapper_clean(data_PATH, map_PATH, map_name_list):
    """清理反应模板 DATA 文件并整理力场信息。
    
    Args:
        data_PATH: 数据文件的路径。
        map_PATH: 映射文件的路径。
    """
    data_files = []
    for map_name in map_name_list:
        data_files = data_files + [f'pre_{map_name}.data'] + [f'post_{map_name}.data']
    data_files = data_files + ['system.data']
    tmp_PATH = os.getcwd()
    os.chdir(data_PATH)

    # 直接调用Python函数而不是通过shell脚本
    try:
        logging.info("清理反应模板 DATA 文件...")
        run_automapper_clean(
            directory=".",  # 使用当前目录(已切换到data_PATH)
            data_files=data_files,
            coeff_file="system.in.settings"
        )
        logging.info("整理力场信息完成")
    except Exception as e:
        logging.error(f"清理文件力场出错: {e}")
        raise
    finally:
        # 复制清理后的文件到映射目录
        for map_name in map_name_list:
            os.system(f'cp "{data_PATH}cleanedpre_{map_name}.data" "{map_PATH}clean_pre_{map_name}.data"')
            os.system(f'cp "{data_PATH}cleanedpost_{map_name}.data" "{map_PATH}clean_post_{map_name}.data"')
        
        # 恢复原始工作目录
        os.chdir(tmp_PATH)


# def exec_AutoMapper2(map_PATH, reaction_index_dicts_list, atom_ele_type_str):
#     """生成映射模板文件。
    
#     Args:
#         map_PATH: 映射文件的路径。
#         data_PATH: 数据文件的路径。
#     """

#     str_list = []
#     for index, reaction_index_dict in enumerate(reaction_index_dicts_list):
#         # 反应前后的连接原子
#         bond_str = '--ba ' + str(reaction_index_dict['N_r']) + ' ' + str(reaction_index_dict['C_r']) + ' ' + str(reaction_index_dict['N_p']) + ' ' + str(reaction_index_dict['C_p']) + ' '
#         del_str = '--da '  + str(reaction_index_dict['Cl_d']) + ' ' + str(reaction_index_dict['H_d']) + ' ' + str(reaction_index_dict['Cl_p']) + ' ' + str(reaction_index_dict['H_p']) + ' '
#         elem_str = '--ebt ' + atom_ele_type_str  + ' '
#         str_list.append(f'AutoMapper.py {map_PATH} map cleanedpre.data cleanedpost_{index}.data --save_name pre_mol.data post_{index}_mol.data ' + bond_str + del_str + elem_str + '\n')
#         str_list.append(f'mv {map_PATH}automap.data {map_PATH}automap_{index}.data' + '\n')
    
#     with open(map_PATH + 'AutoMapper.sh', 'w') as file:
#         file.writelines(str_list)
    
#     tmp_PATH = os.getcwd()
#     os.chdir(map_PATH)
#     try:
#         logging.info("自动映射，产生反应模板文件中...")
#         os.system(f'. {map_PATH}AutoMapper.sh')
#     except Exception as e:
#         logging.error(f"自动映射出错: {e}")
#         raise
#     finally:
#         os.chdir(tmp_PATH)

def exec_AutoMapper(map_PATH, map_name_list, reaction_index_dicts_list, atom_ele_type_str):
    """生成映射模板文件。
    
    Args:
        map_PATH: 映射文件的路径。
        map_name_list: 反应模板名称列表。
        reaction_index_dicts_list: 反应索引字典列表。
        atom_ele_type_str: 原子元素类型字符串。
    """
    # 保存当前工作目录
    tmp_PATH = os.getcwd()
    
    try:
        # 切换到映射目录
        os.chdir(map_PATH)
        
        for reaction_index_dict, map_name in zip(reaction_index_dicts_list, map_name_list):
            try:
                logging.info(f"自动映射，处理反应 {map_name} 中...")
                run_automapper_map(
                    directory=".",  # 使用当前目录(已切换到map_PATH)
                    pre_reaction_file=f"clean_pre_{map_name}.data",  # 使用exec_AutoMapper_clean生成的文件
                    post_reaction_file=f"clean_post_{map_name}.data",  # 使用exec_AutoMapper_clean生成的文件
                    pre_save_name=f"mol.pre_{map_name}",
                    post_save_name=f"mol.post_{map_name}",
                    bonding_atoms=[
                        str(reaction_index_dict['N_r']),
                        str(reaction_index_dict['C_r']),
                        str(reaction_index_dict['N_p']),
                        str(reaction_index_dict['C_p'])
                    ],
                    elements_by_type=atom_ele_type_str.split(),
                    delete_atoms=[
                        str(reaction_index_dict['Cl_d']),
                        str(reaction_index_dict['H_d']),
                        str(reaction_index_dict['Cl_p']),
                        str(reaction_index_dict['H_p'])
                    ],
                    debug=False,
                    map_file_name=f"txt.{map_name}"
                )
                logging.info(f"处理反应 {map_name} 完成")
            except Exception as e:
                logging.error(f"处理反应 {map_name} 时出错: {e}")
                raise
                
    finally:
        # 确保无论是否出错都恢复原始工作目录
        os.chdir(tmp_PATH)