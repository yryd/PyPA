# execute.py

import logging
import os
import glob
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


def exec_ltemplify(data_PATH, lt_PATH, name_list):
    """根据分子 DATA 文件生成分子 LT 分子模板文件并指认 GAFF2 力场。
    
    Args:
        data_PATH: 分子 DATA 文件的路径。
        lt_PATH: 输出 LAMMPS 模板文件的路径。
        name: 分子名称。
    """
    data_str = '.data'
    lt_str = '.lt'
    str_list = []
    
    for name in name_list:
        str_list.append(f'ltemplify.py -name "{name} inherits GAFF2" -molid "1"  -ignore-coeffs -ignore-angles -ignore-bond-types -ignore-masses {data_PATH}{name}{data_str} > {lt_PATH}{name}{lt_str}\n')
    # 将命令写入到脚本文件
    with open(lt_PATH + 'ltemplify_my.sh', 'w') as file:
        file.writelines(str_list)

    try:
        os.system(f'. {lt_PATH}ltemplify_my.sh > /dev/null 2>&1')
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

def exec_moltemplate(lt_PATH, sys_PATH, data_PATH, loop = 1):
    """使用 moltemplate 通过 LT 文件，构建 LAMMPS DATA 文件。
    
    Args:
        lt_PATH: LAMMPS 模板文件的路径。
        sys_PATH: 系统文件的路径。
        data_PATH: 输出数据文件的路径。
    """
    # 使用 str_list 列表构建命令内容
    str_list = []
    str_list.append(f'moltemplate.sh -atomstyle "full" {lt_PATH}pre.lt\n')
    for index in range(loop):
        str_list.append(f'moltemplate.sh -atomstyle "full" {lt_PATH}post_{index}.lt\n')
    str_list.append(f'moltemplate.sh -pdb {sys_PATH}system.pdb -atomstyle "full" {lt_PATH}system.lt\n')

    # 将命令写入到脚本文件
    with open(data_PATH + 'run_moltemplate.sh', 'w') as file:
        file.writelines(str_list)

    tmp_PATH = os.getcwd()
    os.chdir(data_PATH)

    # 执行脚本并丢弃输出信息
    try:
        logging.info("等待构建反应模板 DATA 文件...")
        os.system(f'. {data_PATH}run_moltemplate.sh > /dev/null 2>&1')
    except Exception as e:
        logging.error(f"构建反应模板出错: {e}")
        raise

    # 删除不需要的文件
    os.remove('pre.in')
    os.remove('pre.in.init')
    os.remove('pre.in.settings')
    # os.remove('run.in.EXAMPLE')

    for index in range(loop):
        os.remove(f'post_{index}.in')
        os.remove(f'post_{index}.in.init')
        os.remove(f'post_{index}.in.settings')

    os.chdir(tmp_PATH)
    os.system(f'rm -rf {data_PATH}output_ttree')
    return

def exec_AutoMapper_clean(data_PATH, map_PATH):
    str_list = []
    str_list2 = []
    
    # 获取匹配的 post 文件，使用 glob 搜索文件
    post_files = glob.glob(os.path.join(data_PATH, 'post_*.data'))
    post_files_str = ' '.join(os.path.basename(file) for file in post_files)
    str_list.append(f'AutoMapper.py {data_PATH} clean pre.data {post_files_str} system.data --coeff_file system.in.settings\n')

    with open(data_PATH + 'AutoMapper_clean.sh', 'w') as file:
        file.writelines(str_list)
    
    tmp_PATH = os.getcwd()
    os.chdir(data_PATH)

    # 执行脚本并丢弃输出信息
    try:
        logging.info("清理反应模板 DATA 文件...")
        os.system(f'. {data_PATH}AutoMapper_clean.sh')
        logging.info("整理力场信息完成")
    except Exception as e:
        logging.error(f"清理文件力场出错: {e}")
        raise
    finally:
        cleaned_files = glob.glob(os.path.join(data_PATH, 'cleaned*.data'))
        for cleaned_file in cleaned_files:
            str_list2.append(f'cp {cleaned_file} {map_PATH}\n')

        with open(data_PATH + 'cp_clean_file.sh', 'w') as file:
            file.writelines(str_list2)
        
        os.system(f'. {data_PATH}cp_clean_file.sh > /dev/null 2>&1')

        os.chdir(tmp_PATH)


def exec_AutoMapper(map_PATH, reaction_index_dicts_list, atom_ele_type_str):
    """生成映射模板文件。
    
    Args:
        map_PATH: 映射文件的路径。
        data_PATH: 数据文件的路径。
    """

    str_list = []
    for index, reaction_index_dict in enumerate(reaction_index_dicts_list):
        # 反应前后的连接原子
        bond_str = '--ba ' + str(reaction_index_dict['N_r']) + ' ' + str(reaction_index_dict['C_r']) + ' ' + str(reaction_index_dict['N_p']) + ' ' + str(reaction_index_dict['C_p']) + ' '
        del_str = '--da '  + str(reaction_index_dict['Cl_d']) + ' ' + str(reaction_index_dict['H_d']) + ' ' + str(reaction_index_dict['Cl_p']) + ' ' + str(reaction_index_dict['H_p']) + ' '
        elem_str = '--ebt ' + atom_ele_type_str  + ' '
        str_list.append(f'AutoMapper.py {map_PATH} map cleanedpre.data cleanedpost_{index}.data --save_name pre_mol.data post_{index}_mol.data ' + bond_str + del_str + elem_str + '\n')
        str_list.append(f'mv {map_PATH}automap.data {map_PATH}automap_{index}.data' + '\n')
    
    with open(map_PATH + 'AutoMapper.sh', 'w') as file:
        file.writelines(str_list)
    
    tmp_PATH = os.getcwd()
    os.chdir(map_PATH)
    try:
        logging.info("自动映射，产生反应模板文件中...")
        os.system(f'. {map_PATH}AutoMapper.sh')
    except Exception as e:
        logging.error(f"自动映射出错: {e}")
        raise
    finally:
        os.chdir(tmp_PATH)
