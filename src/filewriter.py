# filewriter.py

def write_file(file_PATH, file_name, str_list):
    with open(file_PATH + file_name, 'w') as file:
        file.writelines(str_list)

def write_sys_lt_str(lt_PATH, name_list, num_list, box_len = 0):
    """生成 LT 系统模板文件的字符串内容。
    
    Args:
        name_list: 模板分子的名称列表。
    
    Returns:
        str_list: 文件的多行字符串列表。
    """
    str_list = []
    str_list.append('import "gaff2.lt"\n')
    for i in name_list:
        str_list.append(f'import "{lt_PATH}{i}.lt"\n')
    str_list.append('\n')
    for i,j in zip(name_list, num_list):
        str_list.append(f'mol_{i} = new {i} [{j}]\n')
    
    if box_len:
        str_list.append('write_once("Data Boundary") {\n')
        str_list.append(f'    0.0 {box_len} xlo xhi\n')
        str_list.append(f'    0.0 {box_len} xlo xhi\n')
        str_list.append(f'    0.0 {box_len} zlo zhi\n')
        str_list.append('}\n')
    return str_list

def write_sys_lt(lt_PATH, name_list, num_list, box_len = 0):
    """生成 LT 系统模板文件。
    
    Args:
        name_list: 模板分子的名称列表。
    
    Returns:
        str_list: 文件的多行字符串列表。
    """
    str_list = write_sys_lt_str(lt_PATH, name_list, num_list, box_len)
    with open(lt_PATH + 'system.lt', 'w') as file:
        file.writelines(str_list)
    return

def write_template_lt(lt_PATH, pre_template_list, post_template_list, index = 0):
    """根据给定的模板列表生成 LAMMPS 模板文件。
    
    Args:
        lt_PATH: 输出 LAMMPS 模板文件的路径。
        pre_template_list: 反应物模板分子名称列表。
        post_template_list: 产物模板分子名称列表。
    """
    pre_str = write_sys_lt_str(lt_PATH, pre_template_list, [1,1])
    post_str = write_sys_lt_str(lt_PATH, post_template_list, [1,1])
    write_file(lt_PATH, 'pre.lt', pre_str)
    write_file(lt_PATH, f'post_{index}.lt', post_str)

def write_packmol_inp(sys_PATH, mol_PATH, mol_list, mol_num_list, box_len, file_type = 'pdb', output_name = 'system'):
    """生成 Packmol 输入文件。
    
    Args:
        sys_PATH: Packmol 输入文件的路径。
        mol_PATH: 分子文件的路径。
        mol_list: 分子名称列表。
        mol_num_list: 每种分子的数量列表。
        box_len: 系统的盒子长度。
    """
    str_list = []
    str_list.append('nloop0 1000\n')
    str_list.append('tolerance 2.0\n')
    str_list.append(f'filetype {file_type}\n')
    str_list.append(f'output {sys_PATH}{output_name}.{file_type}\n')
    str_list.append(f'add_box_sides {box_len}\n')
    
    for i,j in zip(mol_list, mol_num_list):
        # 数量必须大于 0 才不报错
        if j:
            str_list.append(f'structure {mol_PATH}{i}.{file_type}\n')
            str_list.append(f'  number {j}\n')
            str_list.append(f'  inside box 0. 0. 0. {box_len}. {box_len}. {box_len}.\n')
            str_list.append(f'end structure\n')
    
    write_file(sys_PATH, 'system.inp', str_list)
    return


def write_lammps_in(data_PATH, map_PATH, run_PATH):
    """生成 LAMMPS 输入文件以进行模拟。
    
    Args:
        data_PATH: 数据文件的路径。
        map_PATH: 映射文件的路径。
        run_PATH: LAMMPS 输入文件的输出路径。
    """
    # 使用 str_list 列表构建输入文件内容
    str_list = []
    str_list.append("# ----------------- Init Section -----------------\n")
    str_list.append("units           real\n")
    str_list.append("atom_style      full\n")
    str_list.append("bond_style      hybrid harmonic\n")
    str_list.append("angle_style     hybrid harmonic\n")
    str_list.append("dihedral_style  hybrid fourier\n")
    str_list.append("improper_style  hybrid cvff\n")
    str_list.append("pair_style      hybrid lj/charmm/coul/long 9.0 10.0 10.0\n")
    str_list.append("kspace_style    pppm 0.0001\n\n")

    str_list.append("pair_modify     mix arithmetic\n")
    str_list.append("special_bonds   amber\n")
    str_list.append("# ----------------- Atom Definition Section -----------------\n")
    str_list.append(f'read_data "{data_PATH}cleanedsystem.data"\n')
    str_list.append("# ----------------- Settings Section -----------------\n")
    str_list.append(f'include "{data_PATH}cleanedsystem.in.settings"\n')
    str_list.append("# ----------------- Simulation -----------------\n\n")

    str_list.append(f'molecule        pre {map_PATH}pre_mol.data\n')
    str_list.append(f'molecule        post {map_PATH}post_mol.data\n\n')

    str_list.append("neighbor        2.5 bin\n")
    str_list.append("neigh_modify    every 1 delay 0 check yes\n\n")

    str_list.append("group           MPD molecule <> 1 100 \n")
    str_list.append("group           TMC molecule <> 101 250 \n\n")

    str_list.append("velocity        all create 300.0 4928459 rot yes dist gaussian \n\n")

    str_list.append("min_style       cg\n")
    str_list.append("minimize        1e-05 1e-05 10000 100000\n\n")

    str_list.append("timestep        1.0\n")
    str_list.append("fix             relax_md all nvt temp 298.15 298.15 100.0 \n")
    str_list.append("run             3000\n")
    str_list.append("unfix relax_md\n\n")

    str_list.append("dump            mydump1 all xtc 100 traj_npt.xtc\n\n")

    str_list.append("timestep        1.0\n")
    str_list.append("fix             npt_md all npt temp 298.15 298.15 100.0 iso 1.0 1.0 1000.0 \n")
    str_list.append("run             5000\n")
    str_list.append("unfix npt_md\n\n")

    str_list.append("fix             xlink_fix all bond/react stabilization yes statted_grp .03 &\n")
    str_list.append("                    react rxn1 all 100 0.0 10.0 pre post {map_PATH}automap.data stabilize_steps 100\n")
    str_list.append("fix             nvt_md all nvt temp 298.15 298.15 100.0\n")
    str_list.append("run             100000\n\n")

    str_list.append("thermo          100\n\n")

    str_list.append("write_restart   system_after_npt.rst\n")
    str_list.append("write_data      system.data\n\n")

    str_list.append("# write_dump      all custom pysimm.dump.tmp id q x y z vx vy vz\n")
    str_list.append("quit\n")

    # 将构建的内容写入 LAMMPS 输入文件
    write_file(run_PATH, 'in.system', str_list)