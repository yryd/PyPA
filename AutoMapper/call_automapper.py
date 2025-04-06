from typing import List, Optional
from AutoMapper.LammpsUnifiedCleaner import file_unifier
from AutoMapper.LammpsToMolecule import lammps_to_molecule
from AutoMapper.MapProcessor import map_processor

def run_automapper_clean(
    directory: str,
    data_files: List[str],
    coeff_file: str
) -> None:
    """
    调用 AutoMapper 的 clean 功能
    
    参数:
        directory: 工作目录路径
        data_files: 要处理的数据文件列表
        coeff_file: 系数文件路径
    """
    print(f'数据文件列表: {data_files}')
    file_unifier(directory, coeff_file, data_files)

def run_automapper_molecule(
    directory: str,
    data_file: str,
    save_name: str
) -> None:
    """
    调用 AutoMapper 的 molecule 功能
    
    参数:
        directory: 工作目录路径
        data_file: 输入数据文件路径
        save_name: 输出分子文件名
    """
    lammps_to_molecule(directory, data_file, save_name)

def run_automapper_map(
    directory: str,
    pre_reaction_file: str,
    post_reaction_file: str,
    pre_save_name: str,
    post_save_name: str,
    bonding_atoms: List[str],
    elements_by_type: List[str],
    delete_atoms: Optional[List[str]] = None,
    create_atoms: Optional[List[str]] = None,
    debug: bool = False
) -> None:
    """
    调用 AutoMapper 的 map 功能
    
    参数:
        directory: 工作目录路径
        pre_reaction_file: 反应前数据文件路径
        post_reaction_file: 反应后数据文件路径
        pre_save_name: 反应前分子保存名称
        post_save_name: 反应后分子保存名称
        bonding_atoms: 成键原子ID列表 (至少4个)
        elements_by_type: 按类型的元素符号列表
        delete_atoms: 要删除的原子ID列表 (可选)
        create_atoms: 要创建的原子ID列表 (可选)
        debug: 是否启用调试模式 (默认为False)
    """
    # 分割成键原子为前两个和后两个
    pre_bonding_atoms = bonding_atoms[:2]
    post_bonding_atoms = bonding_atoms[2:]
    
    map_processor(
        directory, 
        pre_reaction_file, 
        post_reaction_file, 
        pre_save_name, 
        post_save_name, 
        pre_bonding_atoms, 
        post_bonding_atoms, 
        delete_atoms, 
        elements_by_type, 
        create_atoms, 
        debug
    )

# 示例用法
if __name__ == "__main__":
    # 示例1: 调用 clean 功能
    run_automapper_clean(
        directory=".",
        data_files=["pre-reaction.data", "post-reaction.data"],
        coeff_file="system.in.settings"
    )
    
    # 示例2: 调用 molecule 功能
    run_automapper_molecule(
        directory=".",
        data_file="cleanedpre-reaction.data",
        save_name="pre-molecule.data"
    )
    
    # 示例3: 调用 map 功能
    run_automapper_map(
        directory=r"PyPA/tmp/6/map",
        pre_reaction_file="cleanedpre.data",
        post_reaction_file="cleanedpost_0.data",
        pre_save_name="pre_mol.data",
        post_save_name="post_0_mol.data",
        bonding_atoms=["3", "18", "4", "5"],
        elements_by_type=["C", "C", "H", "H", "H", "H", "Cl", "N", "N", "O"],
        delete_atoms=["19", "11", "39", "40"],
        debug=True
    )