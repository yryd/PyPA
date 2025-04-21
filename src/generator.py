# generator.py
# 该模块处理反应物 SMILES 字符串生成
from src.molecular import MolecularModule
from src.optimizer import init_mol_prop
from rdkit import Chem
import logging
import os

# 设置日志
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def are_smiles_different(smiles1, smiles2):
    """
    检查两个 SMILES 是否代表不同的分子结构。
    
    Args:
        smiles1: 第一个 SMILES 表示法
        smiles2: 第二个 SMILES 表示法
    
    Returns:
        bool: 如果两个 SMILES 代表不同的分子结构，则返回 True，否则返回 False。
    """
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        raise ValueError("无效的 SMILES 表示法")

    # 生成标准化的 SMILES
    canonical_smiles1 = Chem.MolToSmiles(mol1, canonical=True)
    canonical_smiles2 = Chem.MolToSmiles(mol2, canonical=True)

    return canonical_smiles1 != canonical_smiles2

def generate_reaction_smile(r1, r2, g_type_smile = 'C(=O)N'):
    """
    生成反应物的 SMILES 字符串。
    
    Args:
        r1: 第一个反应物的 MolecularModule 对象
        r2: 第二个反应物的 MolecularModule 对象
        g_type_smile: 产物中要识别的官能团 SMILES 字符串，默认为酰胺基团
    
    Returns:
        tuple: 包含所有种产物 SMILES 列表、副产物 SMILES 和反应索引字典列表的元组
    """

    # print(r1.smiles)
    # print(r2.smiles)
    # 找到 r1 中所有可反应的胺基索引
    r1_n_indices = r1.mol_functional_group
    if not r1_n_indices:
        logging.warning(f"反应物1 {r1.smiles} 未找到可反应的胺基")
        return [], "[H]Cl", []
    
    # 找到 r2 中所有可反应的酰氯基团索引
    r2_c_indices = r2.mol_functional_group
    if not r2_c_indices:
        logging.warning(f"反应物2 {r2.smiles} 未找到可反应的酰氯基团")
        return [], "[H]Cl", []

    logging.info(f"找到反应物1的胺基数量: {len(r1_n_indices)}, 反应物2的酰氯基团数量: {len(r2_c_indices)}")
    
    product_smiles_list = []
    atom_indices_dict_list = []
    
    # 遍历 r1 中的每个 N 原子
    for n_index in r1_n_indices:
        # 遍历 r2 中的每个酰氯基团
        for c_index in r2_c_indices:
            atom_indices_dict = {}
            # 反应逻辑
            # 1. 创建新分子添加两个分子
            combined_mol = Chem.CombineMols(r1.mol, r2.mol)
            link_mol = Chem.RWMol(combined_mol)
            # 2. 找到 N 和 C 的索引
            n_atom_idx = n_index[0]  # N 的索引
            c_atom_idx = r1.mol.GetNumAtoms() + c_index[0]  # C 原子的索引（在 link_mol 中的位置）
            
            # 反应前需要键合的原子，注意原子序号从 1 开始，rdkit 索引从 0 开始
            atom_indices_dict['N_r'] = n_atom_idx + 1
            # 注意先后顺序，制作反应模板时也是未反应含两个反应物的体系
            atom_indices_dict['C_r'] = c_atom_idx + 1
            
            # 反应前标记需要删除的原子
            atom_indices_dict['Cl_d'] = r1.mol.GetNumAtoms() + c_index[2] + 1
            atom_indices_dict['H_d'] = n_index[1] + 1
            
            # 3. 形成酰胺键
            link_mol.AddBond(n_atom_idx, c_atom_idx, Chem.BondType.SINGLE)

            # 4. 去掉 Cl 原子，Cl 是在 r2 基团元组中的第三个元素
            cl_atom_idx = r1.mol.GetNumAtoms() + c_index[2] # Cl 原子的索引（在 link_mol 中的位置）
            # 从 link_mol 中去掉 Cl
            link_mol.RemoveAtom(cl_atom_idx)
            
            h_atom_idx = n_index[1]
            # 从 link_mol 中去掉一个 N 上的 H 原子
            link_mol.RemoveAtom(h_atom_idx)

            try:
                # 清除所有原子的原子映射号
                for atom in link_mol.GetAtoms():
                    atom.SetAtomMapNum(0)
                # 5. 生成新的 SMILES
                product_smiles = Chem.MolToSmiles(link_mol, canonical=True)
                if not product_smiles:
                    logging.error("生成产物 SMILES 失败")
                    continue

                product_mol = MolecularModule.mol_from_smiles(product_smiles)
                if not product_mol:
                    logging.error(f"从 SMILES {product_smiles} 创建分子对象失败")
                    continue
                
                # 重新计算产物基团索引
                p_group_matches = product_mol.GetSubstructMatches(Chem.MolFromSmiles(g_type_smile))
                if not p_group_matches:
                    logging.error(f"在产物 {product_smiles} 中未找到目标官能团")
                    continue
                p_group_index = list(p_group_matches)[0]
            except Exception as e:
                logging.error(f"处理产物时发生错误: {str(e)}")
                continue
            # 反应后键合的原子序号
            atom_indices_dict['N_p'] = p_group_index[2] + 1
            atom_indices_dict['C_p'] = p_group_index[0] + 1

            # 反应后删除的原子序号
            atom_indices_dict['Cl_p'] = product_mol.GetNumAtoms() + 1
            atom_indices_dict['H_p'] = product_mol.GetNumAtoms() + 2


            # 将产物 SMILES 添加到列表中（避免重复）
            if product_smiles not in product_smiles_list:
                product_smiles_list.append(product_smiles)
                atom_indices_dict_list.append(atom_indices_dict)
    
    # 副产物固定为HCl
    byproduct_smiles = "[H]Cl"

    return (product_smiles_list, byproduct_smiles, atom_indices_dict_list)


def init_product_info(path_dict, r1_mol, r2_mol, is_second_reaction=False):
    """
    初始化产物信息，将产物的各种信息填充到path_dict中
    
    Args:
        path_dict: 包含基本目录结构、反应信息的字典
        r1_mol: 第一个反应物的MolecularModule对象
        r2_mol: 第二个反应物的MolecularModule对象
        is_second_reaction: 是否为二次反应，默认为False
        
    Returns:
        path_dict: 更新后的字典，包含产物信息
        product_mols: 产物分子对象列表
    """
    logging.info('初始化产物分子信息...')
    
    # 生成反应产物SMILES
    product_smiles_list, byproduct_smile, reaction_index_dicts_list = generate_reaction_smile(r1_mol, r2_mol)
    
    product_mols = []
    # 创建产物文件名列表，格式为r1_r2_index
    for i in range(len(product_smiles_list)):
        # 生成当前产物的文件名
        if is_second_reaction:
            # 二次反应产物命名格式
            p_file_name = f"{path_dict['file_name'][0]}_{path_dict['file_name'][1]}_{i+1}_{path_dict['file_name'][1]}_2nd_{i+1}"
        else:
            p_file_name = f"{path_dict['file_name'][0]}_{path_dict['file_name'][1]}_{i+1}"
        
        # 更新path_dict中的file_name列表，添加产物
        path_dict['file_name'].append(p_file_name)
        
        # 更新type字典
        path_dict['type'][p_file_name] = 'p'
        
        # 更新smiles字典
        path_dict['smiles'][p_file_name] = product_smiles_list[i]
        
    # 添加反应索引信息
    if not is_second_reaction:
        path_dict['reaction_index_dicts_list'] = reaction_index_dicts_list
        path_dict['p_smiles_num'] = len(product_smiles_list)
    else:
        path_dict['reaction_index_dicts_list'] += reaction_index_dicts_list
        path_dict['p_smiles_num_2nd'] = len(product_smiles_list)
    # 处理副产物相关信息，只在第一次反应时添加
    if not is_second_reaction:
        byp_file_name = 'byp'
        path_dict['file_name'].append(byp_file_name)
        path_dict['type'][byp_file_name] = 'byp'
        path_dict['smiles'][byp_file_name] = byproduct_smile
    
    # 实例化产物分子
    if is_second_reaction:
        # 二次反应只处理产物，不处理副产物
        p_file_names = path_dict['file_name'][-len(product_smiles_list):]
        path_dict, p_mols = init_mol_prop(path_dict, p_file_names, mol_ff = None)
        product_mols.extend(p_mols)
    else:
        # 第一次反应处理产物和副产物
        p_file_names = path_dict['file_name'][3:-1] if 'byp' in path_dict['file_name'] else path_dict['file_name'][3:]
        byp_file_names = [path_dict['file_name'][-1]] if 'byp' in path_dict['file_name'] else []
        path_dict, p_mols = init_mol_prop(path_dict, p_file_names, mol_ff = None)
        if byp_file_names:
            path_dict, byp_mols = init_mol_prop(path_dict, byp_file_names, mol_ff = {'Cl': 'Cl', 'H': 'hx'})
            product_mols.extend(byp_mols)
        product_mols.extend(p_mols)
    
    logging.info('初始化产物分子信息完成')
    
    # 检查是否需要进行二次反应
    # 由于反应改变了N原子的原子类型，所有的键、角、二面角等信息都会改变，需要邻居原子间隔三个以上的原子类型都不会改变
    # 因此通过检查group_min_distances是否小于等于3来判断是否需要进行二次反应，添加模拟所需的额外力场参数
    if not is_second_reaction:  # 防止无限递归
        # 检查反应物的group_min_distances是否小于等于3
        for i, p_mol in enumerate(p_mols):
            # 检查group_min_distances是否小于等于3
            if r1_mol.group_min_distances is not None and r1_mol.group_min_distances <= 3:
                # 使用当前产物作为新的反应物，与酰氯再次反应
                p_mol.mol_type = 'r1'
                # 计算分子属性，重新计算反应基团
                p_mol.cal_mol_prop()
                path_dict, second_product_mols = init_product_info(path_dict, p_mol, r2_mol, is_second_reaction=True)
                
                # 将二次反应产物添加到产物列表中
                product_mols.extend(second_product_mols)
                break  # 只对第一个满足条件的产物进行二次反应
    
    return path_dict, product_mols


if __name__ == "__main__":
    smiles1 = "CC1C=CC(N)=CC=1N"
    smiles2 = "O=C(Cl)C1=CC=C(C2=CC(C(=O)Cl)=CC(C(=O)Cl)=C2)C=C1"
    product_smiles_list, byproduct_smiles, atom_indices_dict_list = generate_reaction_smile(smiles1, smiles2)
    print(product_smiles_list, byproduct_smiles, atom_indices_dict_list)

