# generator.py
# 该模块处理反应物 SMILES 字符串生成
from src.molecular import MolecularModule
from rdkit import Chem

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

def generate_reaction_smile(mol1_smiles, mol2_smiles, g_type_smile = 'C(=O)N'):
    """
    生成反应物的 SMILES 字符串。
    
    Args:
        mol1_smiles: 第一个反应物的 SMILES 表示法
        mol2_smiles: 第二个反应物的 SMILES 表示法
    
    Returns:
        tuple: 包含所有种产物 SMILES 列表和副产物 SMILES 的元组
    """
    r1 = MolecularModule(mol1_smiles, 'r1')
    r2 = MolecularModule(mol2_smiles, 'r2')

    # print(r1.smiles)
    # print(r2.smiles)
    # 找到 r1 中所有可反应的胺基索引
    r1_n_indices = r1.mol_functional_group
    
    # 找到 r2 中所有可反应的酰氯基团索引
    r2_c_indices = r2.mol_functional_group


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

            # 4. 生成新的 SMILES
            product_smiles = Chem.MolToSmiles(link_mol, canonical=True)

            product_mol = MolecularModule.mol_from_smiles(product_smiles)
            
            # 重新计算产物基团索引
            p_group_index = list(product_mol.GetSubstructMatches(Chem.MolFromSmiles(g_type_smile)))[0]
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




if __name__ == "__main__":
    smiles1 = "CC1C=CC(N)=CC=1N"
    smiles2 = "O=C(Cl)C1=CC=C(C2=CC(C(=O)Cl)=CC(C(=O)Cl)=C2)C=C1"
    product_smiles_list, byproduct_smiles, atom_indices_dict_list = generate_reaction_smile(smiles1, smiles2)
    print(product_smiles_list, byproduct_smiles, atom_indices_dict_list)

