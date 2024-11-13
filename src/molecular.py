# molecular.py
# 小分子对象与自定义分子指纹生成

from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem import rdmolops, Descriptors
from rdkit.Chem import AllChem
import numpy as np
import logging

# 设置日志配置
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class MolecularModule:
    def __init__(self, smiles : str, g_type : str):
        """
        初始化分子模块，创建分子对象并计算特征属性。
        
        Args:
            smiles: 分子的 SMILES 表示法
            g_type: 分子类型:
                    反应物: r1(氨基), r2(酰氯)
                    溶剂: sol
                    产物: p
                    副产物: byp
        """
        self.smiles = smiles
        self.mol_type = g_type
        self.mol = self.mol_from_smiles(smiles)

        self.group_smiles = self.get_functional_group(self.mol_type)
        self.mol_functional_group = self.functional_group_index(self.mol, self.group_smiles)
        
        self.rings = None
        self.count_group = None
        self.group_min_distances, self.group_max_distances = (None, None)
        self.ar_sp3_balance = None
        self.aromatic_rings = None
        # 生成分子指纹
        self.molfinger = None

    def cal_mol_prop(self):
        self.rings = self.count_rings(self.mol)
        self.count_group = len(self.mol_functional_group)
        self.group_min_distances, self.group_max_distances = self.calculate_distances(self.mol, self.mol_functional_group, self.count_group)
        self.ar_sp3_balance = self.calculate_ar_sp3_balance(self.mol)
        self.aromatic_rings = self.count_aromatic_rings(self.mol)
        # 生成分子指纹
        self.molfinger = self.generate_bit_fingerprint(self.rings, self.count_group, self.group_min_distances, self.group_max_distances, self.aromatic_rings)

    @staticmethod
    def get_functional_group(mol_type):
        """
        根据 mol_type 识别分子官能团。
        分子类型:
        反应物: r1(N), r2(Cl)
        溶剂: sol
        产物: p
        副产物: byp
        
        Returns:
            group_smiles: 分子基团的 SMILES 表示法
        """
        if mol_type == "r1":
            return "N"
        elif mol_type == "r2":
            return "C(=O)Cl"
        elif mol_type == "sol":
            return "C"
        elif mol_type == "p":
            return "C(=O)N"
        elif mol_type == "byp":
            return "Cl"
        else:
            logging.error(f"分子类型错误！输入r1, r2, sol, p, byp以内的字符")
            raise

    @staticmethod
    def functional_group_index(mol, group_smiles):
        """
        查找分子中指定基团的索引。
        
        Returns:
            基团索引列表
        """
        group_mol = Chem.MolFromSmiles(group_smiles)
        if not group_mol:
            logging.error(f"无效的基团 SMILES 表示法")
            raise

        indices = []
        # 处理 CH3 的特殊情况
        if group_smiles == "C":
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 6:  # 碳原子
                    index_tuple = (atom.GetIdx(),)
                    hydrogen_count = 0
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 1:  # 氢原子
                            hydrogen_count += 1
                            index_tuple += (neighbor.GetIdx(),)
                    if hydrogen_count == 3:
                        indices.append(index_tuple)

        # 处理 NH 和 NH2 的情况
        elif group_smiles == "N":
            for atom in mol.GetAtoms():
                if atom.GetAtomicNum() == 7:  # 氮原子
                    index_tuple = (atom.GetIdx(),)
                    hydrogen_count = 0
                    for neighbor in atom.GetNeighbors():
                        if neighbor.GetAtomicNum() == 1:  # 氢原子
                            hydrogen_count += 1
                            index_tuple += (neighbor.GetIdx(),)
                    if hydrogen_count in (1, 2):
                        indices.append(index_tuple)
        else:
            indices = mol.GetSubstructMatches(group_mol)
            indices = list(indices)

        return indices
    
    @staticmethod
    def calculate_distances(rkmol, mol_functional_group, count_group):
        """
        计算分子中指定基团之间的最小和最大距离。

        Returns:
            最小距离和最大距离
        """
        if count_group < 2:
            return -1, -1
        distances = []
        for i in range(count_group):
            for j in range(i + 1, count_group):
                shortest_path = rdmolops.GetShortestPath(rkmol, mol_functional_group[i][0], mol_functional_group[j][0])
                distances.append(len(shortest_path) - 1)
        return min(distances), max(distances)

    @staticmethod
    def count_rings(rkmol):
        """
        计算分子中的环的数量。
        
        Returns:
            环的数量
        """
        return len(Chem.GetSSSR(rkmol))
    
    
    @staticmethod
    def calculate_ar_sp3_balance(rkmol):
        """
        计算芳香(sp2)原子数与 sp2 + sp3 总碳数的比值。
        
        Returns:
            ar / (sp3 + ar) 比值
        """
        num_aromatic_carbon = len(rkmol.GetAromaticAtoms())
        num_sp3_carbon = sum(1 for atom in rkmol.GetAtoms() if str(atom.GetHybridization()) == 'SP3' and atom.GetSymbol() == 'C')
        if (num_sp3_carbon + num_aromatic_carbon) != 0:
            balance = num_aromatic_carbon / (num_sp3_carbon + num_aromatic_carbon)
        else:
            balance = 0
        return balance

    @staticmethod
    def count_aromatic_rings(rkmol):
        """
        计算芳香环的数量。
        
        Returns:
            芳香环的数量
        """
        return Descriptors.NumAromaticRings(rkmol)

    @staticmethod
    def generate_bit_fingerprint(rings, count_group, group_min_distances, group_max_distances, aromatic_rings):
        """
        生成比特位描述的分子指纹，每四个比特位描述一个特征。

        Args:
            rings: 环的数量
            count_group: 基团数量
            group_min_distances: 基团最小距离
            group_max_distances: 基团最大距离
            aromatic_rings: 芳香环的数量

        Returns:
            指纹数组
        """
        # 初始化比特位数组
        fingerprint = [0] * 20  # 5个特征，每个特征4个比特位

        # 设置特征值
        features = {
            'num_rings': rings,                           # 环的数量
            'num_functional_groups': count_group,         # 基团数量
            'min_distance': group_min_distances,          # 基团最小距离
            'max_distance': group_max_distances,          # 基团最大距离
            'num_aromatic_rings': aromatic_rings          # 芳香环的数量
        }

        # 将特征值转换为比特位
        for i, (key, value) in enumerate(features.items()):
            # 确保每个特征值小于16，超出范围则设为15
            if value >= 16:
                value = 15
            
            # 将值转换为4位二进制
            for j in range(4):
                fingerprint[i * 4 + j] = (value >> (3 - j)) & 1  # 取出相应的比特位

        return fingerprint


    @staticmethod
    def rkmol_print(mol, smiles, filepath='./r1'):
        """
        根据 SMLIES 字符串生成分子的 PDB 文件、PNG 图片和 MOL 文件。
        注意 rdkit 产生的对象不包含 H 原子，并且没有初始XYZ坐标
        本程序使用 radonpy 来产生mol文件
        
        Args:
            mol: RDKit分子对象
            smiles: 分子的 SMILES 表示法
            filepath: 分子的保存路径，包含文件名但不包含后缀
        """
        # 不含 H 的图片输出
        simlpe_mol = Chem.MolFromSmiles(smiles)
        Draw.MolToFile(simlpe_mol, filepath + ".png")
        
        # 生成图片保存成文件，添加索引编号
        IPythonConsole.ipython_useSVG = True
        atoms = mol.GetNumAtoms()
        for idx in range(atoms):
            mol.GetAtomWithIdx(idx).SetProp('molAtomMapNumber', str(mol.GetAtomWithIdx(idx).GetIdx()))
        Draw.MolToFile(mol, filepath + "_index.png")
        
        # 写入 PDB 文件
        writer = Chem.PDBWriter(filepath + ".pdb")
        writer.write(mol)
        writer.close()
        
        # 写入 MOL 文件
        mol_block = Chem.MolToMolBlock(mol)
        with open(filepath + ".mol", 'w') as file:
            file.write(mol_block)

    @staticmethod
    def mol_from_smiles(smiles, coord=True, version=2, ez='E', chiral='S'):
        """从SMILES字符串生成RDKit分子对象。该方法源于包：radonpy.core.utils，作者：yhayashi1986/Yoshihiro Hayashi

        Args:
            smiles (str): SMILES表示的分子结构。
            coord (bool): 是否生成3D坐标，默认为True。
            version (int): ETKDG算法的版本，默认为2。
            ez (str): 控制双键的立体化学，'E'表示优先级高的取代基在双键的对面，'Z'表示在同侧。
            chiral (str): 控制手性中心的配置，'S'表示S配置，'R'表示R配置。

        Returns:
            Chem.Mol: 生成的RDKit分子对象，如果转换失败则返回None。
        """
        
        # 统计连接原子的数量
        n_conn = smiles.count('[*]') + smiles.count('*') + smiles.count('[3H]')
        smi = smiles.replace('[*]', '[3H]').replace('*', '[3H]')

        # 选择ETKDG算法的版本
        if version == 3:
            etkdg = AllChem.ETKDGv3()
        elif version == 2:
            etkdg = AllChem.ETKDGv2()
        else:
            etkdg = AllChem.ETKDG()
        
        # 设置ETKDG的参数
        etkdg.enforceChirality = True
        etkdg.useRandomCoords = False
        etkdg.maxAttempts = 100

        # 从SMILES转换为RDKit分子对象
        try:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)  # 添加氢原子
        except Exception as e:
            logging.error(f'无法从 {smiles} 转换为 RDKit 分子对象 : {e}')
            return None

        # 指定立体化学
        Chem.AssignStereochemistry(mol)

        # 获取聚合物主链的原子、键和二面角
        backbone_atoms = []
        backbone_bonds = []
        backbone_dih = []

        if n_conn == 2:
            link_idx = []
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == "H" and atom.GetIsotope() == 3:
                    link_idx.append(atom.GetIdx())
            backbone_atoms = Chem.GetShortestPath(mol, link_idx[0], link_idx[1])

            for i in range(len(backbone_atoms) - 1):
                bond = mol.GetBondBetweenAtoms(backbone_atoms[i], backbone_atoms[i + 1])
                backbone_bonds.append(bond.GetIdx())
                # 检查双键的立体化学
                if bond.GetBondTypeAsDouble() == 2 and str(bond.GetStereo()) == 'STEREONONE' and not bond.IsInRing():
                    backbone_dih.append((backbone_atoms[i - 1], backbone_atoms[i], backbone_atoms[i + 1], backbone_atoms[i + 2]))

        # 列出未指定立体化学的双键（不包括聚合物主链的键和环结构中的键）
        db_list = []
        for bond in mol.GetBonds():
            if bond.GetBondTypeAsDouble() == 2 and str(bond.GetStereo()) == 'STEREONONE' and not bond.IsInRing():
                if n_conn == 2 and bond.GetIdx() in backbone_bonds:
                    continue
                else:
                    db_list.append(bond.GetIdx())

        # 枚举立体异构体
        opts = Chem.EnumerateStereoisomers.StereoEnumerationOptions(unique=True, tryEmbedding=True)
        isomers = tuple(Chem.EnumerateStereoisomers.EnumerateStereoisomers(mol, options=opts))

        if len(isomers) > 1:
            # logging.info('%i 个立体异构体候选已生成', len(isomers))
            chiral_num_max = 0
            
            for isomer in isomers:
                ez_flag = False
                chiral_flag = 0

                Chem.AssignStereochemistry(isomer)

                # 控制未指定双键的立体化学
                ez_list = []
                for idx in db_list:
                    bond = isomer.GetBondWithIdx(idx)
                    if str(bond.GetStereo()) == 'STEREOANY' or str(bond.GetStereo()) == 'STEREONONE':
                        continue
                    elif ez == 'E' and (str(bond.GetStereo()) == 'STEREOE' or str(bond.GetStereo()) == 'STEREOTRANS'):
                        ez_list.append(True)
                    elif ez == 'Z' and (str(bond.GetStereo()) == 'STEREOZ' or str(bond.GetStereo()) == 'STEREOCIS'):
                        ez_list.append(True)
                    else:
                        ez_list.append(False)

                if len(ez_list) > 0:
                    ez_flag = np.all(np.array(ez_list))
                else:
                    ez_flag = True

                # 控制未指定的手性
                chiral_list = np.array(Chem.FindMolChiralCenters(isomer))
                if len(chiral_list) > 0:
                    chiral_centers = chiral_list[:, 0]
                    chirality = chiral_list[:, 1]
                    chiral_num = np.count_nonzero(chirality == chiral)
                    if chiral_num == len(chiral_list):
                        chiral_num_max = chiral_num
                        chiral_flag = 2
                    elif chiral_num > chiral_num_max:
                        chiral_num_max = chiral_num
                        chiral_flag = 1
                else:
                    chiral_flag = 2

                if ez_flag and chiral_flag:
                    mol = isomer
                    # logging.info('更新最优构象')
                    if chiral_flag == 2:
                        break

        # 生成3D坐标
        if coord:
            try:
                enbed_res = AllChem.EmbedMolecule(mol, etkdg)
            except Exception as e:
                logging.error(f'无法生成 %s 的3D坐标：{e}', smiles)
                return None
            if enbed_res == -1:
                etkdg.useRandomCoords = True
                enbed_res = AllChem.EmbedMolecule(mol, etkdg)
                if enbed_res == -1:
                    logging.error(f'无法生成 %s 的3D坐标：{e}', smiles)
                    return None

        # 将聚合物主链中未指定双键的二面角修改为180度
        if len(backbone_dih) > 0:
            for dih_idx in backbone_dih:
                Chem.rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), dih_idx[0], dih_idx[1], dih_idx[2], dih_idx[3], 180.0)

                for na in mol.GetAtomWithIdx(dih_idx[2]).GetNeighbors():
                    na_idx = na.GetIdx()
                    if na_idx != dih_idx[1] and na_idx != dih_idx[3]:
                        break
                Chem.rdMolTransforms.SetDihedralDeg(mol.GetConformer(0), dih_idx[0], dih_idx[1], dih_idx[2], na_idx, 0.0)

        return mol


if __name__ == "__main__":
    # smiles = "O=C(Cl)C1=CC(C(=O)Cl)=CC(C2=CC(C(=O)Cl)=CC(C(=O)Cl)C2)=C1"
    # smiles = "C1CCC(C)C1"
    smiles = "C1CC(C2CCNCC2)CCN1"
    mol = MolecularModule(smiles, 'r1')
    print(mol.mol_functional_group)
    print(mol.group_min_distances, mol.group_max_distances)
    print(mol.count_group)
    print(mol.aromatic_rings)
    print(mol.rings)
    print(mol.ar_sp3_balance)
    print(mol.molfinger)
    
    