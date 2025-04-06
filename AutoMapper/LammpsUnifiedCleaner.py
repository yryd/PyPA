##############################################################################
# 开发者: Matthew Bone
# 最后更新: 30/07/2021
# 更新者: Matthew Bone
#
# 联系方式:
# Bristol Composites Institute (BCI)
# Department of Aerospace Engineering - University of Bristol
# Queen's Building - University Walk
# Bristol, BS8 1TR
# U.K.
# Email - matthew.bone@bristol.ac.uk
#
# 文件描述:
# 处理任意数量的LAMMPS 'read_data'输入文件和一个系数文件，
# 将它们统一处理以将系数值减少到最低可能值。
# 最初设计用于处理Moltemplate的system.data和system.in.settings文件。
# 如果只指定一个数据文件，此功能类似于cleanup_moltemplate.sh。

# 假设条件:
# LAMMPS原子类型为full
# 单元大小在文件之间相同
# 对系数使用数字定义而不是通配符*
##############################################################################

import os
from natsort import natsorted
from itertools import combinations_with_replacement
from AutoMapper.LammpsTreatmentFuncs import clean_data, clean_settings, add_section_keyword, save_text_file
from AutoMapper.LammpsSearchFuncs import get_data, get_coeff, find_sections, get_header, convert_header

def file_unifier(directory, coeffsFile, dataList):
    # 切换到文件目录
    os.chdir(directory)

    # 从dataList加载文件，整理并初始化类对象
    lammpsData = []
    for dataFile in dataList:
        with open(dataFile, 'r') as f:
            data = f.readlines()

        # 整理数据
        data = clean_data(data)

        headerDict = get_header(data)

        # 初始化数据类
        data = Data(data, headerDict)
        lammpsData.append(data)

    def union_types(typeAttr, lammpsData=lammpsData):
        lammpsTypes = []
        for data in lammpsData:
            # 获取用于确定类型的属性
            getTypeAttr = 'get_' + typeAttr
            func = getattr(data, getTypeAttr)
            # 运行函数并追加到列表
            lammpsTypes.append(func())

        # 合并集合以移除重复项并按数字顺序排序
        types = natsorted(set().union(*lammpsTypes))
        numTypes = (typeAttr, str(len(types))) # 使用元组以便稍后可以在字典中访问类型
        
        return types, numTypes

    # 合并集合并为每种类型创建排序列表
    atomTypes, numAtomTypes = union_types('atom_types')
    bondTypes, numBondTypes = union_types('bond_types')
    angleTypes, numAngleTypes = union_types('angle_types')
    dihedralTypes, numDihedralTypes = union_types('dihedral_types')
    improperTypes, numImproperTypes = union_types('improper_types')

    # 更新各节
    for data in lammpsData:
        bondDict = data.change_section_types(bondTypes, 'bonds')
        angleDict = data.change_section_types(angleTypes, 'angles')
        dihedralDict = data.change_section_types(dihedralTypes, 'dihedrals')
        improperDict = data.change_section_types(improperTypes, 'impropers')
        massDict = data.change_mass_types(atomTypes)
        data.change_atom_types(massDict)

    # 更新头部 - 将删除多行注释并只保留第一个
    sectionTypeCounts = [numAtomTypes, numBondTypes, numAngleTypes, numDihedralTypes, numImproperTypes]
    for data in lammpsData:
        data.change_header(sectionTypeCounts)

    # 保存数据文件
    for index, data in enumerate(lammpsData):
        # 将所有不同的数据源合并为一个列表
        combinedData = [data.header, data.masses, data.atoms, data.bonds, data.angles, data.dihedrals, data.impropers]
        # 将列表的列表扁平化一级
        combinedData = [val for sublist in combinedData for val in sublist]

        # 保存为文本文件
        save_text_file('cleaned' + dataList[index], combinedData)
    
    ####设置部分####

    # 将数据文件加载为列表的列表
    with open(coeffsFile, 'r') as f:
        settings = f.readlines()
    
    # 整理设置并分割
    settings = clean_settings(settings)
    settings = [line.split() for line in settings]

    # 创建原始原子类型pair_coeff对
    originalPairTuples = list(combinations_with_replacement(atomTypes, 2))
    # 获取所有pair_coeffs
    pairCoeff = get_coeff("pair_coeff", settings)
    
    # 找到分子所需的有效的pair_coeff对
    validPairCoeff = []
    for pair in originalPairTuples:
        for coeff in pairCoeff:
            if coeff[1] == pair[0] and coeff[2] == pair[1]:
                validPairCoeff.append(coeff)

    # 使用massDict更新pair_coeffs中的原子类型
    for pair in validPairCoeff:
        pair[1] = massDict[pair[1]]
        pair[2] = massDict[pair[2]]

    def valid_coeffs(coeffType, updateDict, settingsData=settings):
        # 获取系数行
        coeffs = get_coeff(coeffType, settingsData)

        # 从updateDict的键中找到有效的系数
        validCoeffs = []
        for key in updateDict.keys():
            for coeff in coeffs:
                if coeff[1] == key:
                    validCoeffs.append(coeff)
                    break

        # 使用updateDict的值更新系数
        for coeff in validCoeffs:
            coeff[1] = updateDict[coeff[1]]

        return validCoeffs

    # 更新系数值
    validBondCoeff = valid_coeffs('bond_coeff', bondDict)
    validAngleCoeff = valid_coeffs('angle_coeff', angleDict)
    validDihedralCoeff = valid_coeffs('dihedral_coeff', dihedralDict)
    validImproperCoeff = valid_coeffs('improper_coeff', improperDict)

    # 合并所有系数源
    combinedCoeffs = [validPairCoeff, validBondCoeff, validAngleCoeff, validDihedralCoeff, validImproperCoeff]
    # 将列表的列表扁平化一级
    combinedCoeffs = [val for sublist in combinedCoeffs for val in sublist]

    # 保存系数文件
    save_text_file('cleaned' + coeffsFile, combinedCoeffs)

# 处理LAMMPS数据的类
class Data:
    def __init__(self, data, headerDict):
        self.data = data
        self.sectionIndexList = find_sections(self.data)

        # 头部数据
        self.header = headerDict

        # 节数据
        self.atoms = get_data('Atoms', self.data, self.sectionIndexList)
        self.masses = get_data('Masses', self.data, self.sectionIndexList)
        self.bonds = get_data('Bonds', self.data, self.sectionIndexList)
        self.angles = get_data('Angles', self.data, self.sectionIndexList)
        self.dihedrals = get_data('Dihedrals', self.data, self.sectionIndexList)
        self.impropers = get_data('Impropers', self.data, self.sectionIndexList)
    
    def get_atom_types(self):
        atom_types = {atom[2] for atom in self.atoms}
        return atom_types
    
    def get_bond_types(self):
        bond_types = {bond[1] for bond in self.bonds}
        return bond_types

    def get_angle_types(self):
        angle_types = {angle[1] for angle in self.angles}
        return angle_types

    def get_dihedral_types(self):
        dihedral_types = {dihedral[1] for dihedral in self.dihedrals}
        return dihedral_types

    def get_improper_types(self):
        improper_types = {improper[1] for improper in self.impropers}
        return improper_types

    def change_mass_types(self, unioned_atom_types):
        """更新质量类型"""
        # 获取原子中使用的质量
        valid_masses = [mass for mass in self.masses if mass[0] in unioned_atom_types]
        # 获取原子类型
        mass_types = [mass[0] for mass in valid_masses]
        
        # 创建原始原子类型键和新类型值的字典
        new_index_range = list(range(1, len(mass_types)+1))
        new_index_range = [str(val) for val in new_index_range]
        mass_zip = zip(mass_types, new_index_range)
        mass_change_dict = dict(mass_zip)
        
        # 将原子类型更改为新类型
        for massList in valid_masses:
            massList[0] = mass_change_dict[massList[0]]
            # 如果存在注释，在开头添加空格以匹配moltemplate
            if len(massList) >2:
                massList[2] = massList[2].rjust(2)
        
        # 更新self并添加节关键字
        self.masses = add_section_keyword('Masses', valid_masses)

        return mass_change_dict
    
    def change_atom_types(self, mass_change_dict):
        """更新原子类型"""
        for atomList in self.atoms:
            atomList[2] = mass_change_dict[atomList[2]]

        add_section_keyword('Atoms', self.atoms)

    def change_section_types(self, unioned_types, data_section):
        """更新节类型
        
        此函数与change_mass_types不同，因为它不会移除行
        """
        # 使用旧类型键和新类型值构建新字典
        new_index_range = list(range(1, len(unioned_types) + 1))
        new_index_range = [str(val) for val in new_index_range]
        type_zip = zip(unioned_types, new_index_range)
        type_change_dict = dict(type_zip)

        # 更新数据
        sectionData = getattr(self, data_section)
        for entryList in sectionData:
            entryList[1] = type_change_dict[entryList[1]]

        # 添加节关键字
        sentenceCaseSection = data_section.capitalize()
        add_section_keyword(sentenceCaseSection, sectionData)

        return type_change_dict

    def change_header(self, typeList):
        """更新头部信息"""
        # 遍历类型元组并更新头部
        for typeData in typeList:
            self.header[typeData[0]] = [typeData[1]] # 必须是列表，否则>1位数的类型会被空格分隔
        
        # 将列表值转换回字符串
        stringHeader = convert_header(self.header)
        
        # 将头部恢复为字符串列表的列表
        self.header = stringHeader