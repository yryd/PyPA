##############################################################################
# 开发者: Matthew Bone
# 最后更新: 03/08/2021
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
# 将LAMMPS的'read_data'输入文件转换为LAMMPS分子格式文件。
# 应该能读取任何有效格式的LAMMPS输入文件；已移除对头部的假设。
##############################################################################

import os
from AutoMapper.LammpsTreatmentFuncs import clean_data, add_section_keyword, refine_data, save_text_file, format_comment
from AutoMapper.LammpsSearchFuncs import get_data, find_sections, get_header, convert_header

def lammps_to_molecule(directory, fileName, saveName, bondingAtoms: list =None, deleteAtoms=None, validIDSet=None, renumberedAtomDict=None):
    # 切换到文件目录
    os.chdir(directory)

    # 将文件加载到Python中作为列表的列表
    with open(fileName, 'r') as f:
        lines = f.readlines()
    
    # 整理输入
    tidiedLines = clean_data(lines)
    
    # 构建sectionIndexList
    sectionIndexList = find_sections(tidiedLines)

    # 获取原子数据
    atoms = get_data('Atoms', tidiedLines, sectionIndexList)
    atoms = refine_data(atoms, 0, validIDSet, renumberedAtomDict)

    # 获取键数据
    bonds = get_data('Bonds', tidiedLines, sectionIndexList)
    bonds = refine_data(bonds, [2, 3], validIDSet, renumberedAtomDict)
    bondInfo = ('bonds', len(bonds))
    bonds = add_section_keyword('Bonds', bonds)

    # 获取角度数据
    angles = get_data('Angles', tidiedLines, sectionIndexList)
    angles = refine_data(angles, [2, 3, 4], validIDSet, renumberedAtomDict)
    angleInfo = ('angles', len(angles))
    angles = add_section_keyword('Angles', angles)

    # 获取二面角数据
    dihedrals = get_data('Dihedrals', tidiedLines, sectionIndexList)
    dihedrals = refine_data(dihedrals, [2, 3, 4, 5], validIDSet, renumberedAtomDict)
    dihedralInfo = ('dihedrals', len(dihedrals))
    dihedrals = add_section_keyword('Dihedrals', dihedrals)

    # 获取非正常二面角数据
    impropers = get_data('Impropers', tidiedLines, sectionIndexList)
    impropers = refine_data(impropers, [2, 3, 4, 5], validIDSet, renumberedAtomDict)
    improperInfo = ('impropers', len(impropers))
    impropers = add_section_keyword('Impropers', impropers)

    # 重新排列原子数据以获取类型、电荷、坐标 - 假设原子类型非常重要
    types = [[atom[0], atom[2]] for atom in atoms]
    typeInfo = ('atoms', len(types))
    types = add_section_keyword('Types', types)

    charges = [[atom[0], atom[3]] for atom in atoms]
    charges = add_section_keyword('Charges', charges)

    coords = [[atom[0], atom[4], atom[5], atom[6]] for atom in atoms]
    coords = add_section_keyword('Coords', coords)

    # 获取并更改头部值
    header = get_header(tidiedLines)
    
    # 如果提供了新的ID，则用新的数据长度更新数字
    if validIDSet is not None:
        for info in [typeInfo, bondInfo, angleInfo, dihedralInfo, improperInfo]:
            header[info[0]] = [info[1]]

    # 创建键原子注释
    commentString = []
    if bondingAtoms is not None:
        commentString = format_comment(bondingAtoms, 'Bonding_Atoms ')
    if deleteAtoms is not None:
        deleteAtomComment = format_comment(deleteAtoms, 'Delete_Atoms')
        commentString.extend([deleteAtomComment])
    header['comment'].extend(commentString)

    # 移除不必要的头部元素
    keepList = ['comment', 'atoms', 'bonds', 'angles', 'dihedrals', 'impropers']
    cutHeader = {key: header[key] for key in keepList}

    # 将头部转换回字符串列表的列表
    header = convert_header(cutHeader)

    # 合并为一个长的输出列表
    outputList = []
    totalList = [header, types, charges, coords, bonds, angles, dihedrals, impropers]
    
    for keyword in totalList:
        outputList.extend(keyword)
        
    # 输出为文本文件
    save_text_file(saveName, outputList)

