# readdata.py
# 该模块用于读取 LAMMPS 数据文件

import re

def clean_data(lines):
    """
    清理读取的行数据，去除空行、注释和多余的空白字符。

    Args:
        lines (list of str): 原始行数据列表。

    Returns:
        list of str: 清理后的行数据列表。
    """
    # 移除空行
    lines = [line for line in lines if line != '\n']

    # 移除注释 - 保留质量标签中的注释，例如 # C_3
    lines = [re.sub(r'(?<!\d\s\s)#(.*)', '' ,line) for line in lines]

    # 移除换行符
    lines = [re.sub(r'\n', '', line) for line in lines]

    # 移除由于注释被移除而导致的空字符串
    lines = [line for line in lines if line != '']

    # 移除尾部空白字符
    lines = [re.sub(r'\s+$', '', line) for line in lines]

    return lines

def get_data(sectionName, lines, sectionIndexList, useExcept=True):
    """
    获取指定节的数据。

    Args:
        sectionName (str): 要获取的节的名称。
        lines (list of str): 清理后的行数据列表。
        sectionIndexList (list of int): 各节的索引列表。
        useExcept (bool): 是否使用异常处理检查节的存在，默认为 True。

    Returns:
        list of list of str: 指定节的数据，以列表形式返回。
    """
    if useExcept:  # 检查节名称是否存在于 LAMMPS 数据中
        try:
            startIndex = lines.index(sectionName)
        except ValueError:
            # 如果不存在，返回空列表，可以正常添加到主列表中
            data = []
            return data
    else:  # 允许后续的 try/except 块捕获缺失的节名称
        startIndex = lines.index(sectionName)

    endIndex = sectionIndexList[sectionIndexList.index(startIndex) + 1]
    
    data = lines[startIndex + 1:endIndex]  # +1 表示节名称不被包含
    data = [val.split() for val in data]

    return data

def find_sections(lines):
    """
    查找节关键词的索引。

    Args:
        lines (list of str): 清理后的行数据列表。

    Returns:
        list of int: 节的索引列表，包括文件末尾的索引。
    """
    # 查找节关键词的索引 - isalpha 确保节关键词中没有空格、换行符或标点符号
    sectionIndexList = [lines.index(line) for line in lines if line.isalpha()]

    # 将文件末尾的索引添加为最后一个索引
    sectionIndexList.append(len(lines))

    return sectionIndexList

def read_data_atomtype(filename):
    """
    从指定文件读取 LAMMPS 数据，并返回原子类型、原子和键的数据。

    Args:
        filename (str): 要读取的 LAMMPS 数据文件的路径。

    Returns:
        list of str: 原子类型数据。
    """
    with open(filename, 'r') as f:
        lines = f.readlines()
    
    # 清理输入数据
    tidiedLines = clean_data(lines)
    # 构建节索引列表
    sectionIndexList = find_sections(tidiedLines)
    # 获取原子类型数据
    types = get_data('Masses', tidiedLines, sectionIndexList)
    return types

def get_map_type_str(filename):
    type_data = read_data_atomtype(filename)
    # 将质量转换为元素符号，以下为 GAFF2 力场所有元素
    elements = {
        12.01: 'C',
        1.008: 'H',
        35.45: 'Cl',
        14.01: 'N',
        16.00: 'O',
        19.00: 'F',
        79.90: 'Br',
        126.9: 'I',
        30.97: 'P',
        32.06: 'S'
    }
    # 创建一个空列表来存储元素
    element_list = []

    for type in type_data:
        mass = float(type[1])
        if mass in elements:
            element_list.append(elements[mass])

    # 将元素列表转换为用空格分隔的字符串
    element_string = ' '.join(element_list)
    return element_string

if __name__ == "__main__":
    str0 = get_map_type_str('tmp/1/data/cleanedpre.data')
    print(str0)