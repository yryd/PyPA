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
# 一系列用于预处理、搜索和保存 LAMMPS 'read_data' 输入文件的函数。
# 为 AutoMapper 处理而构建，但可应用于广泛的问题
##############################################################################

import re # 用于 clean_data, clean_settings
from collections import Counter # 用于 refine_data
from operator import itemgetter # 用于 refine_data
from natsort import natsorted # 用于 refine_data

# 函数可能稍后移至通用函数文件
def clean_data(lines):
    # 移除空行
    lines = [line for line in lines if line != '\n']

    # 移除注释 - 负向后视断言意味着质量中的标签注释被保留，例如 # C_3
    lines = [re.sub(r'(?<!\d\s\s)#(.*)', '' ,line) for line in lines]

    # 移除换行符
    lines = [re.sub(r'\n', '', line) for line in lines]

    # 移除由注释移除导致的列表中的空字符串
    lines = [line for line in lines if line != '']

    # 移除尾随空格
    lines = [re.sub(r'\s+$', '', line) for line in lines]

    return lines

def clean_settings(lines):
    # 移除换行符
    lines = [re.sub(r'\n', '', line) for line in lines]

    # 移除制表符
    lines = [re.sub(r'\t', '', line) for line in lines]

    # 将多个空格替换为一个
    lines = [re.sub(r'\s{2,}', ' ', line) for line in lines]

    return lines

def refine_data(data, searchIndex: list, IDset=None, newAtomIDs=None):
    '''
    搜索多个索引以匹配 atomID 值。
    如果找到匹配项，则在数据中保留该列表行。
    仅当行出现 len(searchIndex) 次时才输出该行
    这意味着数据行在所有可能的 ID 位置中包含有效的 atomID
    '''
    # 如果未给出 IDSet，则跳过精炼
    if IDset is None:
        return data

    # 方便将整数值转换为列表
    if type(searchIndex) is not list:
        searchIndex = [searchIndex]

    # 将包含有效原子 ID 的行添加到列表
    possibleValidData = []
    for atomID in IDset:
        for row in data:
            for index in searchIndex:
                if row[index] == atomID:
                    possibleValidData.append(row)

    # 在上述搜索中找到的 Lammps ID
    # 例如，对于角度，如果角度 '25' 出现 3 次，则意味着角度 '25' 的组件都在 IDSet 中
    possibleValidIDs = [val[0] for val in possibleValidData]
    IDCount = dict(Counter(possibleValidIDs))
    # 如果 ID 计数器与搜索索引大小相同，则 ID 有效并添加到数据中
    validIDs = [key for key in IDCount.keys() if IDCount[key] == len(searchIndex)]
    validData = []
    for row in possibleValidData:
        if row[0] in validIDs:
            validData.append(row)
            # 防止将来重新发现重复的 ID
            validIDIndex = validIDs.index(row[0])
            del validIDs[validIDIndex]

    # 如果提供了字典，则使用它将原子 ID 从旧的更新为新的
    for rowInd, row in enumerate(validData, start=1):
        if searchIndex != [0]: # 不要为 'atoms' 部分运行此操作
            row[0] = str(rowInd) # 重置 LAMMPS ID，例如键号
        for index in searchIndex:
            row[index] = newAtomIDs[row[index]]

    # 按 ID 重新排序 validData，使用 natsort 因为值是 str 而不是 int
    validData = natsorted(validData, key=itemgetter(0))

    return validData

def add_section_keyword(sectionName, data):
    # 如果列表为空，则不添加关键字 - 空列表意味着文件中没有该部分
    if len(data) == 0:
        return data

    # 在列表开头添加关键字名称
    data.insert(0, '\n')
    data.insert(1, [sectionName])
    data.insert(2, '\n')

    return data

def save_text_file(fileName, dataSource):
    # 保存到文本文件
    with open(fileName, 'w') as f:
        for item in dataSource:
            line = " ".join(item)
            if line != '\n':
                f.write("%s\n" % line)
            else:
                f.write(line)

# 创建带有键原子和边缘原子的注释字符串
def format_comment(IDlist, comment):
    atomList = list(IDlist)
    atomList.insert(0, comment)
    atomString = [' '.join(atomList)] # 必须是列表的列表才能通过后续代码

    return atomString

# 将边缘原子指纹从原子 ID 转换为元素字符串
def edge_atom_fingerprint_strings(edgeAtomFingerprintDict, elementsByTypeDict):
    edgeElementFingerprintDict = {}
    for key, atomList in edgeAtomFingerprintDict.items():
        cutList = [elementsByTypeDict[atom] for atom in atomList]

        # 按字母顺序排序列表
        cutList = natsorted(cutList)

        # 将列表连接为单个字符串
        cutList = ''.join(cutList)

        edgeElementFingerprintDict[key] = cutList

    return edgeElementFingerprintDict
