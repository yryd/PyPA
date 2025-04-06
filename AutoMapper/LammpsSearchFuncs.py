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
# 一系列用于搜索LAMMPS文件信息的函数。
# 这些函数适用于'read_data'文件和'molecule'文件
##############################################################################
from natsort import natsorted
from AutoMapper.LammpsTreatmentFuncs import clean_data

# 获取数据
def get_data(sectionName, lines, sectionIndexList, useExcept = True):
    if useExcept: # 检查LAMMPS数据中是否存在该部分名称
        try:
            startIndex = lines.index(sectionName)
        except ValueError:
            # 如果不存在，返回可以稍后正常添加到主列表的空列表
            data = []
            return data

    else: # 允许后续的try/except块捕获缺失的部分名称
        startIndex = lines.index(sectionName)

    endIndex = sectionIndexList[sectionIndexList.index(startIndex) + 1]
    
    data = lines[startIndex+1:endIndex] # +1表示不包括sectionName
    data = [val.split() for val in data]

    return data

def get_coeff(coeffName, settingsData):
    # 输入预分割的数据
    # 返回所有在[0]索引中包含coeffName的行
    coeffs = [line for line in settingsData if line[0] == coeffName]
    
    return coeffs

def find_sections(lines):
    # 查找部分关键字的索引 - isalpha有效，因为部分关键字中没有空格、换行符或标点符号
    sectionIndexList = [lines.index(line) for line in lines if line.isalpha()]

    # 添加文件结尾作为最后一个索引
    sectionIndexList.append(len(lines))

    return sectionIndexList

# 搜索键对
def pair_search(bond, bondAtom):
    '''
    检查键中的任一原子ID是否是所需的原子ID。
    如果未找到匹配项，将返回None。
    '''
    if bond[2] == bondAtom:
        return bond[3]
    elif bond[3] == bondAtom:
        return bond[2]

# 遍历原子ID、可能的键并找到有效键
def search_loop(bonds, bondAtom):
    nextBondAtomList = []

    for searchAtom in bondAtom:
        for bond in bonds:
            nextAtomID = pair_search(bond, searchAtom)
            if nextAtomID is not None:
                nextBondAtomList.append(nextAtomID)
    
    return nextBondAtomList
        
def edge_atom_fingerprint_ids(edgeAtomList, originalBondList, validAtomSet):
    # 获取边缘原子的邻居
    edgeAtomFingerprintDict = get_neighbours(edgeAtomList, originalBondList, []) # 键原子作为空列表给出，边缘原子永远不会有键原子作为邻居，所以不是问题
    
    # 过滤掉部分结构中有效的原子ID
    filteredFingerprintDict = {}
    for key, atomList in edgeAtomFingerprintDict.items():
        cutList = [atom for atom in atomList if atom not in validAtomSet]
        filteredFingerprintDict[key] = cutList

    return filteredFingerprintDict

def get_neighbours(atomIDList, bondsList):
    '''
    获取atomIDList中每个原子的相邻原子ID

    键原子与其他所有原子以相同方式处理。
    '''

    boundAtomsDict = {atom: list() for atom in atomIDList}

    # 遍历键并在字典中构建绑定原子列表
    for bond in bondsList:
        boundAtomsDict[bond[2]].append(bond[3])
        boundAtomsDict[bond[3]].append(bond[2])

    return boundAtomsDict

def get_additional_neighbours(neighboursDict, searchAtomID, searchNeighbours, bondingAtoms, unique=True):
    ''' 获取给定原子ID的邻居的原子ID。     

        此函数设计用于获取给定原子ID的第二和第三邻居。更远的邻居是可能的，
        但可能会产生意外结果。

        参数:
            unique: 防止搜索返回已经存在于neighboursDict中的原子ID，
                如果指定了searchNeighbours，则也包括在内，以及原子ID本身。 

        返回:
            邻居原子ID的列表
    '''
   
    totalNeighbourSet = set()
    for currentNeighbour in searchNeighbours:
        totalNeighbourSet.update(neighboursDict[currentNeighbour])

    if unique:
        # 如果存在，从totalNeighbourSet中移除原始搜索原子ID
        if searchAtomID in totalNeighbourSet:
            totalNeighbourSet.remove(searchAtomID)

        # 移除键原子 - 不想使用键原子指纹，因为它们在前和后总是不同的
        for bondingAtom in bondingAtoms:
            if bondingAtom in totalNeighbourSet:
                totalNeighbourSet.remove(bondingAtom)

        # 从此搜索中移除邻居
        for currentNeighbour in searchNeighbours:
            if currentNeighbour in totalNeighbourSet:
                totalNeighbourSet.remove(currentNeighbour)
        
        # 如果它们不是指定的searchNeighbours，则从集合中移除初始邻居
        # 这适用于>=第三邻居
        if neighboursDict[searchAtomID] != searchNeighbours:
            for neighbour in neighboursDict[searchAtomID]:
                if neighbour in totalNeighbourSet:
                    totalNeighbourSet.remove(neighbour)

    return list(totalNeighbourSet)

def element_atomID_dict(fileName, elementsByType):
    # 加载分子文件
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # 清理数据并获取电荷
    data = clean_data(lines)
    sections = find_sections(data)
    try: # 尝试从分子文件类型获取类型
        types = get_data('Types', data, sections, useExcept=False)
    except ValueError: # 异常从标准lammps文件类型获取类型
        atoms = get_data('Atoms', data, sections, useExcept=False)
        types = [[atomRow[0], atomRow[2]] for atomRow in atoms]
    typesDict = {row[0]: row[1] for row in types} # 键: ID, 值: 类型

    # 确保elementsByType是大写的
    elementsByTypeDict = {index+1: val.upper() for index, val in enumerate(elementsByType)} # 键: 类型, 值: 元素

    # 断言elementsByType中有足够的类型用于types变量中的最高类型
    largestType = int(natsorted(types, key=lambda x: x[1])[-1][1]) # 类型存储为[原子编号, 类型编号]的列表
    assert len(elementsByType) >= largestType, 'EBT (按类型的元素)缺少值。检查所有类型是否存在并用空格分隔。'

    elementIDDict = {key: elementsByTypeDict[int(val)] for key, val in typesDict.items()}

    return elementIDDict

def get_header(tidiedData):
    '''
    从LAMMPS数据文件的头部提取所有数据。
    返回一个关键字键和列出的数值的字典
    '''
    
    # 通过搜索以字母开头的第一行来查找停止行
    def get_stop_line():
        for index, line in enumerate(tidiedData):
            # 检查以跳过初始注释行:
            if index == 0: continue
            if line[0] == '#': continue

            if line[0].isalpha():
                return index

    headerStopLine = get_stop_line()

    # 构建头部部分的字典，关键字为键，数值列表为值
    headerData = tidiedData[0:headerStopLine]
    headerDict = {'comment': []}
    for line in headerData:
        if line[0].isalpha() or line[0] == '#':
            headerDict['comment'].extend([line])
        else:
            # 通过空格分割行
            cutLine = line.split()
            
            # 搜索行以获取数值 - 由于两个盒子尺寸，使用列表
            valueList = []
            keyList = []
            for element in cutLine:
                # 尝试将值转换为int，失败则转换为float，再失败则跳过
                try:
                    valueList.append(int(element))
                except ValueError:
                    try:
                        valueList.append(float(element))
                    except ValueError:
                        keyList.append(element)

            # 从组装的部分创建字典
            headerDict['_'.join(keyList)] = valueList
    
    return headerDict

def convert_header(header):
    '''将头部字典转换回字符串列表的列表以输出'''
    
    stringHeader = []
    for key, values in header.items():
        headerLine = [' '.join([str(val) for val in values])]
        if key != 'comment':
            headerLine.extend(key.split('_'))
        stringHeader.append(headerLine)
    
    return stringHeader