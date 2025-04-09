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
# AutoMapper 的主要分子和映射生成代码。map_processor 将从常规 LAMMPS 输入文件创建
# 反应前和反应后的分子文件。然后创建完整的映射文件，并决定是否需要创建部分结构。
# 如果需要部分结构，分子文件和映射将被裁剪并重新编号。最后，将生成两个分子文件和
# 一个名为 "automap.data" 的映射文件。
##############################################################################

import os
import logging
import contextlib
from natsort import natsorted
from copy import deepcopy

from AutoMapper.PathSearch import map_from_path
from AutoMapper.LammpsToMolecule import lammps_to_molecule
from AutoMapper.LammpsTreatmentFuncs import save_text_file
from AutoMapper.LammpsSearchFuncs import element_atomID_dict
from AutoMapper.AtomObjectBuilder import build_atom_objects

def map_processor(directory, preDataFileName, postDataFileName, preMoleculeFileName, postMoleculeFileName, preBondingAtoms, postBondingAtoms, deleteAtoms, elementsByType, createAtoms, debug=False, mapFileName='automap.data'):
    # 设置日志级别
    if debug:
        logging.basicConfig(level='DEBUG')
    else:
        logging.basicConfig(level='INFO')
    
    # 分割删除原子列表（如果提供）
    if deleteAtoms is not None:
        assert len(deleteAtoms) % 2 == 0, '错误：为反应前和反应后提供的删除原子ID数量不同。'
        deleteAtomIndex = len(deleteAtoms) // 2
        preDeleteAtoms = deleteAtoms[:deleteAtomIndex]
        postDeleteAtoms = deleteAtoms[deleteAtomIndex:]
    else:
        preDeleteAtoms = None
        postDeleteAtoms = None
    
    # 初始分子创建
    with restore_dir(): # 允许使用相对目录
        lammps_to_molecule(directory, preDataFileName, preMoleculeFileName, preBondingAtoms, deleteAtoms=preDeleteAtoms)
    
    with restore_dir():
        lammps_to_molecule(directory, postDataFileName, postMoleculeFileName, postBondingAtoms, deleteAtoms=postDeleteAtoms)

    # 初始映射创建
    with restore_dir():
        mappedIDList = map_from_path(directory, preMoleculeFileName, postMoleculeFileName, elementsByType, debug, preBondingAtoms, preDeleteAtoms, postBondingAtoms, postDeleteAtoms, createAtoms)

    # 将映射裁剪为最小可能的部分结构
    with restore_dir():
        # 从刚创建的文件加载数据
        os.chdir(directory)
        preElementDict = element_atomID_dict(preMoleculeFileName, elementsByType)
        postElementDict = element_atomID_dict(postMoleculeFileName, elementsByType)

        preAtomObjectDict = build_atom_objects(preMoleculeFileName, preElementDict, preBondingAtoms)
        postAtomObjectDict = build_atom_objects(postMoleculeFileName, postElementDict, postBondingAtoms, createAtoms=createAtoms) 

    # 确定成键原子是否是环的一部分，如果是，确定哪些原子构成环及其邻居
    prePreservedAtomIDs = is_cyclic(preAtomObjectDict, preBondingAtoms, '反应前')
    postPreservedAtomIDs = is_cyclic(postAtomObjectDict, postBondingAtoms, '反应后')
    
    # 确定如果反应是开环反应需要保留的原子
    prePartialAtomsSet, postPartialAtomsSet = is_ring_opening(prePreservedAtomIDs, postPreservedAtomIDs, mappedIDList)

    # 保留距离成键原子最多4个键的原子
    prePartialAtomsSet = keep_all_neighbours(preAtomObjectDict, preBondingAtoms, prePartialAtomsSet)
    postPartialAtomsSet = keep_all_neighbours(postAtomObjectDict, postBondingAtoms, postPartialAtomsSet)

    # 保留删除原子
    if preDeleteAtoms is not None:
        prePartialAtomsSet.update(preDeleteAtoms)
        postPartialAtomsSet.update(postDeleteAtoms)

    # 在后原子集中保留创建原子
    if createAtoms is not None:
        postPartialAtomsSet.update(createAtoms)

    # 查找初始反应前边缘原子
    preEdgeAtoms = find_edge_atoms(preAtomObjectDict, prePartialAtomsSet)

    # 检查边缘是否太靠近类型变化的原子
    preExtendEdgeDict = verify_edge_atoms(preEdgeAtoms, mappedIDList, preAtomObjectDict, postAtomObjectDict)

    # 如果边缘原子需要扩展，则更新映射列表和部分原子集
    mappedIDList, prePartialAtomsSet, postPartialAtomsSet = extend_edge_atoms(preExtendEdgeDict, mappedIDList, preAtomObjectDict, postAtomObjectDict, prePartialAtomsSet, postPartialAtomsSet)
    
    # 在可能的扩展后重新查找边缘原子
    preEdgeAtoms = find_edge_atoms(preAtomObjectDict, prePartialAtomsSet)

    # 检查并获取不是deleteIDs的副产物原子
    postAtomByproducts = get_byproducts(postAtomObjectDict, postBondingAtoms)
    if postAtomByproducts is not None:
        logging.debug(f'找到副产物。副产物是 {postAtomByproducts}（后ID）')
        postPartialAtomsSet.update(postAtomByproducts)

    # 按preAtomID排序mappedIDList
    mappedIDList = natsorted(mappedIDList, key=lambda x: x[0])

    # 创建空的partialMappedIDList以填充返回值
    partialMappedIDList = []

    # 如果部分结构与完整结构长度不同，则重新编号映射
    # 如果它们相等，则只输出映射，不需要更改
    if len(prePartialAtomsSet) != len(preAtomObjectDict):
        logging.debug(f'创建部分映射。')
        # 构建部分映射并获取重新编号字典
        mappedIDList, preRenumberdAtomDict, postRenumberedAtomDict, partialMappedIDList = create_partial_map(mappedIDList, prePartialAtomsSet, postPartialAtomsSet)

        # 重新编号分子创建和输出的关键特征
        preBondingAtoms = renumber(preBondingAtoms, preRenumberdAtomDict)
        preDeleteAtoms = renumber(preDeleteAtoms, preRenumberdAtomDict)

        postBondingAtoms = renumber(postBondingAtoms, postRenumberedAtomDict)
        postDeleteAtoms = renumber(postDeleteAtoms, postRenumberedAtomDict)

        # 如果有边缘原子，则重新编号
        if preEdgeAtoms is not None:
            preEdgeAtoms = renumber(preEdgeAtoms, preRenumberdAtomDict)

        # 如果包含创建原子，则确定重新编号的原子字典并重新编号
        if createAtoms is not None:
            IDCounter = max([int(val) for val in postRenumberedAtomDict.values()]) # 将计数器初始化为重新编号原子的最高ID

            # 为createAtoms扩展重新编号的原子字典
            for createAtom in createAtoms:
                IDCounter += 1
                postRenumberedAtomDict[createAtom] = str(IDCounter)

            # 使用新的重新编号原子字典重新编号createAtoms
            createAtoms = renumber(createAtoms, postRenumberedAtomDict)


        # 使用部分结构重建分子文件
        with restore_dir():
            lammps_to_molecule(directory, preDataFileName, preMoleculeFileName, preBondingAtoms, deleteAtoms=preDeleteAtoms, validIDSet=prePartialAtomsSet, renumberedAtomDict=preRenumberdAtomDict)

        with restore_dir():
            lammps_to_molecule(directory, postDataFileName, postMoleculeFileName, postBondingAtoms, deleteAtoms=postDeleteAtoms, validIDSet=postPartialAtomsSet, renumberedAtomDict=postRenumberedAtomDict)

    # 输出映射文件
    with restore_dir():
        os.chdir(directory)
        outputData = output_map(mappedIDList, preBondingAtoms, preEdgeAtoms, preDeleteAtoms, createAtoms)
        save_text_file(mapFileName, outputData)

    # 返回mappedIDList供其他函数使用，例如测试
    return [mappedIDList, partialMappedIDList]

def output_map(mappedIDList, preBondingAtoms, preEdgeAtoms, preDeleteAtoms, createAtoms):
    # 成键原子
    bondingAtoms = [['\n', 'BondingIDs', '\n']]
    for atom in preBondingAtoms:
        bondingAtoms.extend([[atom]])
    bondingAtoms.extend(['\n'])

    # 删除原子
    deleteIDCount = []
    deleteAtoms = []
    if preDeleteAtoms is not None:
        deleteIDCount.extend([[str(len(preDeleteAtoms)) + ' deleteIDs']])
        deleteAtoms.extend([['DeleteIDs', '\n']])
        for atom in preDeleteAtoms:
            deleteAtoms.extend([[atom]])
        deleteAtoms.extend(['\n'])

    # 边缘原子
    edgeIDCount = []
    edgeAtoms = []
    if preEdgeAtoms is not None:
        edgeIDCount.extend([[str(len(preEdgeAtoms)) + ' edgeIDs']])
        edgeAtoms.extend([['EdgeIDs', '\n']])
        for atom in preEdgeAtoms:
            edgeAtoms.extend([[atom]])
        edgeAtoms.extend(['\n'])

    # 创建原子
    createIDCount = []
    outputCreateAtoms = []
    if createAtoms is not None:
        createIDCount.extend([[str(len(createAtoms)) + ' createIDs']])
        outputCreateAtoms.extend([['CreateIDs', '\n']])
        for atom in createAtoms:
            outputCreateAtoms.extend([[atom]])
        outputCreateAtoms.extend(['\n'])

    # 等价关系
    equivalences = [['#这是由AutoMapper生成的映射\n'], [str(len(mappedIDList)) + ' equivalences']]
    equivalenceAtoms = [['Equivalences', '\n']]
    for atomPair in mappedIDList:
        equivalenceAtoms.extend([[atomPair[0] + '\t' + atomPair[1]]])

    # 输出数据
    output = []
    totalOutput = [equivalences, deleteIDCount, edgeIDCount, createIDCount, bondingAtoms, deleteAtoms, edgeAtoms, outputCreateAtoms, equivalenceAtoms]

    for section in totalOutput:
        output.extend(section)

    return output

# 用于移动到不同的操作系统路径然后返回原始目录的实用程序
@contextlib.contextmanager
def restore_dir():
    startDir = os.getcwd()
    try:
        yield
    finally:
        os.chdir(startDir)

def bfs(graph, startAtom, endAtom, breakLink=False):
    # 改编自 https://stackoverflow.com/questions/8922060/how-to-trace-the-path-in-a-breadth-first-search
    # 用于跟踪atomID是否已被看到的列表
    discovered = {key: False for key in graph.keys()}

    discovered[startAtom] = True

    newGraph = deepcopy(graph)
     # 如果存在，则断开起始原子和目标原子之间的链接 - 在搜索循环时防止搜索向后进行
    if breakLink:
        newGraph[startAtom].remove(endAtom)

    # 迭代路径同时保持所有路径的记录
    queue = []

    queue.append([startAtom])

    while queue:
        path = queue.pop(0)
        
        # 获取最新的路径元素
        node = path[-1]

        if node == endAtom:
            return path

        for neighbour in newGraph.get(node, []):
            # 防止路径陷入循环
            if discovered[neighbour]:
                continue
            
            # 通过下一个邻居增加路径并添加到队列
            discovered[neighbour] = True
            newPath = list(path)
            newPath.append(neighbour)
            queue.append(newPath)
    
    # 如果到达这里，则未找到路径
    return None

def is_cyclic(atomObjectDict, bondingAtoms, reactionType):
    # 创建相邻键的字典 - 改进：从相邻键中移除H，因为它们不能去任何地方
    moleculeGraph = {atom.atomID: atom.firstNeighbourIDs for atom in atomObjectDict.values()}
    preservedAtomIDs = {} # 相对于每个成键原子

    for bondingAtom in bondingAtoms:
        # 获取起始邻居
        startNeighbours = atomObjectDict[bondingAtom].firstNeighbourIDs

        # 设置preservedAtomIDs
        preservedAtomIDs[bondingAtom] = None
        
        # 迭代邻居直到找到一个循环
        for startAtom in startNeighbours:
            cyclicPath = bfs(moleculeGraph, startAtom, bondingAtom, breakLink=True)

            if cyclicPath is not None:
                logging.debug(f'找到循环：{cyclicPath}。从{startAtom}开始，用于{reactionType}反应。')
                break

        # 找到的路径将被转换为需要保留的atomID集合
        if cyclicPath is not None:
            preservedIDsSet = set()
            for atomID in cyclicPath:
                preservedIDsSet.add(atomID)
                preservedIDsSet.update(atomObjectDict[atomID].firstNeighbourIDs)
            
            preservedAtomIDs[bondingAtom] = preservedIDsSet

    return preservedAtomIDs

def find_mapped_pair(preAtom, mappedIDList):
    # 从映射中获取给定preAtom的post原子
    for pair in mappedIDList:
        if pair[0] == preAtom:
            return pair[1]

def is_ring_opening(prePreservedAtomIDs, postPreservedAtomIDs, mappedIDList):
    '''
    确定反应是否是开环聚合反应。
    返回两个字典，包含成键原子键和保留的原子集。
    '''

    # 初始化用于存储环的字典
    preCyclicAtomsSet = set()
    postCyclicAtomsSet = set()

    for preBondingAtom, prePreservedIDSet in prePreservedAtomIDs.items():
        
        # 如果preBondingAtom是环状的（不是None），获取post成键原子
        if prePreservedIDSet is not None:
            postBondingAtom = find_mapped_pair(preBondingAtom, mappedIDList)

            # 如果post成键原子不是环状的，则假定为开环聚合反应
            if postPreservedAtomIDs[postBondingAtom] is None:
                logging.debug(f'反应被确定为开环反应。')

                # 存储pre键数据
                preCyclicAtomsSet.add(preBondingAtom)
                preCyclicAtomsSet.update(prePreservedIDSet)

                # 从映射获取post键数据并存储
                # 这样做是因为开环反应的postPreservedAtomsIDs将为None
                postCyclicAtomsSet.add(postBondingAtom)
                for preCyclicAtom in prePreservedIDSet:
                    postCyclicAtom = find_mapped_pair(preCyclicAtom, mappedIDList)
                    postCyclicAtomsSet.add(postCyclicAtom)

    return preCyclicAtomsSet, postCyclicAtomsSet

def keep_all_neighbours(atomObjectDict, bondingAtoms, partialAtomSet):
    for bondingAtom in bondingAtoms:
        # 添加成键原子到集合中，以防它尚未存在
        partialAtomSet.add(bondingAtom)

        # 获取原子对象
        atomObject = atomObjectDict[bondingAtom]

        # 将neighbourIDs添加到部分集合中
        partialAtomSet.update(atomObject.firstNeighbourIDs)
        partialAtomSet.update(atomObject.secondNeighbourIDs)
        partialAtomSet.update(atomObject.thirdNeighbourIDs)

    return partialAtomSet

def create_partial_map(mappedIDList, prePartialAtomsSet, postPartialAtomsSet):
    # 移除不在pre和post部分原子集中的ID
    partialMappedIDList = []
    for pair in mappedIDList:
        if pair[0] in prePartialAtomsSet and pair[1] in postPartialAtomsSet:
            partialMappedIDList.append(pair)
        
        # 如果某物在一个集合中存在，它必须在另一个集合中存在
        if pair[0] in prePartialAtomsSet and pair[1] not in postPartialAtomsSet:
            print(f'警告：Pre原子{pair[0]}存在但Post原子{pair[1]}缺失')
        if pair[0] not in prePartialAtomsSet and pair[1] in postPartialAtomsSet:
            print(f'警告：Pre原子{pair[0]}缺失但Post原子{pair[1]}存在')

        # 调试工具
        # if pair[0] not in prePartialAtomsSet:
        #     print(f'Pre原子{pair[0]}不在部分原子中')
        # if pair[1] not in postPartialAtomsSet:
        #     print(f'Post原子{pair[1]}不在部分原子中')

    preRenumberedAtomDict = {}
    postRenumberedAtomDict = {}
    renumberedMappedIDList = []
    # 简单地根据它们在ID列表中的位置重新编号pre和post原子
    for index, pair in enumerate(partialMappedIDList, start=1):
        preRenumberedAtomDict[pair[0]] = str(index)
        postRenumberedAtomDict[pair[1]] = str(index)
        renumberedMappedIDList.append([str(index), str(index)])

    # 断言pre和post字典中的原子数量相同。否则Lammps将失败
    assert len(preRenumberedAtomDict) == len(postRenumberedAtomDict), '在反应前和反应后的部分结构中发现了不同数量的原子。\n请检查您的输入文件，如果问题仍然存在，请在Github上提出问题。'

    return renumberedMappedIDList, preRenumberedAtomDict, postRenumberedAtomDict, partialMappedIDList

def renumber(inputList, renumberedAtomDict):
    '''
    接受数字列表和转换字典作为输入。
    输出相同顺序的转换后的数字列表。如果输入为None，则返回None
    '''
    if inputList is None:
        return None

    outputList = []
    for value in inputList:
        outputList.append(renumberedAtomDict[value])

    return outputList

def find_edge_atoms(atomObjectDict, partialAtomSet):
    edgeAtoms = []
    for atom in partialAtomSet:
        # 跳过H原子，因为它们不能是边缘原子
        if atomObjectDict[atom].element == 'H':
            continue
        
        # 迭代第一个邻居并检查它们是否都在部分原子集中
        for neighbour in atomObjectDict[atom].firstNeighbourIDs:
            # 如果一个邻居不在部分原子集中，则原子必须是边缘
            if neighbour not in partialAtomSet:
                edgeAtoms.append(atom)
                break

    if len(edgeAtoms) > 0:
        return edgeAtoms
    else: # 如果分子中没有边缘原子，则其他函数期望None
        return None

def verify_edge_atoms(preEdgeAtoms, mappedIDList, preAtomObjectDict, postAtomObjectDict):
    # 如果没有给出边缘原子，则返回一个空的扩展列表
    if preEdgeAtoms is None:
        return {}

    # 将mappedIDList转换为mappedIDDict
    mappedIDDict = {pair[0]: pair[1] for pair in mappedIDList}

    # 比较pre和post原子类型是否相同，如果不同则返回True
    def compare_atom_type(preAtom):
        preAtomType = preAtomObjectDict[preAtom].atomType
        pairAtom = mappedIDDict[preAtom]
        postAtomType = postAtomObjectDict[pairAtom].atomType

        if preAtomType != postAtomType:
            return True
        else:
            return False

    # 检查原子类型变化是否太靠近边缘原子
    extendDistanceDict = {}
    for edgeAtom in preEdgeAtoms:
        # 边缘原子
        stopSearch = compare_atom_type(edgeAtom)

        if stopSearch:
            extendDistanceDict[edgeAtom] = 3
            continue

        # 第一个邻居
        preEdgeAtomObject = preAtomObjectDict[edgeAtom]
        for firstNeighbour in preEdgeAtomObject.firstNeighbourIDs:
            
            stopSearch = compare_atom_type(firstNeighbour)
            
            if stopSearch:
                extendDistanceDict[edgeAtom] = 2
                break

        # 第二个邻居
        if stopSearch: continue # 防止第二个邻居运行并覆盖第一个邻居的结果
        for secondNeighbour in preEdgeAtomObject.secondNeighbourIDs:
            stopSearch = compare_atom_type(secondNeighbour)
            
            if stopSearch:
                extendDistanceDict[edgeAtom] = 1
                break

    return extendDistanceDict

def extend_edge_atoms(extendEdgeDict, mappedIDList, preAtomObjectDict, postAtomObjectDict, prePartialAtomsSet, postPartialAtomsSet):
    # 输出扩展的mappedIDList，pre和post部分原子集。重新运行find_edge_atoms以获取新边缘。
    # 如果一个边缘在这个列表中，它至少需要扩展一个
    additionalPreAtoms = []
    additionalPostAtoms = []

    for preEdge, extendDist in extendEdgeDict.items():
        # 对于pre-bond
        additionalPreAtoms.extend(preAtomObjectDict[preEdge].firstNeighbourIDs)

        # 对于post-bond
        postEdge = find_mapped_pair(preEdge, mappedIDList)
        additionalPostAtoms.extend(postAtomObjectDict[postEdge].firstNeighbourIDs)

        # 当进一步邻居被要求
        if extendDist == 2:
            additionalPreAtoms.extend(preAtomObjectDict[preEdge].secondNeighbourIDs)
            additionalPostAtoms.extend(postAtomObjectDict[postEdge].secondNeighbourIDs)
        
        if extendDist == 3:
            additionalPreAtoms.extend(preAtomObjectDict[preEdge].thirdNeighbourIDs)
            additionalPostAtoms.extend(postAtomObjectDict[postEdge].thirdNeighbourIDs)            
        
    # 扩展mappedIDList
    for preAtom in additionalPreAtoms:
        # 可以使用additionalPostAtoms列表完成，但这更容易
        postAtom = find_mapped_pair(preAtom, mappedIDList)
        if [preAtom, postAtom] not in mappedIDList: # 防止添加重复的映射对
            mappedIDList.append([preAtom, postAtom])

    # 更新部分原子集
    prePartialAtomsSet.update(additionalPreAtoms)
    postPartialAtomsSet.update(additionalPostAtoms)

    return mappedIDList, prePartialAtomsSet, postPartialAtomsSet

def get_byproducts(postAtomObjectDict, postBondingAtoms):
    # 确定post结构中是否有从每个post原子到成键原子的路径
    # 如果没有路径存在，则原子必须来自不被删除的副产物
    byproducts = []
    targetBondingAtom = postBondingAtoms[0] # 只需要一个原子，因为它将与另一个原子绑定
    moleculeGraph = {atom.atomID: atom.firstNeighbourIDs for atom in postAtomObjectDict.values()}
    for startAtom in postAtomObjectDict.keys():
        pathToBondingAtom = bfs(moleculeGraph, startAtom, targetBondingAtom, breakLink=False)

        if pathToBondingAtom is None:
            byproducts.append(startAtom)

    if len(byproducts) > 0:
        return byproducts
    else:
        return None
