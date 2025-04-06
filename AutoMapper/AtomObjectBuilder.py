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
# 包含关键的原子对象创建和操作工具。Atom 对象类和构建器是映射创建的基础模块。
##############################################################################

import logging
from collections import Counter

from AutoMapper.LammpsSearchFuncs import get_data, find_sections, get_neighbours, get_additional_neighbours
from AutoMapper.LammpsTreatmentFuncs import clean_data

def build_atom_objects(fileName, elementDict, bondingAtoms, createAtoms=[]):
    # 加载分子文件
    with open(fileName, 'r') as f:
        lines = f.readlines()

    # 清理数据并获取坐标和键
    data = clean_data(lines)
    sections = find_sections(data)
    types = get_data('Types', data, sections)

    atomIDs = [row[0] for row in types]
    bonds = get_data('Bonds', data, sections)
    
    # 构建邻居字典
    neighboursDict = get_neighbours(atomIDs, bonds)

    # 移除 createAtoms 作为邻居 - 会混淆映射且不需要
    if createAtoms is not None:
        for keyAtom, neighbours in neighboursDict.items():
            updatedList = []
            for atom in neighbours:
                if atom not in createAtoms:
                    updatedList.append(atom) # 这种方式是因为不能在迭代时更新列表
            
            neighboursDict[keyAtom] = updatedList

    def get_elements(neighbourIDs, elementDict):
        return [elementDict[atomID]for atomID in neighbourIDs]

    atomObjectDict = {}
    for index, atomID in enumerate(atomIDs):
        # 防止 createAtoms 进入对象字典
        if createAtoms is not None:
            if atomID in createAtoms:
                continue

        # 获取原子类型
        atomType = types[index][1]

        # 建立所有邻居关系
        neighbours = neighboursDict[atomID]
        secondNeighbours = get_additional_neighbours(neighboursDict, atomID, neighbours, bondingAtoms)
        thirdNeighbours = get_additional_neighbours(neighboursDict, atomID, secondNeighbours, bondingAtoms)

        neighbourElements = get_elements(neighbours, elementDict)
        secondNeighbourElements = get_elements(secondNeighbours, elementDict)
        thirdNeighbourElements = get_elements(thirdNeighbours, elementDict)

        # 检查原子是否是成键原子，返回布尔值
        if atomID in bondingAtoms:
            bondingAtom = True
        else:
            bondingAtom = False

        atom = Atom(atomID, atomType, elementDict[atomID], bondingAtom, neighbours, secondNeighbours, thirdNeighbours, neighbourElements, secondNeighbourElements, thirdNeighbourElements)
        atomObjectDict[atomID] = atom
    
    return atomObjectDict

def compare_symmetric_atoms(postNeighbourAtomObjectList, preNeighbourAtom, outputType, allowInference=True):
    # 邻居比较 - 无推断
    def compare_neighbours(neighbourLevel):
        neighbourComparison = [getattr(atomObject, neighbourLevel) for atomObject in postNeighbourAtomObjectList]
        neighbourFingerprint = [''.join(sorted(elements)) for elements in neighbourComparison] # 排序以获得字母顺序的指纹

        # 移除重复指纹
        countFingerprints = Counter(neighbourFingerprint)
        tuppledFingerprints = [(index, fingerprint) for index, fingerprint in enumerate(neighbourFingerprint) if countFingerprints[fingerprint] == 1]

        # 如果任何指纹为空(即原子没有X邻居)返回None
        for _, fingerprint in tuppledFingerprints:
            if fingerprint == '':
                return None

        # 任何潜在的后续邻居匹配前原子的指纹，返回后续邻居
        preNeighbourFingerprint = ''.join(sorted(getattr(preNeighbourAtom, neighbourLevel)))
        for index, fingerprint in tuppledFingerprints:
            if preNeighbourFingerprint == fingerprint:
                logging.debug(f'前: {preNeighbourAtom.atomID}, 后: {postNeighbourAtomObjectList[index].atomID} 通过 {neighbourLevel} 找到')
                if outputType == 'index':
                    return index
                elif outputType == 'atomID':
                    return postNeighbourAtomObjectList[index].atomID
                else:
                    print('compare_symmetric_atoms 指定了无效的输出类型')

    # 第一邻居比较
    symmetryResult = compare_neighbours('firstNeighbourElements')

    # 第二邻居比较
    if symmetryResult is None:
        symmetryResult = compare_neighbours('secondNeighbourElements')

    # 第三邻居比较
    if symmetryResult is None:
        symmetryResult = compare_neighbours('thirdNeighbourElements')

    # 如果通过了所有这些检查，猜测分配并警告用户
    if symmetryResult is not None:
        return symmetryResult
    else:
        if allowInference: # 仅在推断开启时
            # 通过将 postNeighbourAtomList 分解为与前原子元素匹配的原子来找到所有潜在选择
            possibleChoices = []
            for index, postNeighbourAtom in enumerate(postNeighbourAtomObjectList):
                if postNeighbourAtom.element == preNeighbourAtom.element:
                    possibleChoices.append((index, postNeighbourAtom.atomID))

            # 让用户知道已进行推断     
            logging.debug(f'前: {preNeighbourAtom.atomID}, 后: {possibleChoices[0][1]} 通过对称推断找到')
            print(
                f'注意: 前键原子ID {preNeighbourAtom.atomID} 已通过推断分配给后键原子ID {possibleChoices[0][1]}。潜在选择是 {[atom[1] for atom in possibleChoices]}。请检查是否正确。'
            )
            if outputType == 'index':
                    return possibleChoices[0][0]
            elif outputType == 'atomID':
                return possibleChoices[0][1]
            else:
                print('compare_symmetric_atoms 指定了无效的输出类型')


class Atom():
    def __init__(self, atomID, atomType, element, bondingAtom, neighbourIDs, secondNeighbourIDs, thirdNeighbourIDs, neighbourElements, secondNeighbourElements, thirdNeighbourElements):
        self.atomID = atomID
        self.atomType = atomType
        self.element = element
        self.bondingAtom = bondingAtom

        # 邻居
        self.mappedNeighbourIDs = neighbourIDs # 根据映射更改
        self.firstNeighbourIDs = neighbourIDs.copy() # 在整个映射过程中固定
        self.secondNeighbourIDs = secondNeighbourIDs
        self.thirdNeighbourIDs = thirdNeighbourIDs

        self.mappedNeighbourElements = neighbourElements # 根据映射更改
        self.firstNeighbourElements = neighbourElements.copy() # 在整个映射过程中固定
        self.secondNeighbourElements = secondNeighbourElements
        self.thirdNeighbourElements = thirdNeighbourElements

    def check_mapped(self, mappedIDs, searchIndex, elementDict):
        """更新邻居ID。

        通过移除已映射的ID来更新neighbourIDs。
        这将在所有邻居映射尝试之前调用，以防止原子被多次映射。

        参数:
            mappedIDs: 此时映射的ID总列表。这将包含前和后原子ID
            searchIndex: 决定使用前还是后原子ID

        返回:
            更新现有的类变量 self.NeighbourIDs
        """
        searchIndexMappedIDs = [row[searchIndex] for row in mappedIDs]
        
        self.mappedNeighbourIDs = [ID for ID in self.mappedNeighbourIDs if ID not in searchIndexMappedIDs]
        self.mappedNeighbourElements = [elementDict[atomID]for atomID in self.mappedNeighbourIDs]


    def map_elements(self, atomObject, preAtomObjectDict, postAtomObjectDict):
        """通过比较相邻元素符号将前原子ID映射到后原子ID。

        比较前原子邻居的字符串化学元素符号的出现与已知后原子邻居的出现。
        基于已知前后原子对的邻居创建新的映射ID列表和任何缺失的ID列表。
        依赖compare_symmetric atoms处理出现次数>1的非H原子。

        参数:
            self: 这是一个已成功映射的前原子
            atomObject: 已映射到前原子的已知后原子
            preAtomObjectDict: 分子中所有前原子的字典
            postAtomObjectDict: 分子中所有后原子的字典

        返回:
            部分映射ID列表，部分缺失的前后原子列表以及队列的额外原子
        """

        # 输出变量
        mapList = []
        missingPreAtoms = []
        queueAtoms = []

        def allowed_maps(preAtom, postAtom):
            '''提前填充缺失原子，而不是基于元素出现做出误导性映射'''
            # 检查元素在前和后原子中出现的次数是否相同
            # 如果不同，不允许进行映射，原子被移动到缺失列表
            preElementOccurences = Counter(preAtom.mappedNeighbourElements)
            postElementOccurences = Counter(postAtom.mappedNeighbourElements)

            allowedMapDict = {}
            for element, count in preElementOccurences.items():
                if count == postElementOccurences[element]:
                    allowedMapDict[element] = True
                else:
                    allowedMapDict[element] = False
                
            # 强制所有H为True，因为氢可以通过推断映射
            if 'H' in allowedMapDict:
                allowedMapDict['H'] = True

            return allowedMapDict

        allowedMapDict = allowed_maps(self, atomObject)

        # 匹配函数
        def matchNeighbour(preAtom, postAtom, preAtomIndex, postAtomIndex, mapList, queueList):
            # 将前后原子ID附加到映射
            mapList.append([preAtom.mappedNeighbourIDs[preAtomIndex], postAtom.mappedNeighbourIDs[postAtomIndex]])
            
            # 将所有非氢原子原子ID添加到队列
            if preAtom.mappedNeighbourElements[preAtomIndex] != 'H':
                queueList.append([preAtom.mappedNeighbourIDs[preAtomIndex], postAtom.mappedNeighbourIDs[postAtomIndex]])

            # 从映射ID和映射元素原子对象值中移除后原子ID
            postAtom.mappedNeighbourIDs.pop(postAtomIndex)
            postAtom.mappedNeighbourElements.pop(postAtomIndex)

        # 循环遍历前原子的邻居并与后原子的邻居比较
        for preIndex, neighbour in enumerate(self.mappedNeighbourElements):
            elementOccurence = atomObject.mappedNeighbourElements.count(neighbour)

            # 检查是否允许与邻居元素的映射，如果不允许则将当前元素添加到缺失列表
            if allowedMapDict[neighbour] == False:
                missingPreAtoms.append(self.mappedNeighbourIDs[preIndex])
                continue

            # 如果在后原子列表中没有匹配，则是缺失的前原子
            if elementOccurence == 0:
                missingPreAtoms.append(self.mappedNeighbourIDs[preIndex])
            
            # 如果只有一个匹配元素，则分配原子ID - 如果元素移动而相同元素占据其位置，这会出错吗？
            elif elementOccurence == 1:
                postIndex = atomObject.mappedNeighbourElements.index(neighbour)
                logging.debug(f'前: {self.mappedNeighbourIDs[preIndex]}, 后: {atomObject.mappedNeighbourIDs[postIndex]} 通过单元素出现找到')
                matchNeighbour(self, atomObject, preIndex, postIndex, mapList, queueAtoms)

            # 多个匹配元素需要额外考虑
            elif elementOccurence > 1:
                if neighbour == 'H': # H可以简单处理，因为在这种情况下所有H彼此等效 - 忽略手性
                    postHydrogenIndexList = [index for index, element in enumerate(atomObject.mappedNeighbourElements) if element == 'H']
                    postIndex = postHydrogenIndexList.pop()
                    logging.debug(f'前: {self.mappedNeighbourIDs[preIndex]}, 后: {atomObject.mappedNeighbourIDs[postIndex]} 通过氢对称推断找到')
                    matchNeighbour(self, atomObject, preIndex, postIndex, mapList, queueAtoms)
                    
                else:
                    # 获取邻居后原子对象
                    postNeighbourIndices = [index for index, val in enumerate(atomObject.mappedNeighbourElements) if val == neighbour]
                    postNeighbourAtomIDs = [atomObject.mappedNeighbourIDs[i] for i in postNeighbourIndices]
                    postNeighbourAtomObjects = [postAtomObjectDict[atomID] for atomID in postNeighbourAtomIDs] 

                    # 获取可能的前原子对象
                    preNeighbourAtomObject = preAtomObjectDict[self.mappedNeighbourIDs[preIndex]]

                    # 为当前前原子找到后原子ID
                    postNeighbourAtomID = compare_symmetric_atoms(postNeighbourAtomObjects, preNeighbourAtomObject, 'atomID')
                    if postNeighbourAtomID is not None:
                        postIndex = atomObject.mappedNeighbourIDs.index(postNeighbourAtomID)
                        matchNeighbour(self, atomObject, preIndex, postIndex, mapList, queueAtoms)
                    else:
                        # 如果找不到后原子，将前原子添加到缺失原子列表
                        print(f'无法为前原子 {self.mappedNeighbourIDs[preIndex]} 找到对称对')
                        missingPreAtoms.append(self.mappedNeighbourIDs[preIndex])                                             


        # 在mapList中搜索缺失的后原子
        mappedPostAtomList = [row[1] for row in mapList]
        missingPostAtoms = [neighbour for neighbour in atomObject.mappedNeighbourIDs if neighbour not in mappedPostAtomList]

        return mapList, missingPreAtoms, missingPostAtoms, queueAtoms
