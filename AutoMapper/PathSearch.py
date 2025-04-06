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
# 这是主要的映射代码，使用自定义路径搜索来确定键合前后分子文件中相似的原子。
# 该代码只能从MapProcessor调用，不能独立创建映射。
##############################################################################

import os
import logging
import sys

from AutoMapper.LammpsSearchFuncs import element_atomID_dict
from AutoMapper.AtomObjectBuilder import build_atom_objects, compare_symmetric_atoms
from AutoMapper.QueueFuncs import Queue, queue_bond_atoms, run_queue

def map_delete_atoms(preDeleteAtoms, postDeleteAtoms, mappedIDList):
    """将删除的原子添加到映射列表中"""
    if preDeleteAtoms is not None:
        assert postDeleteAtoms is not None, '在前键文件中找到删除原子但在后键文件中未找到'
        assert len(preDeleteAtoms) == len(postDeleteAtoms), '前后键文件中的删除原子数量不同'
        for index, preAtom in enumerate(preDeleteAtoms):
            mappedIDList.append([preAtom, postDeleteAtoms[index]])
            logging.debug(f'前: {preAtom}, 后: {postDeleteAtoms[index]} 通过用户指定的删除原子找到')

def get_missing_atom_objects(missingAtomList, atomObjectDict):
    """获取缺失原子的对象列表"""
    missingAtomObjects = []
    for atom in missingAtomList:
        atomObject = atomObjectDict[atom]
        missingAtomObjects.append(atomObject)

    return missingAtomObjects

def map_missing_atoms(missingPreAtomObjects, missingPostAtomObjects, mappedIDList, queue, allowInference):
    """映射缺失的原子"""
    missingCheckCounter = 1
    while missingCheckCounter < 4 and len(missingPostAtomObjects) > 0:
        mappedPreAtomIndex = []
        for preIndex, preAtom in enumerate(missingPreAtomObjects):
            missingPostAtomElements = [atom.element for atom in missingPostAtomObjects]
            elementOccurence = missingPostAtomElements.count(preAtom.element)

            if elementOccurence == 0:
                print(f"错误: 无法在后缺失原子中找到 {preAtom.atomID} 的匹配项。请重试或手动映射该原子")

            elif elementOccurence == 1:
                postIndex = missingPostAtomElements.index(preAtom.element)
                logging.debug(f'前: {preAtom.atomID}, 后: {missingPostAtomObjects[postIndex].atomID} 通过单一元素出现找到')
                match_missing(preAtom, postIndex, missingPostAtomObjects, mappedIDList, queue, preIndex, mappedPreAtomIndex)
                
            elif elementOccurence > 1:
                if preAtom.element == 'H':
                    postHydrogenIndexList = [index for index, element in enumerate(missingPostAtomElements) if element == 'H']
                    postIndex = postHydrogenIndexList.pop()
                    logging.debug(f'前: {preAtom.atomID}, 后: {missingPostAtomObjects[postIndex].atomID} 通过氢对称推断找到')
                    match_missing(preAtom, postIndex, missingPostAtomObjects, mappedIDList, queue, preIndex, mappedPreAtomIndex)
                else:
                    potentialPostAtomObjects = [atomObject for atomObject in missingPostAtomObjects if atomObject.element == preAtom.element]
                    postIndex = compare_symmetric_atoms(potentialPostAtomObjects, preAtom, 'index', allowInference=allowInference)
                    if postIndex is not None:
                        match_missing(preAtom, postIndex, potentialPostAtomObjects, mappedIDList, queue, preIndex, mappedPreAtomIndex)
                        logging.debug('上述原子ID对通过缺失原子对称比较找到')

        # 刷新missingPreAtomObjects以避免在后续循环中打印不必要的错误消息
        for index in sorted(mappedPreAtomIndex, reverse=True):
            del missingPreAtomObjects[index]        
        
        missingCheckCounter += 1

def match_missing(preAtom, postAtomMissingIndex, missingPostAtomObjects, mapList, queue, preIndex, mappedPreAtomIndex):
    """匹配缺失的原子并更新映射列表和队列"""
    postAtom = missingPostAtomObjects[postAtomMissingIndex]

    mapList.append([preAtom.atomID, postAtom.atomID])
    
    if preAtom.element != 'H':
        queue.add([[preAtom, postAtom]]) # 绕过add_to_queue()

    missingPostAtomObjects.pop(postAtomMissingIndex)
    mappedPreAtomIndex.append(preIndex)

def update_missing_list(missingAtomList, mappedIDList, mapIndex):
    """更新缺失原子列表"""
    mappedAtoms = [pair[mapIndex] for pair in mappedIDList]
    newMissingAtomList = [atom for atom in missingAtomList if atom not in mappedAtoms]

    return newMissingAtomList

def map_from_path(directory, preFileName, postFileName, elementsByType, debug, preBondingAtoms, preDeleteAtoms, postBondingAtoms, postDeleteAtoms, createAtoms):
    """主路径搜索映射函数"""
    if debug:
        logging.basicConfig(level='DEBUG')
    else:
        logging.basicConfig(level='INFO')

    os.chdir(directory)
    preElementDict = element_atomID_dict(preFileName, elementsByType)
    postElementDict = element_atomID_dict(postFileName, elementsByType)
    elementDictList = [preElementDict, postElementDict]

    preAtomObjectDict = build_atom_objects(preFileName, preElementDict, preBondingAtoms)
    postAtomObjectDict = build_atom_objects(postFileName, postElementDict, postBondingAtoms, createAtoms=createAtoms)

    if createAtoms is None:
        assert len(preAtomObjectDict) == len(postAtomObjectDict), f'前后键文件中的原子数量不同。前: {len(preAtomObjectDict)}, 后: {len(postAtomObjectDict)}'

    missingPreAtomList = []
    missingPostAtomList = []
    mappedIDList = []

    queue = Queue()

    queue_bond_atoms(preAtomObjectDict, preBondingAtoms, postAtomObjectDict, postBondingAtoms, mappedIDList, queue)

    map_delete_atoms(preDeleteAtoms, postDeleteAtoms, mappedIDList)

    run_queue(queue, mappedIDList, preAtomObjectDict, postAtomObjectDict, missingPreAtomList, missingPostAtomList, elementDictList)

    missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)

    timeoutCounter = 1
    inference = False
    while len(missingPreAtomList) > 0 and timeoutCounter < 11:
        missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)
        missingPostAtomList = update_missing_list(missingPostAtomList, mappedIDList, 1)

        missingPreAtomObjects = get_missing_atom_objects(missingPreAtomList, preAtomObjectDict)
        missingPreAtomCount = len(missingPostAtomList)

        mappedPostAtoms = [pair[1] for pair in mappedIDList]
        totalPostAtomList = list(postAtomObjectDict.keys())
        unfoundMissingPostAtoms = [atomID for atomID in totalPostAtomList if atomID not in mappedPostAtoms and atomID not in missingPostAtomList]
        missingPostAtomList.extend(unfoundMissingPostAtoms)

        missingPostAtomObjects = get_missing_atom_objects(missingPostAtomList, postAtomObjectDict)

        map_missing_atoms(missingPreAtomObjects, missingPostAtomObjects, mappedIDList, queue, inference)

        missingPreAtomList = update_missing_list(missingPreAtomList, mappedIDList, 0)
        missingPostAtomList = update_missing_list(missingPostAtomList, mappedIDList, 1)

        run_queue(queue, mappedIDList, preAtomObjectDict, postAtomObjectDict, missingPreAtomList, missingPostAtomList, elementDictList)
        logging.debug(f'循环 {timeoutCounter} 后的缺失前原子: {missingPreAtomList}') 

        if missingPreAtomCount == len(missingPreAtomList):
            inference = True
        else:
            inference = False

        timeoutCounter += 1

    if len(missingPreAtomList) > 0:
        print('错误: 缺失原子搜索超时。映射中将缺少原子。如果问题持续存在，请在GitHub上提交问题。')
        sys.exit()

    return mappedIDList
