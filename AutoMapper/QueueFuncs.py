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
# 这是一个队列控制工具的集合，用于支持映射中使用的自定义路径搜索。
##############################################################################

import logging
from collections import deque

# 用于搜索的类和函数
class Queue:
    def __init__(self):
        self.elements = deque()
    
    def empty(self) -> bool:
        return not self.elements
    
    def add(self, x):
        for pair in x:
            self.elements.append(pair)

    def get(self):
        return self.elements.popleft()


def add_to_queue(queue, queueAtoms, preAtomObjectDict, postAtomObjectDict):
    """将原子对添加到队列中"""
    queueAtomObjects = []
    for pair in queueAtoms:
        preAtom = preAtomObjectDict[pair[0]]
        postAtom = postAtomObjectDict[pair[1]]
        queueAtomObjects.append([preAtom, postAtom])
    queue.add(queueAtomObjects)

def queue_bond_atoms(preAtomObjectDict, preBondingAtoms, postAtomObjectDict, postBondingAtoms, mappedIDList, queue):
    """将键合原子添加到队列和映射列表中"""
    for index, preBondAtom in enumerate(preBondingAtoms):
        preAtomObject = preAtomObjectDict[preBondAtom]
        postAtomObject = postAtomObjectDict[postBondingAtoms[index]]
        queue.add([[preAtomObject, postAtomObject]])
        mappedIDList.append([preBondAtom, postBondingAtoms[index]])
        logging.debug(f'前: {preBondAtom}, 后: {postBondingAtoms[index]} 通过用户指定的键合原子找到')

def run_queue(queue, mappedIDList, preAtomObjectDict, postAtomObjectDict, missingPreAtomList, missingPostAtomList, elementDictList):
    """运行队列处理"""
    while not queue.empty():
        currentAtoms = queue.get()
        for mainIndex, atom in enumerate(currentAtoms):
            atom.check_mapped(mappedIDList, mainIndex, elementDictList[mainIndex])
        
        newMap, missingPreAtoms, missingPostAtoms, queueAtoms = currentAtoms[0].map_elements(currentAtoms[1], preAtomObjectDict, postAtomObjectDict)

        # 将队列原子转换为原子类对象并添加到队列
        add_to_queue(queue, queueAtoms, preAtomObjectDict, postAtomObjectDict)

        # 扩展缺失列表
        missingPreAtomList.extend(missingPreAtoms)
        missingPostAtomList.extend(missingPostAtoms)

        # 将新对添加到映射ID列表
        mappedIDList.extend(newMap)