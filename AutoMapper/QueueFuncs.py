##############################################################################
# Developed by: Matthew Bone
# Last Updated: 30/07/2021
# Updated by: Matthew Bone
#
# Contact Details:
# Bristol Composites Institute (BCI)
# Department of Aerospace Engineering - University of Bristol
# Queen's Building - University Walk
# Bristol, BS8 1TR
# U.K.
# Email - matthew.bone@bristol.ac.uk
#
# File Description:
# This is a collection of the queue control tools that facilitate the custom
# path search used in mapping.
##############################################################################

import logging
from collections import deque

# Classes and functions for search
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
    queueAtomObjects = []
    for pair in queueAtoms:
        preAtom = preAtomObjectDict[pair[0]]
        postAtom = postAtomObjectDict[pair[1]]
        queueAtomObjects.append([preAtom, postAtom])
    queue.add(queueAtomObjects)

def queue_bond_atoms(preAtomObjectDict, preBondingAtoms, postAtomObjectDict, postBondingAtoms, mappedIDList, queue):
    # Loop through bonding atoms, getting atom objects and adding them to queue and mapped list
    for index, preBondAtom in enumerate(preBondingAtoms):
        preAtomObject = preAtomObjectDict[preBondAtom]
        postAtomObject = postAtomObjectDict[postBondingAtoms[index]]
        queue.add([[preAtomObject, postAtomObject]])
        mappedIDList.append([preBondAtom, postBondingAtoms[index]])
        logging.debug(f'Pre: {preBondAtom}, Post: {postBondingAtoms[index]} found with user specified bond atom')

def run_queue(queue, mappedIDList, preAtomObjectDict, postAtomObjectDict, missingPreAtomList, missingPostAtomList, elementDictList):
    while not queue.empty():
        currentAtoms = queue.get()
        for mainIndex, atom in enumerate(currentAtoms):
            atom.check_mapped(mappedIDList, mainIndex, elementDictList[mainIndex])
        
        newMap, missingPreAtoms, missingPostAtoms, queueAtoms = currentAtoms[0].map_elements(currentAtoms[1], preAtomObjectDict, postAtomObjectDict)

        # Convert queue atoms to atom class objects and add to queue
        add_to_queue(queue, queueAtoms, preAtomObjectDict, postAtomObjectDict)

        # Extend missing lists
        missingPreAtomList.extend(missingPreAtoms)
        missingPostAtomList.extend(missingPostAtoms)

        # Add new pairs to mapped ID list
        mappedIDList.extend(newMap)