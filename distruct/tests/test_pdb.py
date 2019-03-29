#!/usr/bin/env python3
#####################################
#
# Filename : test_pdb.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Mon 30 Apr 2018 05:13:44 PM CEST
#
# Last Modified : Fri 10 Aug 2018 06:47:28 PM CEST
#
#####################################

# TODO test termini
# TODO test RNA and protein and DNA and mixed complexes
# TODO test parse xray and nmr generated structures
# TODO test all functions


from Bio.PDB import PDBParser


from distruct.tools.pdb import get_contacts
from distruct import config
testFilePath = config.data_path + 'tests/'


def test_get_contacts():

    code = '1ptq'
    fileName = testFilePath + code + '.pdb'

    structure = PDBParser().get_structure(code, fileName)

    contacts = get_contacts(structure[0])

    for contact in contacts:
        fullID1, fullID2 = contact[0]
        entity1 = entity2 = structure[0]

        for id1, id2 in zip(fullID1[2:], fullID2[2:]):

            if entity1.get_level() == 'R' and entity2.get_level() == 'R':
                entity1 = entity1[id1[0]]
                entity2 = entity2[id2[0]]
            elif entity1.get_level() != 'R' and entity2.get_level() != 'R':
                entity1 = entity1[id1]
                entity2 = entity2[id2]
            else:
                raise ValueError("Level mismatch of entities in contact while\
                        iterating through fullID")
            pass

        # if fullID1[2] in structure[0] and fullID2 in structure[0]:
        #     chain1 = structure[0][fullID1[2]]
        #     chain2 = structure[0][fullID2[2]]
        #     if fullID1[3] in chain1 and fullID2[3] in chain2:
        #         residue1 = chain1[fullID1[3]]
        #         residue2 = chain2[fullID2[3]]
        #         if fullID1[4][0] in residue1 and fullID2[4][0] in residue2:
        #             pass
        #         pass
        #     pass
        # pass
    return
