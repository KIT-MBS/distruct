#!/usr/bin/env python3
#####################################
#
# Filename : pdbparser.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Thu 11 May 2017 16:35:51 CEST
#
# Last Modified : Mon 26 Jun 2017 18:04:11 CEST
#
#####################################

from Bio import PDB


def readPDB(PDBCode, fileName):
    # TODO

    parser = None

    if fileName.endswith('.pdb'):
        parser = PDB.PDBParser()
    elif fileName.endswith('.cif'):
        parser = PDB.MMCIFParser()
        pass
    structure = parser.get_structure(PDBCode, fileName)

    seq = []

    # for model in structure:
    #     for chain in model:
    #         for residue in chain:
    #             for atom in residue:
    #                 pass
    #             pass
    #         pass
    #     pass

    return structure, seq
