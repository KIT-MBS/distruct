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
# Last Modified : Thu 03 Aug 2017 11:21:29 CEST
#
#####################################

from Bio import PDB


def readPDB(PDBCode, fileName):

    parser = None

    if fileName.endswith('.pdb'):
        parser = PDB.PDBParser()
    elif fileName.endswith('.cif'):
        parser = PDB.MMCIFParser()
        pass
    structure = parser.get_structure(PDBCode, fileName)

    return structure
