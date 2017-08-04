#!/usr/bin/env python3
#####################################
#
# Filename : structurewriter.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Fri 04 Aug 2017 11:07:00 CEST
#
# Last Modified : Fri 04 Aug 2017 11:09:54 CEST
#
#####################################

from Bio import PDB


# TODO for large structures, split up in chains?
def write_solution(fileName, structure):
    io = PDB.PDBIO()
    io.set_structure(structure)
    io.save(fileName)
    return
