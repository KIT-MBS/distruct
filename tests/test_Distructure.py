#!/usr/bin/env python3
#####################################
#
# Filename : test_Distructure.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Mon 11 Jun 2018 04:04:59 PM CEST
#
# Last Modified : Thu 28 Mar 2019 11:14:02 AM CET
#
#####################################

from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser

from distruct import Distructure
from distruct import config


testFilePath = config.data_path + "tests/"


def test_Distructure():

    code = "1ptq"
    fileName = testFilePath + code + '.pdb'

    sequences = list()
    with open(fileName, 'r') as f:
        # for record in SeqIO.parse(f, "pdb-seqres"):
        #     print(record)
        #     sequences.append(record.seq)
        #     pass
        sequences = [r.seq for r in SeqIO.parse(f, "pdb-seqres")]
        pass

    # create distruct
    ds = Distructure(code, sequences)

    # create primary contacts
    ds.generate_primary_contacts()
    # generate coordinates
    ds.run()
    return

if __name__ == '__main__':
    test_Distructure()
