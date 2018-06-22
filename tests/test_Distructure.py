#!/usr/bin/env python3
#####################################
#
# Filename : test_Distructure.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Mon 11 Jun 2018 04:04:59 PM CEST
#
# Last Modified : Fri 22 Jun 2018 04:53:38 PM CEST
#
#####################################

from Bio import SeqIO
from Bio.PDB.PDBParser import PDBParser

from MOBi import Distructure
from MOBi import data
from MOBi import config
from MOBi import Superimposer


testFilePath = config.data_path + "tests/"


def test_Distructure():

    code = "1ptq"
    fileName = testFilePath + code + '.pdb'

    sequences = list()
    with open(fileName, 'rU') as f:
        # for record in SeqIO.parse(f, "pdb-seqres"):
        #     print(record)
        #     sequences.append(record.seq)
        #     pass
        sequences = [r.seq for r in SeqIO.parse(f, "pdb-seqres")]
        pass

    # create distruct
    ds = Distructure(code, sequences)

    # create primary edges
    ds.generate_primary_edges()
    # generate coordinates
    ds.run()
    return
