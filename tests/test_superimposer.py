#!/usr/bin/env python3
#####################################
#
# Filename : test_Superimposer.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Fri 18 May 2018 06:28:53 PM CEST
#
# Last Modified : Wed 20 Jun 2018 03:15:18 PM CEST
#
#####################################


import numpy as np
from pytest import approx


from MOBi import config

testFilePath = config.data_path + "tests/"


def test_superimposer():
    x = np.array([[51.65, -1.90, 50.07],
         [50.40, -1.23, 50.65],
         [50.68, -0.04, 51.54],
         [50.22, -0.02, 52.85]], 'f')

    y = np.array([[51.30, -2.99, 46.54],
         [51.09, -1.88, 47.58],
         [52.36, -1.20, 48.03],
         [52.71, -1.18, 49.38]], 'f')
    # x = np.array([[-1., 0., 5.],
    #              [1., 0., 5.]], 'f')
    # y = np.array([[0., 1., 0.],
    #              [0., -1., 0.]], 'f')

    sup = MOBi.Superimposer()
    sup.set_coords(x, y)
    # TODO is this really that bad??
    assert sup.rms == approx(0., abs=1e-2)
    return


def test_superimposer_molecule():

    from Bio.PDB.PDBParser import PDBParser

    code = '1ptq'
    fileName = testFilePath + code + '.pdb'

    fixedS = PDBParser().get_structure(code, fileName)
    movingS = PDBParser().get_structure(code, fileName)

    # TODO transform moving

    sup = MOBi.Superimposer()
    sup.set_coords(fixed, moving)

    assert sup.rms == approx(0.)
    return


def test_compare():

    # code = '1ptq'
    # fileName = testFilePath + code + '.pdb'

    # refStructure = MOBi.tools.pdbparser.read_PDB(code, fileName, topDB)

    # resStructure = MOBi.tools.pdbparser.read_PDB(code, fileName, topDB)

    # TODO compare "my" superimposer to the one implemented in biopython

    assert 0

    return
