#!/usr/bin/env python3
#####################################
#
# Filename : test_Superimposer.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Fri 18 May 2018 06:28:53 PM CEST
#
# Last Modified : Thu 23 Aug 2018 04:49:50 PM CEST
#
#####################################


import numpy as np
from pytest import approx


from distruct import config
from distruct import Superimposer
from distruct import Distructure
from distruct.tools.pdb import cull_atoms

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

    sup = Superimposer()
    sup.set_coords(x, y)
    # TODO is this really that bad??
    assert sup.rms == approx(0., abs=1e-2)
    return


def test_superimposer_atoms():

    from Bio.PDB.PDBParser import PDBParser

    code = '1ptq'
    fileName = testFilePath + code + '.pdb'

    fixedS = PDBParser().get_structure(code, fileName)
    movingS = PDBParser().get_structure(code, fileName)

    # TODO transform moving

    sup = Superimposer()
    sup.set_atoms(list(fixedS.get_atoms()), list(movingS.get_atoms()))

    assert sup.rms == approx(0.)
    return


def test_superimposer_structure():

    from Bio import SeqIO
    from Bio.PDB import PDBParser

    code = '1ptq'
    fileName = testFilePath + code + '.pdb'

    refStructure = PDBParser().get_structure(code, fileName)

    sequences = []
    with open(fileName, 'rU') as f:
        sequences = [r.seq for r in SeqIO.parse(f, "pdb-seqres")]
        pass

    ds = Distructure('test', sequences, [[r.get_id() for r in refStructure.get_residues() if r.get_id()[0] == ' ']])
    ds.generate_primary_contacts()
    ds.run()

    refStructure = PDBParser().get_structure(code, fileName)

    sup = Superimposer()
    sup.set_structures(refStructure, ds)
    return


def test_compare():
    """
    Compare the result of the diSTruct superimposer to the biopython one.
    """

    from Bio import SeqIO
    from Bio.PDB import Superimposer as BPSuperimposer
    from Bio.PDB import PDBParser

    from distruct.tools.pdb import get_contacts

    code = '1ptq'
    fileName = testFilePath + code + '.pdb'

    refStructure = PDBParser().get_structure(code, fileName)
    contacts = get_contacts(refStructure[0], cutOff=5., minSeqDist=0)

    sequences = []
    with open(fileName, 'rU') as f:
        sequences = [r.seq for r in SeqIO.parse(f, "pdb-seqres")]
        pass

    ds = Distructure(
            'test',
            sequences,
            [
                [r.get_id() for r in c if r.get_id()[0] == ' ']
                for c in refStructure[0]
            ]
    )
    ds.generate_primary_contacts()
    ds.set_tertiary_contacts(contacts)
    ds.run()

    refStructure = PDBParser().get_structure(code, fileName)
    tempStructure = ds.copy()

    refAtoms = list(cull_atoms(refStructure.get_atoms(), ds))
    resAtoms = list(cull_atoms(tempStructure.get_atoms(), refStructure))

    assert len(refAtoms) > 3
    assert len(refAtoms) == len(resAtoms)

    dssup = Superimposer()
    dssup.set_atoms(refAtoms, resAtoms)
    dsRMSD = dssup.rms

    bpsup = BPSuperimposer()
    bpsup.set_atoms(refAtoms, resAtoms)
    bpRMSD = bpsup.rms

    for atom in resAtoms:
        atom.set_coord(-1 * atom.get_coord())
        pass

    bpsup.set_atoms(refAtoms, resAtoms)
    if bpsup.rms < bpRMSD:
        bpRMSD = bpsup.rms
        pass

    assert dsRMSD == approx(bpRMSD)

    return
