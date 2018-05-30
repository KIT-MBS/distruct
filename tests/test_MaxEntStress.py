#!/usr/bin/env python3
#####################################
#
# Filename : test_MaxEntStress.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Tue 17 Apr 2018 04:25:49 PM CEST
#
# Last Modified : Mon 28 May 2018 10:11:46 AM CEST
#
#####################################


# TODO remove this
import numpy as np

import Bio


import MOBi


testFilePath = MOBi.config.data_path + 'tests/'


def test_BMS():
    topDB = MOBi.data.defaultTopologyDB
    code = '1ptq'
    fileName = testFilePath + code + '.pdb'

    refStructure = MOBi.tools.pdbparser.read_PDB(code, fileName, topDB)

    # NOTE get sequence out of reference structure
    sequences = []
    offsets = []
    for c in refStructure[0]:
        sequences.append([])
        offset = None
        for r in c:
            if not offset:
                offset = r.get_id()[1]
                offsets.append(offset)
                pass
            sequences[-1].append(r.get_resname())
            pass
        pass

    primaryEdges = []
    primaryDistances = []
    primaryWeights = []

    # resStructure = MOBi.tools.pdbparser.build_structure(code, sequences, topDB, offsets=offsets)
    resStructure = MOBi.tools.pdbparser.read_PDB(code, fileName, topDB)

    # TODO beautify edge generator

    for c in resStructure[0]:
        chainPrimaryEdges, chainPrimaryDistances, chainPrimaryWeights = MOBi.tools.pdbparser.get_primary_edges(c, topDB, useStructureDistances = True)
        primaryEdges += chainPrimaryEdges
        primaryDistances += chainPrimaryDistances
        primaryWeights += chainPrimaryWeights
        pass

    tertiaryRefEdges, tertiaryDistances, tertiaryWeights, IDpairs = MOBi.tools.pdbparser.get_tertiary_edges(refStructure[0], 5., 5, getContacts=True)
    tertiaryEdges = MOBi.tools.pdbparser.translate_to_edges(IDpairs, resStructure[0])

    edges = primaryEdges + tertiaryEdges
    distances = primaryDistances + tertiaryDistances
    weights = primaryWeights + tertiaryWeights

    alpha = 1.
    q = 0.

    coordinates = MOBi.doublyWrappedMaxent(len(list(resStructure.get_atoms())), alpha, q, 300, edges, distances, weights)

    # TODO beautify structure extraction

    # save coord to structure object
    j = 0
    for c in resStructure[0]:
        for r in c:
            for a in r:
                a.set_coord(np.array(coordinates[j]))
                j += 1
                pass
            pass
        pass

    # TODO put in error measures in MOBi
    # TODO put in procrustes with mirror
    # TODO put in fixed starting conditions
    sup = MOBi.Superimposer()
    # RMSDStructure = MOBi.tools.pdbparser.read_PDB(code, fileName, topDB)
    # for c in RMSDStructure[0]:
    #     for r in c:
    #         for a in r:
    #             coord = resStructure[0][c.get_id()][r.get_id()][a.get_id()].get_coord()
    #             a.set_coord(coord)
    #             pass
    #         pass
    #     pass

    refAtoms = list(refStructure.get_atoms())
    resAtoms = list(resStructure.get_atoms())
    # RMSDAtoms = list(RMSDStructure.get_atoms())

    # sup.set_atoms(list(refStructure.get_atoms()), list(RMSDStructure.get_atoms()))

    sup.set_atoms(refAtoms, resAtoms)
    RMSD = sup.rms
    sup.apply(resStructure.get_atoms())

    # io = Bio.PDB.PDBIO()
    # io.set_structure(resStructure)
    # io.save('test.pdb')

    # TODO compare my RMSD to Biopython RMSD

    from pytest import approx
    # TODO check duplicate edges and cutoffs
    # TODO put in actual value
    # assert RMSD == approx(0.)
    assert RMSD < 1.0  # TODO find cause for deterioration

    return
