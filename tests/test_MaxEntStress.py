#!/usr/bin/env python3
#####################################
#
# Filename : test_MaxEntStress.py
#
# Projectname :
#
# Author : Oskar Taubert
#
# Creation Date : Tue 17 Apr 2018 04:25:49 PM CEST
#
# Last Modified : Fri 04 May 2018 06:13:40 PM CEST
#
#####################################


# TODO remove this
import numpy as np

import Bio


import MOBi

testFilePath = MOBi.config.data_path + 'tests/'


def test_BMS():
    topDB = MOBi.fileio.parse_primary_edge_database("amber99sb-ildn_protein")
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

    resStructure = MOBi.tools.pdbparser.build_structure(code, sequences, topDB, offsets=offsets)

    # TODO beautify edge generator

    for c in resStructure[0]:
        chainPrimaryEdges, chainPrimaryDistances, chainPrimaryWeights = MOBi.tools.pdbparser.get_primary_edges(c, topDB, useStructureDistances = False)
        primaryEdges += chainPrimaryEdges
        primaryDistances += chainPrimaryDistances
        primaryWeights += chainPrimaryWeights
        pass

    tertiaryRefEdges, tertiaryDistances, tertiaryWeights, IDpairs = MOBi.tools.pdbparser.get_tertiary_edges(refStructure[0], 5., 2, getContacts=True)
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
    sup = Bio.PDB.Superimposer()
    RMSDStructure = MOBi.tools.pdbparser.read_PDB(code, fileName, topDB)
    for c in RMSDStructure[0]:
        for r in c:
            for a in r:
                coord = resStructure[0][c.get_id()][r.get_id()][a.get_id()].get_coord()
                a.set_coord(coord)
                pass
            pass
        pass
    sup.set_atoms(list(refStructure.get_atoms()), list(RMSDStructure.get_atoms()))
    RMSD = sup.rms

    from pytest import approx
    # TODO check duplicate edges and cutoffs
    assert RMSD == approx(0.)

    return
