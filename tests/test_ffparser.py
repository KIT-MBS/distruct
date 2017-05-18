#!/usr/bin/env python3
#####################################
#
# Filename : test_ffparser.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Thu 11 May 2017 14:12:24 CEST
#
# Last Modified : Thu 18 May 2017 18:20:53 CEST
#
#####################################

import MOBi

testFilePath = MOBi.data_path + 'tests/'
# TODO find the gromacs data path automatically
# should be $GMXLIB/share/gromacs/top/ but $GMXLIB is not defined
ffPath = '/usr/share/gromacs/top/'


def test_readAtoms():
    line = "   CA   CH1     0.00000     1"
    result = MOBi.tools.ffparser.readAtoms(line)
    assert result['CA'] == 'CH1'
    pass


def test_readBonds():
    line = "   CA    CB    gb_27   "
    result = MOBi.tools.ffparser.readBonds(line)
    assert result[('CA', 'CB')] == 'gb_27'
    pass


def test_readAngles():
    line = "   CB    CA     C     ga_13   "
    result = MOBi.tools.ffparser.readAngles(line)
    assert result[('CB', 'CA', 'C')] == 'ga_13'
    pass


def test_readImpropers():
    # TODO test for different kinds
    line = "   CA     N     C    CB     gi_2    "
    result = MOBi.tools.ffparser.readImpropers(line)
    assert result[('CA', 'N', 'C', 'CB')] == 'gi_2'
    pass


def test_readDihedrals():
    line = "    N    CA     C    +N     gd_45   "
    result = MOBi.tools.ffparser.readDihedrals(line)
    assert result[('N', 'CA', 'C', '+N')] == 'gd_45'
    pass

# TODO parse force field parameter functios


def test_parseResidueTopology():
    # TODO data packaging
    fileName = testFilePath + 'test.rtp'
    result = MOBi.tools.ffparser.parseResidueTopology(fileName, buildingBlocks=['ALA'])
    # TODO
    assert 'ALA' in result
    pass


# def test_parseForceFieldParams():
#     fileName = ""
#     result = MOBi.tools.ffparser.parseForceField(fileName)
#     # TODO
#     assert ('CA', 'CB') in result
#     pass
#
#
# def test_ffparse():
#     testForceField = 'amber99sb-ildn'
#     result = MOBi.tools.ffparser.ffparse(ffPath, testForceField)
#
#     assert 'ALA' in result
#     for residue in result:
#         # TODO may be dropped for more generic situations where not only AAs present
#         assert 'CA' in residue['atoms']
#         assert ('CA', 'C') in residue['bonds']
#         assert ('CA', 'C', 'O') in residue['angles']
#         # TODO test for default improper dihedral in backbone
#     pass
