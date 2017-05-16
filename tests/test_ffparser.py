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
# Last Modified : Tue 16 May 2017 10:28:00 CEST
#
#####################################

import MOBi

testFilePath = MOBi.data_path + 'test/'
# TODO find the gromacs data path automatically
# should be $GMXLIB/share/gromacs/top/ but $GMXLIB is not defined
ffPath = '/usr/share/gromacs/top/'


def test_readAtoms():
    line = "   CA   CH1     0.00000     1"
    result = MOBi.tools.readAtoms(line)
    assert result[0] == 'C'
    assert result[1] == 'C'
    pass


def test_readBonds():
    line = "   CA    CB    gb_27   "
    result = MOBi.tools.readBonds(line)
    assert result[0] == 'CA'
    assert result[1] == 'CB'
    assert result[2] == 'gb_27'
    pass


def test_readAngles():
    line = "   CB    CA     C     ga_13   "
    result = MOBi.tools.readAngles(line)
    assert result[0] == 'CB'
    assert result[1] == 'CA'
    assert result[2] == 'C'
    assert result[3] == 'ga_13'
    pass


def test_readImpropers():
    # TODO test for different kinds
    line = "   CA     N     C    CB     gi_2    "
    result = MOBi.tools.readImpropers(line)
    assert result[0] == 'CA'
    assert result[1] == 'N'
    assert result[2] == 'C'
    assert result[3] == 'CB'
    assert result[4] == 'gi_2'
    pass


def test_readDihedrals():
    line = "    N    CA     C    +N     gd_45   "
    result = MOBi.tools.readDihedrals(line)
    assert result[0] == 'N'
    assert result[1] == 'CA'
    assert result[2] == 'C'
    assert result[3] == '+N'
    assert result[4] == 'gd_45'
    pass


def test_parseResidueTopology():
    # TODO data packaging
    fileName = ""
    result = MOBi.tools.parseResidueTopology(fileName)
    # TODO
    assert 'ALA' in result
    pass


# def test_parseForceFieldParams():
#     fileName = ""
#     result = MOBi.tools.parseForceField(fileName)
#     # TODO
#     assert ('CA', 'CB') in result
#     pass
#
#
# def test_ffparse():
#     testForceField = 'amber99sb-ildn'
#     result = MOBi.tools.ffparse(ffPath, testForceField)
#
#     assert 'ALA' in result
#     for residue in result:
#         # TODO may be dropped for more generic situations where not only AAs present
#         assert 'CA' in residue['atoms']
#         assert ('CA', 'C') in residue['bonds']
#         assert ('CA', 'C', 'O') in residue['angles']
#         # TODO test for default improper dihedral in backbone
#     pass
