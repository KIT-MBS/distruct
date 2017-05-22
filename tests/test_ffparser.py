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
# Last Modified : Mon 22 May 2017 14:07:27 CEST
#
#####################################

import MOBi

testFilePath = MOBi.data_path + 'tests/'


# NOTE example lines taken from a (?) gromos ff

def test_read_atoms():
    line = "   CA   CH1     0.00000     1"
    result = MOBi.tools.ffparser.read_atoms(line)
    assert result['CA'] == 'CH1'
    pass


def test_read_bonds():
    line = "   CA    CB    gb_27   "
    result = MOBi.tools.ffparser.read_bonds(line)
    assert result[('CA', 'CB')] == 'gb_27'
    pass


def test_read_angles():
    line = "   CB    CA     C     ga_13   "
    result = MOBi.tools.ffparser.read_angles(line)
    assert result[('CB', 'CA', 'C')] == 'ga_13'
    pass


def test_read_impropers():
    # TODO test for different kinds
    line = "   CA     N     C    CB     gi_2    "
    result = MOBi.tools.ffparser.read_impropers(line)
    assert result[('CA', 'N', 'C', 'CB')] == 'gi_2'
    pass


def test_read_dihedrals():
    line = "    N    CA     C    +N     gd_45   "
    result = MOBi.tools.ffparser.read_dihedrals(line)
    assert result[('N', 'CA', 'C', '+N')] == 'gd_45'
    pass

# TODO parse force field parameter functios


def test_residue_topology():
    fileName = testFilePath + 'test.ff/aminoacids.rtp'
    result = MOBi.tools.ffparser.parse_residue_topology(fileName, buildingBlocks=['ALA'])
    # TODO
    assert 'ALA' in result
    pass


def test_parse_forcefield_params():
    fileName = testFilePath + 'test.ff/ffbonded.itp'
    params, macros = MOBi.tools.ffparser.parse_forcefield_params(fileName)
    # TODO put in something sensible for biophysics context
    assert ('CA', 'CB') in params['bondtypes']
    pass


def test_generate_chemical_primary_edge_database():
    testForceField = 'test'
    result = MOBi.tools.ffparser.generate_chemical_primary_edge_database(testForceField, topPath=testFilePath, buildingBlocks=['ALA'])

    assert 'ALA' in result
    for residue in result:
        assert 'CA' in result[residue]['atoms']
        assert ('CA', 'C') in result[residue]['bonds']
        assert ('CA', 'C', 'O') in result[residue]['angles']
        # TODO test for default improper dihedral in backbone
    pass
