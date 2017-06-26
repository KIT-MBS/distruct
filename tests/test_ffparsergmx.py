#!/usr/bin/env python3
#####################################
#
# Filename : test_ffparsergmx.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Thu 11 May 2017 14:12:24 CEST
#
# Last Modified : Fri 23 Jun 2017 12:59:47 CEST
#
#####################################

import math

import MOBi

testFilePath = MOBi.config.data_path + 'tests/'

# TODO test on gromos and amber ffs


def test_read_atoms():
    line = "CA    CT           0.03370     3"
    result = MOBi.tools.ffparsergmx.read_atoms(line)
    assert result['CA'] == 'CT'
    pass


def test_read_bonds():
    line = "CA    CB"
    result = MOBi.tools.ffparsergmx.read_bonds(line)
    assert result['CA', 'CB'] is None
    assert result['CB', 'CA'] is None
    pass


def test_read_angles():
    line = "N   CA  C"
    result = MOBi.tools.ffparsergmx.read_angles(line)
    assert result['N', 'CA', 'C'] is None
    assert result['C', 'CA', 'N'] is None
    pass


def test_read_impropers():
    line = "-C    CA     N     H"
    result = MOBi.tools.ffparsergmx.read_impropers(line)
    assert result == {('-C', 'CA', 'N', 'H'): None}
    pass


# def test_read_dihedrals():
#     # TODO
#     assert 1
#     pass


def test_parse_residue_topology():
    fileName = testFilePath + 'test.ff/buildingblocks.rtp'
    macros = {'A1-4_improper': '0.0'}
    buildingBlocks = ['BB']
    ignoredDirectives = ['bondedtypes']
    result = MOBi.tools.ffparsergmx.parse_residue_topology(
            fileName,
            macros,
            buildingBlocks,
            ignoredDirectives)
    assert result['BB']['atoms'] == {'A1': 'A1', 'A2': 'A2', 'A3': 'A3', 'A4': 'A4'}
    assert result['BB']['bonds'] == {('A1', 'A2'): None, ('A2', 'A3'): None, ('A3', 'A4'): None, ('A2', 'A1'): None, ('A3', 'A2'): None, ('A4', 'A3'): None}
    assert result['BB']['angles'] == {('A1', 'A2', 'A3'): None, ('A2', 'A3', 'A4'): None, ('A3', 'A2', 'A1'): None, ('A4', 'A3', 'A2'): None}
    assert result['BB']['impropers'] == {('A1', 'A2', 'A3', 'A4'): 0.0}
    assert 'dihedrals' not in result['BB']
    return


def test_read_bondtypes():
    line = "CT  N     1   0.14490  282001.6"
    result = MOBi.tools.ffparsergmx.read_bondtypes(line)
    assert math.isclose(result[('CT', 'N')], 1.449)
    assert math.isclose(result[('N', 'CT')], 1.449)
    return


def test_read_angletypes():
    line = "C   N   CT            1 121.900   418.400"
    result = MOBi.tools.ffparsergmx.read_angletypes(line)
    assert math.isclose(result[('C', 'N', 'CT')], 121.9)
    assert math.isclose(result[('CT', 'N', 'C')], 121.9)
    return


def test_read_dihedraltypes():
    line = "CT  O   C  OH      4     180.00   4.23154   2    "
    result = MOBi.tools.ffparsergmx.read_dihedraltypes(line)
    assert math.isclose(result[('CT', 'O', 'C', 'OH')], 180.00)
    line = "C   N   CT  C      9       0.00  1.23523    3    "
    result = MOBi.tools.ffparsergmx.read_dihedraltypes(line)
    assert result == {}
    return


def test_parse_forcefield_params():
    fileName = testFilePath + 'test.ff/ffbonded.itp'
    params, macros = MOBi.tools.ffparsergmx.parse_forcefield_params(fileName)
    assert macros == {'A1-4_improper': '0.0'}
    assert math.isclose(params['bondtypes'][('A1', 'A2')], 1.2)
    assert math.isclose(params['bondtypes'][('A2', 'A3')], 1.0)
    assert math.isclose(params['bondtypes'][('A3', 'A4')], 1.1)
    assert math.isclose(params['angletypes'][('A1', 'A2', 'A3')], 120.)
    assert math.isclose(params['angletypes'][('A2', 'A3', 'A4')], 115.)
    assert math.isclose(params['dihedraltypes'][('A1', 'A2', 'A3', 'A4')], 180.)
    # TODO test for additional entries
    return


def test_infer_angles():
    # methane
    # TODO proper floating point comparison and more general test
    atoms = {'C': 'CT', 'H1': 'H', 'H2': 'H', 'H3': 'H', 'H4': 'H'}
    bonds = {('C', 'H1'): 1.087, ('C', 'H2'): 1.087, ('C', 'H3'): 1.087, ('C', 'H4'): 1.087}
    angleTypes = {('H', 'CT', 'H'): 109.5}
    result = MOBi.tools.ffparsergmx.infer_angles(atoms, bonds, angleTypes)
    assert result == {
            ('H1', 'C', 'H2'): 109.5,
            ('H1', 'C', 'H3'): 109.5,
            ('H1', 'C', 'H4'): 109.5,
            ('H2', 'C', 'H3'): 109.5,
            ('H2', 'C', 'H4'): 109.5,
            ('H2', 'C', 'H1'): 109.5,
            ('H3', 'C', 'H4'): 109.5,
            ('H3', 'C', 'H1'): 109.5,
            ('H3', 'C', 'H2'): 109.5,
            ('H4', 'C', 'H1'): 109.5,
            ('H4', 'C', 'H2'): 109.5,
            ('H4', 'C', 'H3'): 109.5}
    return


def test_translate_atoms_to_vertices():
    atoms = {'C': 'CT', 'H1': 'H', 'H2': 'H', 'H3': 'H', 'H4': 'H'}
    result = MOBi.tools.ffparsergmx.translate_atoms_to_vertices(atoms)
    assert result == {'C', 'H1', 'H2', 'H3', 'H4'}
    return


def test_translate_bonds_to_edges():
    bonds = {('A1', 'A2'): None, ('A2', 'A3'): 1.2}
    atomTypes = {'A1': 'AT1', 'A2': 'AT2', 'A3': 'AT3'}
    bondTypes = {('AT1', 'AT2'): 1.1}
    result = MOBi.tools.ffparsergmx.translate_bonds_to_edges(bonds, atomTypes, bondTypes)
    assert result == {('A1', 'A2'): 1.1, ('A2', 'A3'): 1.2}
    return


def test_translate_angles_to_edges():
    angles = {('A1', 'A2', 'A3'): 120., ('A2', 'A3', 'A4'): None}
    bondEdges = {('A1', 'A2'): 1.2, ('A2', 'A3'): 1.0, ('A3', 'A4'): 1.1, }
    atomTypes = {'A1': 'AT1', 'A2': 'AT2', 'A3': 'AT3', 'A4': 'AT4'}
    angleTypes = {('AT2', 'AT3', 'AT4'): 115.}
    result = MOBi.tools.ffparsergmx.translate_angles_to_edges(angles, bondEdges, atomTypes, angleTypes)
    print(result)
    assert math.isclose(result[('A1', 'A3')], 1.90787884028338913, rel_tol=1e-5)
    assert math.isclose(result[('A2', 'A4')], 1.7719368430701863, rel_tol=1e-5)
    return


def test_translate_impropers_to_edges():
    impropers = {('A1', 'A2', 'A3', 'A4'): 0.}
    bondEdges = {('A1', 'A2'): 1., ('A2', 'A3'): 1., ('A3', 'A4'): 1.}
    angleEdges = {('A1', 'A3'): math.sqrt(2), ('A2', 'A4'): math.sqrt(2.)}
    atomTypes = {'A1': 'AT1', 'A2': 'AT2', 'A3': 'AT3', 'A4': 'AT4'}
    dihedralTypes = {('AT1', 'AT2', 'AT3', 'AT4'): 180.}

    result = MOBi.tools.ffparsergmx.translate_impropers_to_edges(impropers, angleEdges, bondEdges, atomTypes, dihedralTypes)

    assert math.isclose(result[('A1', 'A4')], 1.)

    impropers = {('A1', 'A2', 'A3', 'A4'): None}
    result = MOBi.tools.ffparsergmx.translate_impropers_to_edges(impropers, angleEdges, bondEdges, atomTypes, dihedralTypes)
    assert math.isclose(result[('A1', 'A4')], math.sqrt(5.))
    return


def test_generate_chemical_primary_edge_database():
    ffname = 'test'
    buildingBlocks = ['BB']
    inferAngles = True
    topPath = testFilePath

    result = MOBi.tools.ffparsergmx.generate_chemical_primary_edge_database(ffname, buildingBlocks, inferAngles, topPath=topPath)
    assert result['BB']['vertices'] == {'A1', 'A2', 'A3', 'A4'}
    assert math.isclose(result['BB']['bondEdges'][('A1', 'A2')], 1.2)
    assert math.isclose(result['BB']['bondEdges'][('A2', 'A3')], 1.0)
    assert math.isclose(result['BB']['bondEdges'][('A3', 'A4')], 1.1)
    assert math.isclose(result['BB']['angleEdges'][('A1', 'A3')], 1.90787884028338913, rel_tol=1e-5)
    assert math.isclose(result['BB']['angleEdges'][('A2', 'A4')], 1.7719368430701863, rel_tol=1e-5)
    assert math.isclose(result['BB']['improperEdges']['A1', 'A4'], 2.065313144262336)
    return
