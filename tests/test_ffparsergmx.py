#!/usr/bin/env python3
#####################################
#
# Filename : test_ffparsergmx.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Thu 11 May 2017 14:12:24 CEST
#
# Last Modified : Mon 30 Jul 2018 02:40:57 PM CEST
#
#####################################

from math import pi
from math import sqrt
from pytest import approx

from distruct.tools import ffparsergmx

from distruct import config

testFilePath = config.data_path + 'tests/'


def test_read_atomtypes():
    line = "CA    CT           0.03370     3"
    result = ffparsergmx.read_atomtypes(line)
    assert result['CA'] == 'CT'
    return

def test_read_atoms():
    line = "CA    CT           0.03370     3"
    result = ffparsergmx.read_atoms(line)
    assert result == 'CA'
    return


def test_read_bonds():
    line = "CA    CB"
    result = ffparsergmx.read_bonds(line)
    assert result['CA', 'CB'] is None
    assert result['CB', 'CA'] is None
    return


def test_read_angles():
    line = "N   CA  C"
    result = ffparsergmx.read_angles(line)
    assert result['N', 'CA', 'C'] is None
    assert result['C', 'CA', 'N'] is None
    return


def test_read_impropers():
    line = "-C    CA     N     H"
    result = ffparsergmx.read_impropers(line)
    assert result == {('-C', 'CA', 'N', 'H'): None}
    return


def test_parse_residue_topology():
    fileName = testFilePath + 'test.ff/buildingblocks.rtp'
    macros = {'A1-4_improper': '0.0'}
    ignoredDirectives = {'bondedtypes'}
    result = ffparsergmx.parse_residue_topology(
            fileName,
            macros,
            ignoredDirectives = ignoredDirectives)
    assert result['BB1']['atoms'] == {'A1': 'A1', 'A2': 'A2', 'A3': 'A3', 'A4': 'A4'}
    assert result['BB1']['bonds'] == {('A1', 'A2'): None, ('A2', 'A3'): None, ('A3', 'A4'): None, ('A2', 'A1'): None, ('A3', 'A2'): None, ('A4', 'A3'): None}
    assert result['BB1']['angles'] == {('A1', 'A2', 'A3'): None, ('A2', 'A3', 'A4'): None, ('A3', 'A2', 'A1'): None, ('A4', 'A3', 'A2'): None}
    assert result['BB1']['impropers'] == {('A1', 'A2', 'A3', 'A4'): 0.0}
    assert 'dihedrals' not in result['BB1']
    return


def test_read_bondtypes():
    line = "CT  N     1   0.14490  282001.6"
    result = ffparsergmx.read_bondtypes(line)
    assert result[('CT', 'N')] == approx(1.449)
    assert result[('N', 'CT')] == approx(1.449)
    return


def test_read_angletypes():
    line = "C   N   CT            1 121.900   418.400"
    result = ffparsergmx.read_angletypes(line)
    assert result[('C', 'N', 'CT')] == approx(121.9 * pi/180.)
    assert result[('CT', 'N', 'C')] == approx(121.9 * pi/180.)
    return


def test_read_dihedraltypes():
    line = "CT  O   C  OH      4     180.00   4.23154   2    "
    result = ffparsergmx.read_dihedraltypes(line)
    assert result[('CT', 'O', 'C', 'OH')] == approx(180.0 * pi/180.)
    line = "C   N   CT  C      9       0.00  1.23523    3    "
    result = ffparsergmx.read_dihedraltypes(line)
    assert result == {}
    return


def test_parse_forcefield_params():
    fileName = testFilePath + 'test.ff/ffbonded.itp'
    params, macros = ffparsergmx.parse_forcefield_params(fileName)
    assert macros == {'A1-4_improper': '0.0'}

    assert params['bondtypes'][('A1', 'A2')] == approx(1.2)
    assert params['bondtypes'][('A2', 'A3')] == approx(1.0)
    assert params['bondtypes'][('A3', 'A4')] == approx(1.1)
    assert params['angletypes'][('A1', 'A2', 'A3')] == approx(120. * pi/180.)
    assert params['angletypes'][('A2', 'A3', 'A4')] == approx(115. * pi/180.)
    assert params['dihedraltypes'][('A1', 'A2', 'A3', 'A4')] == approx(180. * pi/180.)

    # TODO test for additional entries
    return


def test_infer_angles():
    # methane
    atoms = {'C': 'CT', 'H1': 'H', 'H2': 'H', 'H3': 'H', 'H4': 'H'}
    bonds = {('C', 'H1'): 1.087, ('C', 'H2'): 1.087, ('C', 'H3'): 1.087, ('C', 'H4'): 1.087}
    angleTypes = {('H', 'CT', 'H'): 109.5}
    result = ffparsergmx.infer_angles(atoms, bonds, angleTypes)
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


def test_translate_bonds_to_edges():
    bonds = {('A1', 'A2'): None, ('A2', 'A3'): 1.2}
    atomTypes = {'A1': 'AT1', 'A2': 'AT2', 'A3': 'AT3'}
    bondTypes = {('AT1', 'AT2'): 1.1}
    result = ffparsergmx.translate_bonds_to_edges(bonds, atomTypes, bondTypes)
    assert result == {('A1', 'A2'): 1.1, ('A2', 'A3'): 1.2}
    return


def test_translate_angles_to_edges():
    angles = {('A1', 'A2', 'A3'): 120. * pi/180., ('A2', 'A3', 'A4'): None}
    atomTypes = {'A1': 'AT1', 'A2': 'AT2', 'A3': 'AT3', 'A4': 'AT4'}
    bondTypes = {('AT1', 'AT2'): 1.2, ('AT2', 'AT3'): 1.0, ('AT3', 'AT4'): 1.1, }
    angleTypes = {('AT2', 'AT3', 'AT4'): 115. * pi/180.}
    result = ffparsergmx.translate_angles_to_edges(angles, atomTypes, bondTypes, angleTypes)

    assert result[('A1', 'A3')] == approx(1.90787884028338913, rel=1e-5)
    assert result[('A2', 'A4')] == approx(1.7719368430701863, rel=1e-5)
    return


def test_translate_impropers_to_edges():
    impropers = {('A1', 'A2', 'A3', 'A4'): 0.}
    bondEdges = {('A1', 'A2'): 1., ('A2', 'A3'): 1., ('A3', 'A4'): 1.}
    angleEdges = {('A1', 'A3'): sqrt(2), ('A2', 'A4'): sqrt(2.)}
    atomTypes = {'A1': 'AT1', 'A2': 'AT2', 'A3': 'AT3', 'A4': 'AT4'}
    dihedralTypes = {('AT1', 'AT2', 'AT3', 'AT4'): 180.}

    result = ffparsergmx.translate_impropers_to_edges(impropers, angleEdges, bondEdges, atomTypes, dihedralTypes)

    assert result[('A1', 'A4')] == approx(1.)

    impropers = {('A1', 'A2', 'A3', 'A4'): None}
    result = ffparsergmx.translate_impropers_to_edges(impropers, angleEdges, bondEdges, atomTypes, dihedralTypes)
    assert result[('A1', 'A4')] == approx(sqrt(5.))
    return


def test_generate():
    ffname = 'test'
    from Bio import Alphabet
    alphabet = Alphabet.ProteinAlphabet()
    alphabet.size = 3
    alphabet.letters = ['BB1', 'BB2']
    inferAngles = True
    topPath = testFilePath

    result = ffparsergmx.generate(ffname, [alphabet], inferAngles, topPath=topPath)
    assert result['BB1']['vertices'] == [('A1', 'A'), ('A2', 'A'), ('A3', 'A'), ('A4', 'A')]

    assert result['BB1']['bondEdges'][('A1', 'A2')] == approx(1.2)
    assert result['BB1']['bondEdges'][('A2', 'A3')] == approx(1.0)
    assert result['BB1']['bondEdges'][('A3', 'A4')] == approx(1.1)
    assert result['BB1']['angleEdges'][('A1', 'A3')] == approx(1.90787884028338913, rel=1e-5)
    assert result['BB1']['angleEdges'][('A2', 'A4')] == approx(1.7719368430701863, rel=1e-5)
    assert result['BB1']['improperEdges']['A1', 'A4'] == approx(2.065313144262336)

    return
