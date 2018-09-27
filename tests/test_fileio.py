#!/usr/bin/env python3
#####################################
#
# Filename : test_databaseio.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Mon 22 May 2017 14:44:16 CEST
#
# Last Modified : Thu 27 Sep 2018 07:31:31 PM CEST
#
#####################################

import os

from pytest import approx

from distruct.fileio import read_topology_database, write_topology_database
from distruct.tools.ffparsergmx import generate

from distruct import config

testFilePath = config.data_path + 'tests/'
#
#
# def test_parse_primary_edge_database():
#     databaseName = 'test'
#     referenceForceField = 'test'
#
#     database = MOBi.tools.databaseio.parse_primary_edge_database(databaseName, dataDir=testFilePath, fileName="readTest.xml")
#
#     referenceDatabase = MOBi.tools.ffparsergmx.generate_chemical_primary_edge_database(referenceForceField, topPath=testFilePath, buildingBlocks=['ALA'])
#
#     assert database == referenceDatabase
#     pass
#
#
# def test_write_primary_edge_database():
#
#     testFileName = "writeTest.xml"
#     testName = "test"
#
#     database = MOBi.tools.ffparsergmx.generate_chemical_primary_edge_database(testName, topPath=testFilePath, buildingBlocks=['ALA'])
#
#     referenceFileName = "testReference.xml"
#     MOBi.tools.databaseio.write_primary_edge_database(database, testName, dataDir=testFilePath, fileName=testFileName)
#
#     # TODO more efficient comparison? (i.e. store a hash)
#     with open(testFilePath + referenceFileName, 'r') as rf, open(testFilePath + testFileName, 'r') as tf:
#         for l1, l2 in zip(rf, tf):
#             assert l1 == l2
#             pass
#         pass
#     os.remove(testFilePath + testFileName)
#     return

def test_fileio():
    ffName = 'test'
    from Bio import Alphabet
    alphabet = Alphabet.ProteinAlphabet()
    alphabet.size = 3
    alphabet.letters = ['BB1', 'BB2']
    inferAngles = True
    topPath = testFilePath

    testDB = generate(ffName, [alphabet], inferAngles, topPath=topPath)

    write_topology_database(testDB, 'test', [alphabet], outDir=testFilePath)

    result = read_topology_database('test', inDir=testFilePath)
    os.remove(testFilePath + 'test.xml')
    assert result['BB1']['vertices'] == [('A1', 'A'), ('A2', 'A'), ('A3', 'A'), ('A4', 'A')]

    assert result['BB1']['bondEdges'][('A1', 'A2')] == approx(1.2)
    assert result['BB1']['bondEdges'][('A2', 'A3')] == approx(1.0)
    assert result['BB1']['bondEdges'][('A3', 'A4')] == approx(1.1)
    assert result['BB1']['angleEdges'][('A1', 'A3')] == approx(1.90787884028338913, rel=1e-5)
    assert result['BB1']['angleEdges'][('A2', 'A4')] == approx(1.7719368430701863, rel=1e-5)
    assert result['BB1']['improperEdges'][('A1', 'A4')] == approx(2.065313144262336, rel=1e-5)


    return


def test_quotationMarks():
    ffName = "test"
    from Bio import Alphabet
    alphabet = Alphabet.RNAAlphabet()
    alphabet.size = 1
    alphabet.letters = ['A']

    testDB = generate(ffName,  [alphabet], True, topPath=testFilePath)

    write_topology_database(testDB, 'test', [alphabet], outDir=testFilePath)

    result = read_topology_database('test', inDir=testFilePath)

    os.remove(testFilePath + 'test.xml')
    print(result['A']['vertices'])
    print(result['A']['bondEdges'])

    assert result['A']['vertices'] == [('A', 'A'), ("A1'", 'A'), ("A2'", 'A'), ("A1'1", 'A')]

    assert result['A']['bondEdges'][('A', "A1'")] == approx(1.2)
    assert result['A']['bondEdges'][("A1'", "A2'")] == approx(1.0)
    assert result['A']['bondEdges'][("A1'1", "A2'")] == approx(1.1)

    return
