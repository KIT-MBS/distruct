#!/usr/bin/env python3
#####################################
#
# Filename : test_databaseio.py
#
# Projectname :
#
# Author : Oskar Taubert
#
# Creation Date : Mon 22 May 2017 14:44:16 CEST
#
# Last Modified : Tue 08 Aug 2017 15:32:22 CEST
#
#####################################

import os

import MOBi

# TODO redo these

testFilePath = MOBi.config.data_path + 'tests/'
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
