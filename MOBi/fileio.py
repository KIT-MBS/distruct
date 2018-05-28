#!/usr/bin/env python3
#####################################
#
# Filename : fileio.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Tue 15 Aug 2017 11:19:36 AM CEST
#
# Last Modified : Mon 28 May 2018 10:07:16 AM CEST
#
#####################################

from lxml import etree as ET

# TODO handle alphabets

from . import config
defaultDataPath = config.data_path


def read_edge_node(node):
    result = {}

    # TODO conversion for vertices attrib
    # should this be a set?
    result[node.attrib['vertices']] = float(node.attrib['dstance'])

    return result


# NOTE this is not the simplest way to dict2xml but it produces more easily (human) readable files (i feel)(so they are better me-readable, really)
# def parse_primary_edge_database(databaseName, inDir=defaultDataPath, fileName=None):
def read_topology_file(databaseName, inDir=defaultDataPath, fileName=None):
    result = {}
    XMLTree = None
    if fileName == None:
        fileName = databaseName + '.xml'
        pass
    with open(inDir + fileName, 'r') as f:
        XMLTree = ET.parse(f)
        pass

    for buildingBlock in XMLTree.getroot():
        # NOTE overwrites if there are duplicate building block entries in the file
        # if one was added by hand at the end, that one will be used?
        result[buildingBlock.tag] = {}
        for vertices in buildingBlock.findall('vertices'):
            if vertices.tag not in result:
                result[buildingBlock.tag][vertices.tag] = set()
                pass
            for atom in vertices.findall('atom'):
                result[buildingBlock.tag]['vertices'].add(atom.attrib['name'])
                pass
            pass
        for child in buildingBlock:
            if child.tag != 'vertices':
                if child.tag not in result[buildingBlock.tag]:
                    result[buildingBlock.tag][child.tag] = {}
                    pass
                for edge in child:
                    edgeTuple = tuple(edge.attrib['vertices'].strip('()').split(','))
                    edgeTuple = tuple(x.strip(" \'") for x in edgeTuple)
                    result[buildingBlock.tag][child.tag][edgeTuple] = float(edge.attrib['distance'])
                    pass
                pass
            pass
        # TODO iterate over all other tags
        # read bonds and pseudobonds
        # for bonds in buildingBlock.findall('bondsEdges'):
        #     if bonds.tag not in result[buildingBlock.tag]:
        #         result[buildingBlock.tag][bonds.tag] = {}
        #     for edge in bonds.findall('edge'):
        #         # TODO
        #         result[buildingBlock.tag][atom.attrib['vertices']] = float(atom.attrib['distance'])
        #         pass
        #     pass
        # # TODO
        # for pseudoBonds in buildingBlock.findall('pseudobonds'):
        #     if pseudoBonds.tag not in result[buildingBlock.tag]:
        #         result[buildingBlock.tag][pseudoBonds.tag] = {}
        #     for edge in pseudoBonds.findall('edge'):
        #         # TODO
        #         pass
        #     pass

        pass

    # TODO check if test fails correctly
    return result


# NOTE probably a good idea to name the database something like forcefield_buildingblocks.xml
# TODO get the sorting right: vertices, bondEdges, angleEdges, improperEdges
def write_primary_edge_database(
        database,
        databaseName,
        buildingBlocks=[],
        outDir='./',
        fileName=None):
    if not fileName:
        fileName = databaseName + '.xml'
        pass
    if len(buildingBlocks) == 0:
        buildingBlocks = list(database.keys())
        pass

    root = ET.Element(databaseName)
    XMLTree = ET.ElementTree(root)

    for buildingBlock in sorted(database):
        assert buildingBlock in database
        assert 'vertices' in database[buildingBlock]
        bbElement = ET.SubElement(root, buildingBlock)
        atomsElement = ET.SubElement(bbElement, 'vertices')
        for entry in sorted(database[buildingBlock]['vertices']):
            ET.SubElement(atomsElement, 'atom', attrib={'name': entry})
            pass
        directives = sorted(database[buildingBlock].keys())
        directives.remove('vertices')
        for directive in directives:
            directiveElement = ET.SubElement(bbElement, directive)
            for entry in sorted(database[buildingBlock][directive]):
                # TODO is three significant places adequate for all applications?
                distance = None
                if isinstance(database[buildingBlock][directive][entry], float):
                    distance = '{0:.5f}'.format(database[buildingBlock][directive][entry])
                else:
                    # distance = str(database[buildingBlock][directive][entry])
                    raise
                    pass
                ET.SubElement(directiveElement, 'edge', attrib={'vertices': str((entry)), 'distance': distance})
                pass
            pass
        pass

    # TODO add xml version line thingy
    XMLTree.write(outDir + fileName, pretty_print=True)

    return
