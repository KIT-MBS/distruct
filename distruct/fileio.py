#!/usr/bin/env python3
#####################################
#
# Filename : fileio.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Tue 15 Aug 2017 11:19:36 AM CEST
#
# Last Modified : Thu 13 Sep 2018 01:24:14 PM CEST
#
#####################################

from lxml import etree as ET

from . import config
from . import data
defaultDataPath = config.data_path


# NOTE this is not the simplest way to dict2xml but it produces more easily (human) readable files (i feel)(so they are better me-readable, really)
# def parse_primary_edge_database(databaseName, inDir=defaultDataPath, fileName=None):
def read_topology_database(databaseName, inDir=defaultDataPath, fileName=None):
    result = {}
    XMLTree = None
    if fileName == None:
        fileName = databaseName + '.xml'
        pass
    with open(inDir + fileName, 'r') as f:
        XMLTree = ET.parse(f)
        pass

    for buildingBlockNode in XMLTree.getroot():
        # NOTE overwrites if there are duplicate building block entries in the file
        # if one was added by hand at the end, that one will be used?

        if buildingBlockNode.tag == 'alphabets':
            result[buildingBlockNode.tag] = dict()
            for polymerTypeNode in buildingBlockNode:
                polymerType = polymerTypeNode.tag
                result[buildingBlockNode.tag][polymerType] = dict()
                for letterNode in polymerTypeNode:
                    assert len(letterNode.items()) == 1
                    letter, buildingBlock = letterNode.items()[0]
                    result[buildingBlockNode.tag][polymerType][letter] = buildingBlock
                    pass
                pass
            pass
        else:
            result[buildingBlockNode.tag] = dict()
            for vertices in buildingBlockNode.findall('vertices'):
                if vertices.tag not in result:
                    result[buildingBlockNode.tag][vertices.tag] = list()
                    pass
                # TODO check this preserves ordering
                for atom in vertices.findall('atom'):
                    result[buildingBlockNode.tag]['vertices'].append((atom.attrib['name'], atom.attrib['element']))
                    pass
                pass
            for child in buildingBlockNode:
                if child.tag != 'vertices':
                    if child.tag not in result[buildingBlockNode.tag]:
                        result[buildingBlockNode.tag][child.tag] = {}
                        pass
                    for edge in child:
                        s = edge.attrib['vertices'].strip('()')
                        s = s.split(',')
                        for i, x in enumerate(s):
                            if '"' in x:
                                s[i] = x.strip('" ')
                            else:
                                s[i] = x.strip("' ")
                                pass
                            pass
                        # edgeTuple = tuple(edge.attrib['vertices'].strip('()').split(','))
                        # edgeTuple = tuple(x.strip(" \'") for x in edgeTuple)
                        edgeTuple = tuple(s)
                        result[buildingBlockNode.tag][child.tag][edgeTuple] = float(edge.attrib['distance'])
                        pass
                    pass
                pass
            pass
        pass

    return result


def write_topology_database(
        database,
        databaseName,
        alphabets,
        outDir='./',
        fileName=None
    ):
    if fileName is None:
        fileName = databaseName + '.xml'
        pass

    buildingBlocks = list()
    for a in alphabets:
        polymerType = data.polymer_type(a)
        buildingBlocks += [database['alphabets'][polymerType][letter] for letter in a.letters]
        pass

    root = ET.Element(databaseName)
    XMLTree = ET.ElementTree(root)

    if 'alphabets' in database:
        aElement = ET.SubElement(root, 'alphabets')
        for polymerType in database['alphabets']:
            pTElement = ET.SubElement(aElement, polymerType)
            for letter in database['alphabets'][polymerType]:
                buildingBlock = database['alphabets'][polymerType][letter]
                ET.SubElement(pTElement, 'letter', attrib={letter: database['alphabets'][polymerType][letter]})
                pass
            pass
        pass


    for buildingBlock in buildingBlocks:
        assert 'vertices' in database[buildingBlock]
        bbElement = ET.SubElement(root, buildingBlock)
        atomsElement = ET.SubElement(bbElement, 'vertices')
        for entry in database[buildingBlock]['vertices']:
            ET.SubElement(atomsElement, 'atom', attrib={'name': entry[0], 'element': entry[1]})
            pass
        directives = ['bondEdges', 'angleEdges', 'improperEdges', 'dihedralEdges']
        for directive in directives:

            directiveElement = ET.SubElement(bbElement, directive)

            if directive not in database[buildingBlock]:
                continue

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
