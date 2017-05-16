#!/usr/bin/env python3
#####################################
#
# Filename : ffparser.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Thu 11 May 2017 10:54:56 CEST
#
# Last Modified : Tue 16 May 2017 16:20:54 CEST
#
#####################################

import os


def readAtoms(line):
    parts = line.split()
    name = parts[0]
    type = parts[1]
    return (name, type)


def readBonds(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    if len(parts) > 2:
        b0 = parts[3]
    else:
        b0 = None
        pass
    return (i, j, b0)


def readAngles(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    if len(parts) > 3:
        th0 = parts[4]
    else:
        th0 = None
        pass
    return (i, j, k, th0)


def readDihedrals(line):
    # TODO only add, if single minimum
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    l = parts[3]
    if len(parts) > 4:
        phi0 = parts[4]
    else:
        phi0 = None
        pass
    return (i, j, k, l, phi0)


def readImpropers(line):
    # TODO only add, if single minimum
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    l = parts[3]
    if len(parts) > 4:
        q0 = parts[4]
    else:
        q0 = None
        pass
    return (i, j, k, l, q0)


# TODO directives as parameter
# TODO ignore or select molecule list
# NOTE amber type ffs have n- and c- terminus residues, that are used here
# NOTE this still saves not needed junk in the force field
# i.e. for some reason there is urea, some ions etc. in the aminoacids.rtp
def parseResidueTopology(residueTopologyFile):
    result = {}

    # read lines from .rtp file and strip comments
    rtpLines = []
    with open(residueTopologyFile, 'r') as f:
        for line in f:
            line = line.strip()
            if ';' in line:
                line = line.split(';', 1)[0]
                pass
            if len(line) > 0:
                rtpLines.append(line)
                pass
            pass
        pass

    # first line is expected to be bondedtypes directive in a valid .rtp file
    if rtpLines[0].split()[1] != 'bondedtypes':
        # TODO exceptions or just return error?
        raise

    knownDirectives = {'atoms': readAtoms, 'bonds': readBonds, 'exclusions': None, 'angles': readAngles, 'dihedrals': readDihedrals, 'impropers': readImpropers, 'bondedtypes': None}

    # TODO use bondedtypes directive to find default values for function types

    buildingBlock = None  # type of AA or NA usually
    context = None  # see knownDirectives
    for line in rtpLines:
        parts = line.split()
        if parts[0] == '[' and parts[2] == ']':
            directive = parts[1]
            if directive not in knownDirectives:
                buildingBlock = directive
                result[buildingBlock] = {}
            elif buildingBlock:
                context = directive
                result[buildingBlock][context] = []
                pass
        else:
            if buildingBlock and context:
                if knownDirectives[context]:
                    result[buildingBlock][context].append(knownDirectives[context](line))
                    pass
                pass
            pass
        pass

    # TODO remove or ignore not needed molecules
    return result


def readBondTypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    func = parts[2]
    b0 = parts[3]
    return (i, j, func, b0)


def readAngleTypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    func = parts[3]
    th0 = parts[4]
    return (i, j, k, func, th0)


def readDihedralTypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    l = parts[3]
    func = parts[4]
    ph0 = parts[5]
    return (i, j, k, l, func, ph0)


def parseForceFieldParams(ffParamFile):

    itpLines = []

    macros = {}
    with open(ffParamFile, 'r') as f:
        for line in f:
            line = line.strip()
            if ';' in line:
                line = line.split(';', 1)[0]
                pass
            if '#' in line:
                # macros
                parts = line.split(None, 2)
                if parts[0] == '#define':
                    macros[parts[1]] = parts[2]
                    pass
                # ignore macro for the rest of this function
                line = line.split('#', 1)[0]
                print(line)
                pass
            if len(line) > 0:
                itpLines.append(line)
                pass
            pass
        pass

    result = {}
    knownDirectives = {'bondtypes': readBondTypes, 'constrainttypes': None, 'angletypes': readAngleTypes, 'dihedralTypes': readDihedralTypes}
    directive = None
    for line in itpLines:
        parts = line.split()
        if parts[0] == '[' and parts[2] == ']':
            directive = parts[1]
            if directive in knownDirectives and directive not in result:
                result[directive] = []
        else:
            if directive in knownDirectives:
                if knownDirectives[directive]:
                    result[directive].append(knownDirectives[directive](line))
                    pass
            pass
        pass
    return result, macros


# TODO infer angles from bonds
# TODO add list of ignored molecules and/or list of considered molecules (probably better)
# TODO write xml
# TODO rename this
def ffparse(ffname, path, ignore=[], generateAngles=True, generateDihedrals=False):
    # TODO discover folder structure .rtp .itp and so on (look at pdb2gmx)
    residueTopologies = []
    forceFieldParamFiles = []
    if os.path.isdir(path):
        if os.path.isdir(path + ffname + '.ff'):
            residueTopologies = [f for f in os.listdir(path + ffname + '.ff/') if f.endswith('.rtp')]
            forceFieldParamFiles.append(path + ffname + '.ff/ffbonded.itp')
            pass
        pass

    result = {}
    # TODO create function parseResidueTopology()
    for residueTopology in residueTopologies:
        result.update(parseResidueTopology(residueTopology))
        pass

    # parse topology to find building blocks, atoms, bonds, angles and dihedrals

    # TODO check if angles present, if not infer angles from force field and present bonds
    # TODO parse forcefield to find parameters and replace None with those

    ffParams = None
    for forceFieldParamFile in forceFieldParamFiles:
        ffParams = parseForceFieldParams(forceFieldParamFile)
        pass

    print(ffParams)

    for key in result:
        for entry in result[key]:
            # TODO if param == None : read from ff
            # if param == string : read from ff
            # if param == float : leave as is
            # never put in numbers at this stage ? (problem with inter building block edges)
            pass

    return result
