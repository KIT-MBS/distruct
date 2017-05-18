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
# Last Modified : Thu 18 May 2017 17:59:18 CEST
#
#####################################

import os

# import MOBi


# topologyPath = MOBi.gromacs_topology_path


def readAtoms(line):
    parts = line.split()
    name = parts[0]
    type = parts[1]
    return {name: type}


def readBonds(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    if len(parts) > 2:
        b0 = parts[2]
    else:
        b0 = None
        pass

    key = (i, j)
    value = b0
    return {key: value}


def readAngles(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    if len(parts) > 3:
        th0 = parts[3]
    else:
        th0 = None
        pass

    key = (i, j, k)
    value = th0
    return {key: value}


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

    key = (i, j, k, l)
    value = phi0
    return {key: value}


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

    key = (i, j, k, l)
    value = q0
    return {key: value}


# TODO directives as parameter
# TODO ignore or select molecule list
# NOTE amber type ffs have n- and c- terminus residues, that are used here
# NOTE this still saves not needed junk in the force field
# i.e. for some reason there is urea, some ions etc. in the aminoacids.rtp
def parseResidueTopology(residueTopologyFile, macros={}, buildingBlocks=[]):
    result = dict(zip(buildingBlocks, [{}] * len(buildingBlocks)))

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
            # if directive not in knownDirectives:
            if directive in buildingBlocks:
                buildingBlock = directive
                context = None
                result[buildingBlock] = {}
            elif buildingBlock:
                context = directive
                result[buildingBlock][context] = {}
                pass
        else:
            if buildingBlock and context:
                if knownDirectives[context]:
                    # result[buildingBlock][context].append(knownDirectives[context](line))
                    # key, value = knownDirectives[context](line)
                    result[buildingBlock][context].update(knownDirectives[context](line))
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

    key = (i, j)
    value = (func, b0)

    return {key, value}


def readAngleTypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    func = parts[3]
    th0 = parts[4]

    key = (i, j, k)
    value = (func, th0)
    return {key, value}


def readDihedralTypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    l = parts[3]
    func = parts[4]
    ph0 = parts[5]

    key = (i, j, k, l)
    value = (func, ph0)
    return {key, value}


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
                pass
            if len(line.strip()) > 0:
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
                result[directive] = {}
        else:
            if directive in knownDirectives:
                if knownDirectives[directive]:
                    # result[directive].append(knownDirectives[directive](line))
                    # key, value = knownDirectives[directive](line)
                    result[directive].update(knownDirectives[directive](line))
                    pass
            pass
        pass
    return result, macros


# TODO infer angles from bonds
# TODO add list of ignored molecules and/or list of considered molecules (probably better)
# TODO write xml
# TODO rename this
# def ffparse(ffname, ignore=[], generateAngles=True, generateDihedrals=False, topPath=topologyPath):
def ffparse(ffname, ignore=[], generateAngles=True, generateDihedrals=False, topPath=''):
    # TODO discover folder structure .rtp .itp and so on (look at pdb2gmx)
    residueTopologies = []
    forceFieldParamFiles = []
    if os.path.isdir(topPath):
        if os.path.isdir(topPath + ffname + '.ff'):
            residueTopologies = [f for f in os.listdir(topPath + ffname + '.ff/') if f.endswith('.rtp')]
            forceFieldParamFiles.append(topPath + ffname + '.ff/ffbonded.itp')
            pass
        pass

    ffParams = None
    for forceFieldParamFile in forceFieldParamFiles:
        ffParams, macros = parseForceFieldParams(forceFieldParamFile)
        pass

    result = {}
    for residueTopology in residueTopologies:
        result.update(parseResidueTopology(residueTopology, macros))
        pass

    # TODO check if angles present, if not infer angles from force field and present bonds
    # TODO only consider useful dihedrals (planar, with single minimum and so on)
    for residue in result:
        # keys in result are atoms, bonds, angles, impropers, dihedrals
        assert 'atoms' in result[residue]
        # dict from atom name (unique identifier in residue) to type (in ff)
        atoms = {}
        # TODO put macro stuff in topology parse
        for atomEntry in result[residue]['atoms']:
            atoms[atomEntry[0]] = atomEntry[1]
            pass
        for entry in result[residue]['bonds']:
            if entry[3] in macros:
                entry[3] = macros[entry[3]]
            else:
                # get parameter from ff
                atom1 = entry[1]
                atom2 = entry[2]
                atomtype1 = atoms[atom1]
                atomtype2 = atoms[atom2]
                entry[3] = ffParams['bondtypes'][(atomtype1, atomtype2)]
                pass

            pass
        # TODO infer angles
        for entry in result[residue]['angles']:
            pass
        # omega ?
        for entry in result[residue]['dihedrals']:
            pass
        # infer more impropers ??
        for entry in result[residue]['impropers']:
            pass
        pass

    return result
