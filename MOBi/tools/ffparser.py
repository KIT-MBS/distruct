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
# Last Modified : Mon 22 May 2017 13:48:13 CEST
#
#####################################

import os

# import MOBi


# topologyPath = MOBi.gromacs_topology_path


def read_atoms(line):
    parts = line.split()
    name = parts[0]
    type = parts[1]
    return {name: type}


def read_bonds(line):
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


def read_angles(line):
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


def read_dihedrals(line):
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


def read_impropers(line):
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


# NOTE amber type ffs have n- and c- terminus residues, that are used here
# NOTE if buildingBlokcs is empty it considers all the ones it can find
def parse_residue_topology(residueTopologyFile, macros={}, buildingBlocks=[], ignoredDirectives={'exclusions', 'dihedrals'}):
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
                # expand macros here
                # TODO this is not very general, but should suffice for now
                if line.split()[-1] in macros:
                    line = ' '.join([line.split()[:-1]]) + ' ' + macros[line.split()[-1]]
                    pass
                rtpLines.append(line)
                pass
            pass
        pass

    # first line is expected to be bondedtypes directive in a valid .rtp file
    if rtpLines[0].split()[1] != 'bondedtypes':
        # TODO exceptions or just return error?
        raise

    knownDirectives = {
            'atoms': read_atoms,
            'bonds': read_bonds,
            'exclusions': None,
            'angles': read_angles,
            'dihedrals': read_dihedrals,
            'impropers': read_impropers,
            'bondedtypes': None}  # TODO use bondedtypes directive to find default values for function types (check in the ffs if useful)

    buildingBlock = None  # type of AA or NA usually
    context = None  # see knownDirectives
    for line in rtpLines:
        parts = line.split()
        if parts[0] == '[' and parts[2] == ']':
            directive = parts[1]
            if directive in buildingBlocks or (len(buildingBlocks) == 0 and directive not in knownDirectives):
                buildingBlock = directive
                context = None
                result[buildingBlock] = {}
            elif directive in ignoredDirectives:
                context = None
                pass
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

    for buildingBlock in buildingBlocks:
        if buildingBlock not in result:
            # TODO more helpful error messages / warnings
            raise
        pass

    return result


def read_bondtypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    func = parts[2]
    b0 = parts[3]

    key = (i, j)
    value = (func, b0)

    return {key: value}


def read_angletypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    func = parts[3]
    th0 = parts[4]

    key = (i, j, k)
    value = (func, th0)
    return {key: value}


def read_dihedraltypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    l = parts[3]
    func = parts[4]
    ph0 = parts[5]

    key = (i, j, k, l)
    value = (func, ph0)
    return {key: value}


def parse_forcefield_params(ffParamFile):

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
    knownDirectives = {'bondtypes': read_bondtypes, 'constrainttypes': None, 'angletypes': read_angletypes, 'dihedralTypes': read_dihedraltypes}
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


# TODO write xml
# def ffparse(ffname, ignore=[], generateAngles=True, generateDihedrals=False, topPath=topologyPath):
def generate_chemical_primary_edge_database(ffname, buildingBlocks=[], generateAngles=True, generateDihedrals=False, topPath=''):
    residueTopologieFiles = []
    forceFieldParamFiles = []
    if os.path.isdir(topPath):
        if os.path.isdir(topPath + ffname + '.ff'):
            residueTopologieFiles = [f for f in os.listdir(topPath + ffname + '.ff/') if f.endswith('.rtp')]
            forceFieldParamFiles.append(topPath + ffname + '.ff/ffbonded.itp')
            pass
        pass

    ffParams = None
    macros = None
    for forceFieldParamFile in forceFieldParamFiles:
        ffParams, macros = parse_forcefield_params(forceFieldParamFile)
        pass

    result = {}
    for residueTopologyFile in residueTopologieFiles:
        result.update(parse_residue_topology(topPath + ffname + '.ff/' + residueTopologyFile, macros))
        pass

    for residue in result:
        # keys in result are atoms, bonds, angles, impropers, dihedrals
        assert 'atoms' in result[residue]
        for bond, distance in result[residue]['bonds']:
            if not distance:
                ffAtomType1 = result[residue]['atoms'][bond[0]]
                ffAtomType2 = result[residue]['atoms'][bond[1]]
                distance = ffParams['bondtypes'][(ffAtomType1, ffAtomType2)][1]
                pass
            pass
        # TODO infer angles
        for entry in result[residue]['angles']:
            pass
        # # omega ? PRO special casing
        # for entry in result[residue]['dihedrals']:
        #     pass
        # infer more impropers ??
        for entry in result[residue]['impropers']:
            pass
        pass

    return result
