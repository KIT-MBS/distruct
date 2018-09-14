#!/usr/bin/env python3
#####################################
#
# Filename : ffparsergmx.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Thu 11 May 2017 10:54:56 CEST
#
# Last Modified : Fri 14 Sep 2018 06:45:26 PM CEST
#
#####################################

import os

from Bio import Alphabet
from Bio.Data import IUPACData

from math import pi

from . import math as m

from .. import config
from .. import data

topologyPath = config.gromacs_topology_path

# TODO add omega to topologies
# TODO resolve ILE issue: ff: CD, PDB: CD1
# TODO similarly, some atoms in RNAs: OP2/O2P


def read_atomtypes(line):
    parts = line.split()
    name = parts[0]
    type = parts[1]
    return {name: type}


def read_atoms(line):
    parts = line.split()
    name = parts[0]
    return name


def read_bonds(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    b0 = None
    if len(parts) > 2:
        b0 = float(parts[2])
        pass

    return {(i, j): b0, (j, i): b0}


def read_angles(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    th0 = None
    if len(parts) > 3:
        th0 = float(parts[3]) * pi/180.
        pass

    return {(i, j, k): th0, (k, j, i): th0}


def read_dihedrals(line):  # NOTE proper dihedrals are not useful at this point
    # TODO only add, if single minimum
    # NOTE here function type and defaults matter
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    l = parts[3]
    if len(parts) > 4:
        phi0 = float(parts[4]) * pi/180.
    else:
        phi0 = None
        pass

    # NOTE order matters
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
        q0 = float(parts[4]) * pi/180.
    else:
        q0 = None
        pass

    # NOTE order matters
    # TODO more combinations are possible
    key = (i, j, k, l)
    value = q0
    return {key: value}


def parse_residue_topology(
        residueTopologyFile,
        macros={},
        ignoredDirectives={'bondedtypes', 'exclusions', 'dihedrals'},
        knownDirectives={'atoms': read_atomtypes, 'bonds': read_bonds, 'angles': read_angles, 'impropers': read_impropers}):
    """
    Parse force field parameters, translate residue names, cull, translate atom names.
    """

    result = {}

    for directive in ignoredDirectives:
        knownDirectives[directive] = None
        pass

    # read lines from .rtp file and strip comments
    rtpLines = []
    with open(residueTopologyFile, 'r') as f:
        for line in f:
            l = line.strip()
            if ';' in l:
                l = l.split(';', 1)[0]
                pass
            if len(l) > 0:
                # expand macros here
                # TODO this is not very general, but should suffice for now
                if l.split()[-1] in macros:
                    l = ' '.join(l.split()[:-1]) + ' ' + macros[l.split()[-1]]
                    pass
                rtpLines.append(l)
                pass
            pass
        pass

    # first line is expected to be bondedtypes directive in a valid .rtp file
    if rtpLines[0].split()[1] != 'bondedtypes':
        # TODO exceptions or just return error?
        raise

    buildingBlock = None  # type of AA or NA usually
    context = None  # see knownDirectives
    for line in rtpLines:
        parts = line.split()
        if parts[0] == '[' and parts[2] == ']':
            directive = parts[1]
            if directive not in knownDirectives:
                buildingBlock = directive
                context = None
                result[buildingBlock] = {}
            else:
                if directive not in ignoredDirectives:
                    if buildingBlock:
                        context = directive
                        # TODO this is getting too noodly
                        if context == 'atoms':
                            result[buildingBlock]['vertices'] = list()
                        result[buildingBlock][context] = {}
                        pass
                else:
                    context = None
                    pass
                pass
            pass
        else:
            if buildingBlock and context:
                if knownDirectives[context]:
                    result[buildingBlock][context].update(knownDirectives[context](line))
                    # TODO improve, kinda hacky
                    if context == 'atoms':
                        atomName = read_atoms(line)
                        atomType = result[buildingBlock]['atoms'][atomName][0]
                        element = atomType[0]
                        # TODO this does not work e.g. for the Mg in AMBER
                        if len(atomType) > 1:
                            if atomType[1].islower():
                                element += atomType[1].upper()
                            pass
                        result[buildingBlock]['vertices'].append((atomName, element))
                        pass
                    pass
                pass
            pass
        pass

    return result


def read_bondtypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    # func = parts[2]
    b0 = parts[3]

    key = (i, j)
    # key = frozenset([i, j])
    # value = (int(func), 10 * float(b0))
    value = 10 * float(b0)

    return {key: value, key[::-1]: value}


def read_angletypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    # func = parts[3]
    th0 = float(parts[4]) * pi/180.

    return {(i, j, k): th0 , (k, j, i): th0}


# TODO wildcards (X)
def read_dihedraltypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    l = parts[3]
    func = parts[4]
    ph0 = parts[5]
    if int(func) in [2, 4]:  # NOTE these denote impropers. proper dihedrals are not useful at this point

        # TODO check ordering stuff!!
        key = (i, j, k, l)
        # value = (func, ph0)
        value = float(ph0) * pi/180.
        return {key: value}
    else:
        return {}
    pass


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
                parts = line.split(maxsplit=2)
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
    knownDirectives = {'bondtypes': read_bondtypes, 'constrainttypes': None, 'angletypes': read_angletypes, 'dihedraltypes': read_dihedraltypes}
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
                    result[directive].update(knownDirectives[directive](line))
                    pass
            pass
        pass
    return result, macros


def infer_angles(atoms, bonds, angleTypes):
    result = {}

    # NOTE needed to get all angles between residues
    interResidueBonds = [b for b in bonds if '-' in b[0] or '+' in b[0]]

    addedBonds = []
    for bond in interResidueBonds:
        # TODO this should probably not happen for termini
        if '-' in bond[0]:
            bond = [bond[0].strip('-'), '+' + bond[1]]
            pass
        elif '+' in bond[0]:
            bond = [bond[0].strip('+'), '-' + bond[1]]
            pass
        else:
            raise
        if bond not in interResidueBonds:
            addedBonds.append(bond)
        pass

    for atom in atoms:  # NOTE both directions have to be in the dict for each bond
        localBonds = [b for b in bonds if b[0] == atom]

        for bond in addedBonds:
            if atom == bond[0]:
                localBonds.append(bond)
            pass

        while localBonds:
            bond1 = localBonds.pop()
            for bond2 in localBonds:
                i = bond1[1]
                j = atom
                k = bond2[1]

                ffAtomType1 = atoms[i.strip('+-')]
                ffAtomType2 = atoms[j.strip('+-')]
                ffAtomType3 = atoms[k.strip('+-')]

                angle = None
                if (ffAtomType1, ffAtomType2, ffAtomType3) in angleTypes:
                    angle = angleTypes[(ffAtomType1, ffAtomType2, ffAtomType3)]
                elif (ffAtomType3, ffAtomType2, ffAtomType1) in angleTypes:
                    angle = angleTypes[(ffAtomType3, ffAtomType2, ffAtomType1)]
                    pass
                else:
                    print("WARNING!!!: an angle type was not found in the force field parameters:\n", ([i, j, k]), ([ffAtomType1, ffAtomType2, ffAtomType3]))
                    pass
                if angle:
                    result[tuple([i, j, k])] = angle
                    result[tuple([k, j, i])] = angle
                    pass
                pass
            pass
        pass

    return result

# TODO add backbone dihedral omega


def translate_bonds_to_edges(bonds, atomTypes, bondTypes):
    result = {}
    for atomPair in bonds:
        vertices = tuple(sorted(atomPair))
        distance = bonds[atomPair]
        if distance is None:
            ffAtomType1 = atomTypes[atomPair[0].strip('+-')]
            ffAtomType2 = atomTypes[atomPair[1].strip('+-')]
            if (ffAtomType1, ffAtomType2) in bondTypes:
                distance = bondTypes[(ffAtomType1, ffAtomType2)]
            elif(ffAtomType2, ffAtomType1) in bondTypes:
                distance = bondTypes[(ffAtomType2, ffAtomType1)]
            else:
                print("WARNING!!!: a bond type was not found in the force field parameters:\n", atomPair, (ffAtomType1, ffAtomType2))
                continue
            pass
        assert isinstance(distance, float)
        result[vertices] = distance
        pass
    return result


def translate_angles_to_edges(angles, atomTypes, bondTypes, angleTypes):
    # TODO n-termini special casing
    result = {}
    for atomTriplet in angles:
        vertices = tuple(sorted((atomTriplet[0], atomTriplet[2])))
        angle = angles[atomTriplet]

        ffAtomType1 = atomTypes[atomTriplet[0].strip('+-')]
        ffAtomType2 = atomTypes[atomTriplet[1].strip('+-')]
        ffAtomType3 = atomTypes[atomTriplet[2].strip('+-')]

        if angle is None:
            if (ffAtomType1, ffAtomType2, ffAtomType3) in angleTypes:
                angle = angleTypes[(ffAtomType1, ffAtomType2, ffAtomType3)]
            elif (ffAtomType3, ffAtomType2, ffAtomType1) in angleTypes:
                angle = angleTypes[(ffAtomType3, ffAtomType2, ffAtomType1)]
            else:
                print("WARNING!!!: an angle type was not found in the force field parameters:\n", atomTriplet, (ffAtomType1, ffAtomType2, ffAtomType3))
                continue
            pass

        d_ij = None
        d_jk = None
        if (ffAtomType1, ffAtomType2) in bondTypes:
            d_ij = bondTypes[tuple(sorted((ffAtomType1, ffAtomType2)))]
        elif(ffAtomType2, ffAtomType1) in bondTypes:
            d_ij = bondTypes[tuple(sorted((ffAtomType2, ffAtomType1)))]
        else:
            print("WARNING!!!: a bond type was not found in the force field parameters:\n", (atomTriplet[0], atomTriplet[1]), (ffAtomType1, ffAtomType2))
            continue

        if (ffAtomType2, ffAtomType3) in bondTypes:
            d_jk = bondTypes[tuple(sorted((ffAtomType2, ffAtomType3)))]
        elif(ffAtomType3, ffAtomType2) in bondTypes:
            d_jk = bondTypes[tuple(sorted((ffAtomType3, ffAtomType2)))]
        else:
            print("WARNING!!!: a bond type was not found in the force field parameters:\n", (atomTriplet[1], atomTriplet[2]), (ffAtomType2, ffAtomType3))
            continue

        result[vertices] = m.dist_from_angle(angle, d_ij, d_jk)
        pass
    return result


# TODO switch to bondTypes and angleTypes
def translate_impropers_to_edges(impropers, angleEdges, bondEdges, atomTypes, dihedralTypes):
    result = {}
    for atomQuadruplet in impropers:
        # TODO order of residues is important!!!
        # NOTE it is assumed the dihedral angle is 0 or 180 degrees, all other 5 edges are known (3 bonds, 2 angles)

        # check missing edge
        pairs = [(x, y) for x in atomQuadruplet for y in atomQuadruplet if x < y]
        missingEdges = [(x, y) for (x, y) in pairs if (x, y) not in bondEdges if (x, y) not in angleEdges]

        # NOTE remove interresidue bond edges from the missing edges, that are in the next residue in sequence
        interResidueBonds = [b for b in bondEdges if '-' in b[0] or '+' in b[0] or '-' in b[1] or '+' in b[1]]
        addedEdges = []
        for bond in interResidueBonds:
            if '+' in bond[0]:
                addedEdges.append(tuple(sorted((bond[0].strip('+'), '-' + bond[1]))))
            elif '-' in bond[0]:
                addedEdges.append(tuple(sorted((bond[0].strip('-'), '+' + bond[1]))))
            elif '+' in bond[1]:
                addedEdges.append(tuple(sorted(('-' + bond[0], bond[1].strip('+')))))
            elif '-' in bond[1]:
                addedEdges.append(tuple(sorted(('+' + bond[0], bond[1].strip('-')))))
            else:
                raise
            pass

        for bond in addedEdges:
            if bond in missingEdges:
                missingEdges.remove(bond)
                pass
            pass

        if len(missingEdges) > 1:
            print("number of missing edges for improper dihedral is not equal to 1")
            print(atomQuadruplet)
            print(missingEdges)
            continue

        if len(missingEdges) == 0:
            continue

        vertices = tuple(sorted((atomQuadruplet[0], atomQuadruplet[3])))
        if missingEdges[0] != vertices:
            # TODO
            print("order of atoms in the dihedral can not be processed in this version. skipping...")
            print(atomQuadruplet)
            print(missingEdges)
            continue

        improper = impropers[atomQuadruplet]
        if improper is None:
            ffAtomType1 = atomTypes[atomQuadruplet[0]].strip('+-')
            ffAtomType2 = atomTypes[atomQuadruplet[1]].strip('+-')
            ffAtomType3 = atomTypes[atomQuadruplet[2]].strip('+-')
            ffAtomType4 = atomTypes[atomQuadruplet[3]].strip('+-')
            improper = dihedralTypes[(ffAtomType1, ffAtomType2, ffAtomType3, ffAtomType4)]

        ij = tuple(sorted([atomQuadruplet[0], atomQuadruplet[1]]))
        jk = tuple(sorted([atomQuadruplet[1], atomQuadruplet[2]]))
        kl = tuple(sorted([atomQuadruplet[2], atomQuadruplet[3]]))
        ik = tuple(sorted([atomQuadruplet[0], atomQuadruplet[2]]))
        jl = tuple(sorted([atomQuadruplet[1], atomQuadruplet[3]]))

        d_ij = bondEdges[ij]
        d_jk = bondEdges[jk]
        d_kl = bondEdges[kl]
        d_ik = angleEdges[ik]
        d_jl = angleEdges[jl]

        result[vertices] = m.dist_from_improper(improper, d_ij, d_jk, d_kl, d_ik, d_jl)

        pass
    return result

# TODO translate only needed  building blocks (avoids unnecessary warnings and errors)
# def translate(buildingBlockTopologies, ffParams, inferAngles=True, inferImpropers=False):
#     """
#     Translate building block topologies parsed from .rtp files to vertices and edges.
#     """
#     result = {}
#     for b in buildingBlockTopologies:
#         print(b)
#         result[b] = {}
#         result[b]['vertices'] = list()
#         result[b]['bondEdges'] = {}
#         result[b]['angleEdges'] = {}
#         result[b]['improperEdges'] = {}
#         result[b]['dihedralEdges'] = {}
#
#         # result[b]['vertices'] = translate_atoms_to_vertices(
#         #         buildingBlockTopologies[b]['atoms']
#         #         )
#         result[b]['vertices'] = buildingBlockTopologies[b]['vertices']
#
#         if 'bonds' in buildingBlockTopologies[b]:
#             result[b]['bondEdges'].update(
#                     translate_bonds_to_edges(
#                         buildingBlockTopologies[b]['bonds'],
#                         buildingBlockTopologies[b]['atoms'],
#                         ffParams['bondtypes']
#                     )
#             )
#
#             if inferAngles:
#                 inferredAngles = infer_angles(
#                         buildingBlockTopologies[b]['atoms'],
#                         buildingBlockTopologies[b]['bonds'],
#                         ffParams['angletypes'])
#                 result[b]['angleEdges'].update(
#                         translate_angles_to_edges(
#                             inferredAngles,
#                             buildingBlockTopologies[b]['atoms'],
#                             ffParams['bondtypes'],
#                             ffParams['angletypes']
#                         )
#                 )
#                 pass
#             pass
#             if 'angles' in buildingBlockTopologies[b]:
#                 result[b]['angleEdges'].update(translate_angles_to_edges(
#                     buildingBlockTopologies[b]['angles'],
#                     buildingBlockTopologies[b]['atoms'],
#                     ffParams['bondtypes'],
#                     ffParams['angletypes']))
#                 pass
#
#             if inferImpropers:
#                 raise NotImplementedError
#             if 'impropers' in buildingBlockTopologies[b]:
#                 result[b]['improperEdges'].update(
#                         translate_impropers_to_edges(
#                             buildingBlockTopologies[b]['impropers'],
#                             result[b]['angleEdges'],
#                             result[b]['bondEdges'],
#                             buildingBlockTopologies[b]['atoms'],
#                             ffParams['dihedraltypes'])
#                 )
#             pass
#         pass
#     return result


# def generate(
#         ffname,
#         alphabet,
#         inferAngles=True,
#         # inferImpropers=False,  # not implemented
#         topPath=topologyPath):
#
#     """
#     Generate a topology dictionary (residue -> ('vertices' -> list(vertices), edgeTypes -> (atomPair -> distance)))
#     """
#
#     # TODO ordering such that backbone atoms occur first in an ordered dict / order of ff is preserved?
#
#     residueTopologieFiles = []
#     forceFieldParamFiles = []
#
#     if os.path.isdir(topPath):
#         if os.path.isdir(topPath + ffname + '.ff'):
#             residueTopologieFiles = [f for f in os.listdir(topPath + ffname + '.ff/') if f.endswith('.rtp')]
#             forceFieldParamFiles.append(topPath + ffname + '.ff/ffbonded.itp')
#             pass
#         pass
#
#     ffParams = {}
#     macros = {}
#     for forceFieldParamFile in forceFieldParamFiles:
#         ffP, m = parse_forcefield_params(forceFieldParamFile)
#         ffParams.update(ffP)
#         macros.update(m)
#         pass
#
#     buildingBlockTopologies = {}
#     assert(len(residueTopologieFiles) > 0)
#     for residueTopologyFile in residueTopologieFiles:
#         buildingBlockTopologies.update(parse_residue_topology(topPath + ffname + '.ff/' + residueTopologyFile, macros))
#         pass
#
#     # TODO replace this with check on alphabet
#     # for buildingBlock in buildingBlocks:
#     #     if buildingBlock not in buildingBlockTopologies:
#     #         # TODO more helpful error messages / warnings
#     #         raise
#     #     pass
#
#     result = {}
#     # for buildingBlock in buildingBlockTopologies:
#     for letter in alphabet.letters:
#         result[letter] = {}
#
#         result[letter]['vertices'] = list()
#         result[letter]['bondEdges'] = {}
#         result[letter]['angleEdges'] = {}
#         result[letter]['improperEdges'] = {}
#         result[letter]['dihedralEdges'] = {}
#
#         # TODO translate letter to building block and vv
#         # TODO handle termini
#         buildingBlock = None
#
#         if letter in buildingBlockTopologies:
#             buildingBlock = letter
#         else:
#             if isinstance(alphabet, Alphabet.ProteinAlphabet):
#                 buildingBlock = Data.IUPACData.protein_letters_1to3_extended[letter]
#                 if buildingBlock not in buildingBlockTopologies:
#                     buildingBlock = buildingBlock.upper()
#                     if buildingBlock not in buildingBlockTopologies:
#                         # NOTE list is intended to catch meanings: [ beta, gamma, delta, epsilon, zeta, eta, protonated]
#                         # for postfix in ['B', 'G', 'D', 'E', 'Z', 'H', 'P']:
#                         # this is for protonated histidine
#                         # for positivly charged AAs the fully protonated is chosen,
#                         # for negatively charged ones the deprotonated ones
#                         for postfix in ['P']:
#                             if buildingBlock[:2] + postfix in buildingBlockTopologies:
#                                 buildingBlock = buildingBlock[:2] + postfix
#                                 break
#                             pass
#                         pass
#                     pass
#                 pass
#             elif isinstance(alphabet, Alphabet.ThreeLetterProtein):
#                 buildingBlock = letter
#                 for postfix in ['P']:
#                     if buildingBlock[:2] + postfix in buildingBlockTopologies:
#                         buildingBlock = buildingBlock[:2] + postfix
#                         break
#                     pass
#                 pass
#             elif isinstance(alphabet, Alphabet.RNAAlphabet):
#                 pass
#             elif isinstance(alphabet, Alphabet.DNAAlphabet):
#                 pass
#             else:
#                 # TODO error message
#                 assert(False)
#                 pass
#
#         result[letter]['vertices'].update(
#                 translate_atoms_to_vertices(
#                     buildingBlockTopologies[buildingBlock]['atoms']
#                 )
#         )
#
#         if 'bonds' in buildingBlockTopologies[buildingBlock]:
#             result[letter]['bondEdges'].update(
#                     translate_bonds_to_edges(
#                         buildingBlockTopologies[buildingBlock]['bonds'],
#                         buildingBlockTopologies[buildingBlock]['atoms'],
#                         ffParams['bondtypes']
#                     )
#             )
#
#             if inferAngles:
#                 inferredAngles = infer_angles(
#                         buildingBlockTopologies[buildingBlock]['atoms'],
#                         buildingBlockTopologies[buildingBlock]['bonds'],
#                         ffParams['angletypes'])
#                 result[letter]['angleEdges'].update(
#                         translate_angles_to_edges(
#                             inferredAngles,
#                             buildingBlockTopologies[buildingBlock]['atoms'],
#                             ffParams['bondtypes'],
#                             ffParams['angletypes']
#                         )
#                 )
#                 pass
#             pass
#             if 'angles' in buildingBlockTopologies[buildingBlock]:
#                 result[letter]['angleEdges'].update(translate_angles_to_edges(
#                     buildingBlockTopologies[buildingBlock]['angles'],
#                     buildingBlockTopologies[buildingBlock]['atoms'],
#                     ffParams['bondtypes'],
#                     ffParams['angletypes']))
#                 pass
#
#             # if inferImpropers:
#             #     inferredImpropers = infer_impropers()  # TODO
#             #     result[buildingBlock]['improperEdges'].update(translate_impropers_to_edges(
#             #         inferredImpropers,
#             #         result[letter]['angleEdges'],
#             #         result[letter]['bondEdges'],
#             #         buildingBlockTopologies[buildingBlock]['atoms'],
#             #         ffParams['dihedraltypes']))
#             #     pass
#
#             if 'impropers' in buildingBlockTopologies[buildingBlock]:
#                 result[letter]['improperEdges'].update(
#                         translate_impropers_to_edges(
#                             buildingBlockTopologies[buildingBlock]['impropers'],
#                             result[letter]['angleEdges'],
#                             result[letter]['bondEdges'],
#                             buildingBlockTopologies[buildingBlock]['atoms'],
#                             ffParams['dihedraltypes'])
#                 )
#             pass
#         pass
#     return result


# def generate(
#         ffName,
#         alphabets=[],
#         inferAngles=True,
#         inferImpropers=False,
#         topPath=topologyPath):
#     result = {}
#     ffDir = topPath + ffName + ".ff/"
#
#     residueTopologyFiles = [ffDir + f for f in os.listdir(ffDir) if f.endswith('.rtp')]
#     forceFieldParamFiles = [ffDir + "ffbonded.itp"]
#
#     params = {}
#     macros = {}
#     for f in forceFieldParamFiles:
#         p, m = parse_forcefield_params(f)
#         params.update(p)
#         macros.update(m)
#         pass
#
#     buildingBlockTopologies = {}
#     for f in residueTopologyFiles:
#         buildingBlockTopologies.update(parse_residue_topology(f, macros))
#         pass
#
#     result = translate(
#             buildingBlockTopologies,
#             params,
#             inferAngles=inferAngles,
#             inferImpropers=inferImpropers
#             )
#
#     # NOTE translate residue names
#     # NOTE uses the charged version of an AA (where applicable)
#     for letter in alphabet.letters:
#         buildingBlock = letter
#         if letter not in result:
#             if isinstance(alphabet, Alphabet.ProteinAlphabet):
#                 if alphabet.size == 1:
#                     buildingBlock = IUPACData.protein_letters_1to3_extended[letter]
#                     pass
#                 buildingBlock = buildingBlock.upper()
#                 if buildingBlock not in result:
#                     buildingBlock = buildingBlock[:2] + 'P'  # positively charged histidine
#                     pass
#                 if buildingBlock not in result:
#                     raise KeyError("Could not find " + letter + " in the residue database.")
#                 pass
#             elif isinstance(alphabet, Alphabet.DNAAlphabet):
#                 if alphabet.size == 1:
#                     buildingBlock = 'D' + letter
#                     pass
#                 if buildingBlock not in result:
#                     # NOTE should find all unambiguous residues
#                     raise KeyError("Could not find " + letter + " in the residue database.")
#                 raise
#             elif isinstance(alphabet, Alphabet.RNAAlphabet):
#                 if alphabet.size == 1:
#                     buildingBlock = 'R' + letter
#                     pass
#                 if buildingBlock not in result:
#                     # NOTE should find all unambiguous residues
#                     raise KeyError("Could not find " + letter + " in the residue database.")
#                 pass
#             else:
#                 raise Exception("Alphabet too generic.")
#             result[letter] = result[buildingBlock]
#             pass
#         # TODO translate atom names
#         # TODO reorder atoms so backbone comes first
#         pass
#     # NOTE remove unneeded building blocks
#     toDelete = set(result.keys()) - set(alphabet.letters)
#     for l in toDelete:
#         result.pop(l)
#         pass
#     # TODO gap character
#
#     return result


def translate(
        buildingBlock,
        buildingBlockTopologies,
        ffParams,
        inferAngles=True):
    """
    Translate a building block from a force field topology to a distruct topology.
    """

    result = {}

    result['vertices'] = list()
    result['bondEdges'] = {}
    result['angleEdges'] = {}
    result['improperEdges'] = {}
    result['dihedralEdges'] = {}

    result['vertices'] = buildingBlockTopologies[buildingBlock]['vertices']


    b = buildingBlock
    if 'bonds' in buildingBlockTopologies[b]:
        result['bondEdges'].update(
                translate_bonds_to_edges(
                    buildingBlockTopologies[b]['bonds'],
                    buildingBlockTopologies[b]['atoms'],
                    ffParams['bondtypes']
                )
        )

        if inferAngles:
            inferredAngles = infer_angles(
                    buildingBlockTopologies[b]['atoms'],
                    buildingBlockTopologies[b]['bonds'],
                    ffParams['angletypes'])
            result['angleEdges'].update(
                    translate_angles_to_edges(
                        inferredAngles,
                        buildingBlockTopologies[b]['atoms'],
                        ffParams['bondtypes'],
                        ffParams['angletypes']
                    )
            )
            pass
        pass
        if 'angles' in buildingBlockTopologies[b]:
            result['angleEdges'].update(translate_angles_to_edges(
                buildingBlockTopologies[b]['angles'],
                buildingBlockTopologies[b]['atoms'],
                ffParams['bondtypes'],
                ffParams['angletypes']))
            pass

        if 'impropers' in buildingBlockTopologies[b]:
            result['improperEdges'].update(
                    translate_impropers_to_edges(
                        buildingBlockTopologies[b]['impropers'],
                        result['angleEdges'],
                        result['bondEdges'],
                        buildingBlockTopologies[b]['atoms'],
                        ffParams['dihedraltypes'])
            )
        pass
    return result


def generate(
        ffName,
        alphabets=[],
        inferAngles=True,
        topPath=topologyPath):
    """
    Generate a topology database from force field parameters out of a gromacs force field.

    Topology is read from .rtp files. Parameters are read from .itp files.

    NOTE on the structure of this: in the structure any type of residue has to have a unique name.
    That is why I chose to translate the letter from an alphabet to a unique building block name
    instead of putting in one database for every type of polymer.
    Terminal residues will have to be handled slightly differently.
    """

    result = {}
    result['alphabets'] = {}
    for a in alphabets:
        polymerType = data.polymer_type(a)
        result['alphabets'][polymerType] = dict()
        pass

    ffDir = topPath + ffName + '.ff/'

    residueTopologyFiles = [ffDir + f for f in os.listdir(ffDir) if f.endswith('.rtp')]
    forceFieldParamFiles = [ffDir + "ffbonded.itp"]

    params = {}
    macros = {}
    for f in forceFieldParamFiles:
        p, m = parse_forcefield_params(f)
        params.update(p)
        macros.update(m)
        pass

    buildingBlockTopologies = {}
    for f in residueTopologyFiles:
        buildingBlockTopologies.update(parse_residue_topology(f, macros))
        pass

    for a in alphabets:
        polymerType = data.polymer_type(a)
        for letter in a.letters:
            buildingBlock = letter.upper()
            key = buildingBlock

            if polymerType =='AA':
                if a.size == 1:
                    buildingBlock = IUPACData.protein_letters_1to3_extended[letter]
                    buildingBlock = buildingBlock.upper()
                    key = buildingBlock
                    if buildingBlock not in buildingBlockTopologies:
                        buildingBlock = buildingBlock[:2] + 'P'  # positively charged histidine
                        pass
                    if buildingBlock not in buildingBlockTopologies:
                        raise KeyError("Could not find " + letter + " in the residue database.")
                    pass
            else:
                if buildingBlock not in buildingBlockTopologies:
                    if polymerType == 'DNA':
                        buildingBlock = 'D' + buildingBlock
                        key = buildingBlock
                        pass
                    elif polymerType == 'RNA':
                        key = buildingBlock
                        buildingBlock = 'R' + buildingBlock
                        pass
                    else:
                        raise
                    pass
                pass
            result['alphabets'][polymerType][letter] = key
            result[key] = translate(buildingBlock, buildingBlockTopologies, params)
            pass
        pass

    return result
