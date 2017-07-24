#!/usr/bin/env python3
#####################################
#
# Filename : ffparsergmx.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Thu 11 May 2017 10:54:56 CEST
#
# Last Modified : Mon 26 Jun 2017 16:52:17 CEST
#
#####################################

import os

from . import math

from .. import config

topologyPath = config.gromacs_topology_path


def read_atoms(line):
    parts = line.split()
    name = parts[0]
    type = parts[1]
    return {name: type}


def read_bonds(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    b0 = None
    if len(parts) > 2:
        b0 = float(parts[2])
        pass

    # key = (i, j)
    # value = b0
    # return {key: value}
    return {(i, j): b0, (j, i): b0}


def read_angles(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    th0 = None
    if len(parts) > 3:
        th0 = float(parts[3])
        pass

    # key = (i, j, k)
    # value = th0
    # return {key: value}
    return {(i, j, k): th0, (k, j, i): th0}


# # TODO general dihedrals logic
# def read_dihedrals(line):  # NOTE proper dihedrals are not useful at this point
#     # TODO only add, if single minimum
#     # here function type and defaults matter
#     parts = line.split()
#     i = parts[0]
#     j = parts[1]
#     k = parts[2]
#     l = parts[3]
#     if len(parts) > 4:
#         phi0 = float(parts[4])
#     else:
#         phi0 = None
#         pass
#
#     # NOTE order matters
#     key = (i, j, k, l)
#     value = phi0
#     return {key: value}


def read_impropers(line):
    # TODO only add, if single minimum
    parts = line.split()
    i = parts[0]
    j = parts[1]
    k = parts[2]
    l = parts[3]
    if len(parts) > 4:
        q0 = float(parts[4])
    else:
        q0 = None
        pass

    # NOTE order matters
    key = (i, j, k, l)
    value = q0
    return {key: value}


# NOTE amber type ffs have n- and c- terminus residues, that are used here
# NOTE if buildingBlocks is empty it considers all the ones it can find
# NOTE amber type ffs have different AAs for an AA with different protonations (HIS)
# TODO translate names
def parse_residue_topology(
        residueTopologyFile,
        macros={},
        buildingBlocks=[],
        ignoredDirectives={'bondedtypes', 'exclusions', 'dihedrals'},
        knownDirectives={'atoms': read_atoms, 'bonds': read_bonds, 'angles': read_angles, 'impropers': read_impropers}):

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
            if directive in buildingBlocks or (len(buildingBlocks) == 0 and directive not in knownDirectives):
                buildingBlock = directive
                context = None
                result[buildingBlock] = {}
            elif directive not in buildingBlocks:
                if directive in knownDirectives:
                    if directive not in ignoredDirectives:
                        if buildingBlock:
                            context = directive
                            result[buildingBlock][context] = {}
                            pass
                    else:
                        context = None
                        pass
                else:
                    buildingBlock = None
                    pass
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

    return result


# TODO for consistency: put in both directions
def read_bondtypes(line):
    parts = line.split()
    i = parts[0]
    j = parts[1]
    # func = parts[2]
    b0 = parts[3]

    # TODO maybe don't use frozensets, just add the inverted tuple to the dict as well?
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
    th0 = parts[4]

    return {(i, j, k): float(th0), (k, j, i): float(th0)}


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
        value = float(ph0)
        return {key: value}
    else:
        return {}
    pass


# TODO directives as parameter
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
    # TODO n-termini special casing
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


# def infer_impropers():
#     # TODO from rings? probably not useful. peptide bond / backbone dihedral omega?
#     result = {}
#     return result


def translate_atoms_to_vertices(atoms):
    result = set(atoms.keys())
    return result


# TODO these go somewhere else together with all the add_edge_to_graph tools?
def translate_bonds_to_edges(bonds, atomTypes, bondTypes):
    # TODO n-termini special casing
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

        result[vertices] = math.dist_from_angle(angle, d_ij, d_jk)
        pass
    return result


# TODO switch to bondTypes and angleTypes?
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

        result[vertices] = math.dist_from_improper(improper, d_ij, d_jk, d_kl, d_ik, d_jl)

        pass
    return result


def generate_chemical_primary_edge_database(
        ffname,
        buildingBlocks=[],
        inferAngles=True,
        inferImpropers=False,
        topPath=topologyPath):

    residueTopologieFiles = []
    forceFieldParamFiles = []

    if os.path.isdir(topPath):
        if os.path.isdir(topPath + ffname + '.ff'):
            residueTopologieFiles = [f for f in os.listdir(topPath + ffname + '.ff/') if f.endswith('.rtp')]
            forceFieldParamFiles.append(topPath + ffname + '.ff/ffbonded.itp')
            pass
        pass

    ffParams = {}
    macros = {}
    for forceFieldParamFile in forceFieldParamFiles:
        ffP, m = parse_forcefield_params(forceFieldParamFile)
        ffParams.update(ffP)
        macros.update(m)
        pass

    buildingBlockTopologies = {}
    assert(len(residueTopologieFiles) > 0)
    for residueTopologyFile in residueTopologieFiles:
        buildingBlockTopologies.update(parse_residue_topology(topPath + ffname + '.ff/' + residueTopologyFile, macros, buildingBlocks))
        pass

    for buildingBlock in buildingBlocks:
        if buildingBlock not in buildingBlockTopologies:
            # TODO more helpful error messages / warnings
            raise
        pass

    result = {}
    for buildingBlock in buildingBlockTopologies:
        print(buildingBlock)

        result[buildingBlock] = {}

        result[buildingBlock]['vertices'] = set()
        result[buildingBlock]['bondEdges'] = {}
        result[buildingBlock]['angleEdges'] = {}
        result[buildingBlock]['improperEdges'] = {}

        result[buildingBlock]['vertices'].update(
                translate_atoms_to_vertices(
                    buildingBlockTopologies[buildingBlock]['atoms']
                )
        )

        if 'bonds' in buildingBlockTopologies[buildingBlock]:
            result[buildingBlock]['bondEdges'].update(
                    translate_bonds_to_edges(
                        buildingBlockTopologies[buildingBlock]['bonds'],
                        buildingBlockTopologies[buildingBlock]['atoms'],
                        ffParams['bondtypes']
                    )
            )

            if inferAngles:
                inferredAngles = infer_angles(
                        buildingBlockTopologies[buildingBlock]['atoms'],
                        buildingBlockTopologies[buildingBlock]['bonds'],
                        ffParams['angletypes'])
                result[buildingBlock]['angleEdges'].update(
                        translate_angles_to_edges(
                            inferredAngles,
                            buildingBlockTopologies[buildingBlock]['atoms'],
                            ffParams['bondtypes'],
                            ffParams['angletypes']
                        )
                )
                pass
            pass
            if 'angles' in buildingBlockTopologies[buildingBlock]:
                result[buildingBlock]['angleEdges'].update(translate_angles_to_edges(
                    buildingBlockTopologies[buildingBlock]['angles'],
                    buildingBlockTopologies[buildingBlock]['atoms'],
                    ffParams['bondtypes'],
                    ffParams['angletypes']))
                pass

            # if inferImpropers:
            #     inferredImpropers = infer_impropers()  # TODO
            #     result[buildingBlock]['improperEdges'].update(translate_impropers_to_edges(
            #         inferredImpropers,
            #         result[buildingBlock]['angleEdges'],
            #         result[buildingBlock]['bondEdges'],
            #         buildingBlockTopologies[buildingBlock]['atoms'],
            #         ffParams['dihedraltypes']))
            #     pass

            if 'impropers' in buildingBlockTopologies[buildingBlock]:
                result[buildingBlock]['improperEdges'].update(
                        translate_impropers_to_edges(
                            buildingBlockTopologies[buildingBlock]['impropers'],
                            result[buildingBlock]['angleEdges'],
                            result[buildingBlock]['bondEdges'],
                            buildingBlockTopologies[buildingBlock]['atoms'],
                            ffParams['dihedraltypes'])
                )
            pass
        pass
    return result