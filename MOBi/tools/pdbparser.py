#!/usr/bin/env python3
#####################################
#
# Filename : pdbparser.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Thu 11 May 2017 16:35:51 CEST
#
# Last Modified : Sun 31 Dec 2017 04:13:43 PM CET
#
#####################################

import os

import numpy as np

import Bio
from Bio import PDB
# from Bio import SeqUtils
from Bio import Seq
from Bio import SeqIO
# from Bio import Alphabet
# from Bio.Data.IUPACData import protein_letters_1to3

from MOBi import data
from MOBi.tools import math


# NOTE added chemicalDB to check for differing atom naming convention
# TODO add missing residues from a sequence and use distances from force field at the read pdb step?
def read_PDB(PDBCode, fileName, chemicalDB={}, assign_serial_numbers=True):
    # TODO maybe handle hetero stuff if they are present in the force field

    # TODO Glycin Hs are not read correctly!!! (others probably too)

    structure = None

    root, ext = os.path.splitext(fileName)

    def _read_structure_file(structureName, structureFile, mode='.cif'):
        nonlocal assign_serial_numbers
        parser = None
        if mode == '.cif':
            assign_serial_numbers = True
            parser = PDB.MMCIFParser()
        else:
            parser = PDB.PDBParser()
            pass
        return parser.get_structure(structureName, structureFile)

    if ext == '.gz':
        import gzip
        with gzip.open(fileName, 'rt') as f:
            structure = _read_structure_file(PDBCode, f, os.path.splitext(root)[1])
            pass
        pass
    elif ext == '.bz2':
        import bzip
        with bzip.open(fileName, 'rt') as f:
            structure = _read_structure_file(PDBCode, f, os.path.splitext(root)[1])
    else:
        structure = _read_structure_file(PDBCode, fileName, ext)
        pass

    # NOTE this is a bit clunky but it is apparently a bad idea to remove bits of an iterable while iterating over it
    flaggedForRemoval = []
    for model in structure:
        for chain in model:
            for residue in chain:
                if residue.id[0] != ' ':
                    # NOTE this removes e.g. hetero residues
                    flaggedForRemoval.append((chain.id, residue.id))
                    assign_serial_numbers = True
                else:
                    if residue.get_resname() not in chemicalDB:
                        flaggedForRemoval.append((chain.id, residue.id))
                        assign_serial_numbers = True
                        pass
                    # TODO remove this and handle termini properly
                    if residue.has_id('OXT'):
                        residue.detach_child('OXT')
                        pass
                    if chemicalDB:
                        if residue.get_resname() in chemicalDB:
                            for atom in residue:
                                if atom.get_name() not in chemicalDB[residue.get_resname()]['vertices']:
                                    if not residue.has_id(atom.get_name()[:-1]) and atom.get_name()[:-1] in chemicalDB[residue.get_resname()]['vertices']:
                                        # NOTE sometimes atoms in pdbfiles are called e.g. CD1 instead of CD
                                        # print(str(residue), atom.get_name(), " not in ", chemicalDB[residue.get_resname()]['vertices'], " of " + residue.get_resname())
                                        oldID = atom.get_id()
                                        newID = oldID[:-1]
                                        atom.name = newID
                                        atom.id = newID
                                        atom.parent.child_dict[newID] = atom.parent.child_dict[oldID]  # NOTE has_id checks in the parents child_dict. which still has the old id as key. since atoms are not entities for some reason, can't use builtin facilities to manage this
                                        del atom.parent.child_dict[oldID]  # TODO maybe this should be done outside the loop
                                        # print("renaming to ", atom.get_name())
                                        # print(chemicalDB[residue.get_resname()]['bondEdges'])

                                    else:
                                        print(str(residue), atom.get_name(), " not in ", chemicalDB[residue.get_resname()]['vertices'], " of " + residue.get_resname())
                                        # TODO handle this properly
                                        assert False
                                        pass
                                    pass
                                pass
                            pass
                        else:
                            print("unkown residue: ", str(residue), "in chain ", str(chain), " flagged for removal")
                            flaggedForRemoval.append((chain.get_id(), residue.get_id()))
                            assign_serial_numbers = True
                            pass
                    pass
                pass
            pass
        for cr in flaggedForRemoval:
            residue = model[cr[0]][cr[1]]
            model[cr[0]].detach_child(cr[1])
            pass
        flaggedForRemoval = []
        for chain in model:
            if len(chain) == 0:
                flaggedForRemoval.append(chain.id)
                pass
            pass
        for chain in flaggedForRemoval:
            model.detach_child(chain)
        pass

    # NOTE this is necessary after removing some atoms or when loading a .cif, since the structure builder sets None for all serial numbers in that case
    # TODO should happen when serial numbers start at 1 too!!
    if assign_serial_numbers:
        i = 0
        for model in structure:
            for chain in model:
                for residue in chain:
                    for atom in residue:
                        atom.set_serial_number(i)
                        i += 1
                        pass
                    pass
                pass
            pass
        pass

    # aas = {}
    # print(structure)
    # for chain in model:
    #     print(chain)
    #     for residue in chain:
    #         if residue.get_resname() in aas:
    #             aas[residue.get_resname()] += 1
    #         else:
    #             aas[residue.get_resname()] = 1
    #         print('\t' + str(residue.get_id()) + " " + str(residue.get_resname()))
    #         for atom in residue:
    #             print('\t\t' + str(atom.get_serial_number()) + " " + str(atom.get_id()) + " " + str(atom.get_name()))
    #             pass
    #         pass
    #     pass
    # print(aas)
    print(len(list(structure.get_atoms())), ' atoms parsed')

    return structure


def read_sequences_from_PDB(fileName):
    # TODO choose correct alphabet
    sequences = {}
    with open(fileName, 'rU') as f:
        for record in SeqIO.parse(f, 'pdb-seqres'):
            sequence = Seq(record.seq)
            # print(">" + record.id + "\n" + record.seq)
            sequences[record.id] = sequence
            pass
        pass
    return sequences


# # adds Hs and missing atoms.
# # TODO handle missing residues
# def prepare_structure(structure):
#     # TODO
#     result = PDB.structure()
#
#     # NOTE only uses first model
#     model = structure[0]
#     for chain in model:
#         for residue in chain:
#             pass
#         pass
#
#     return result


def get_edge(atom1, atom2):
    edge = tuple(sorted([atom1.get_serial_number(), atom2.get_serial_number()]))
    distance = np.sqrt(np.dot(atom1.get_coord() - atom2.get_coord(), atom1.get_coord() - atom2.get_coord()))
    return edge, distance


# TODO this is only for backbone dihedrals, rename
def get_dihedral_edge(model, chainID, resID, atomID):
    if chainID in model:
        if resID in model[chainID] and resID - 1 in model[chainID]:
            if model[chainID][resID].has_id(atomID) and model[chainID][resID - 1].has_id(atomID):
                atom1 = model[chainID][resID][atomID]
                atom2 = model[chainID][resID - 1][atomID]
                return get_edge(atom1, atom2)
                pass
            pass
        pass
    return None, None


# TODO parameter names
def get_primary_edges(chain, chemicalDB, alphabet=data.PDBReducedProtein, useStructureDistances=True):
    # TODO handle a structure containing chains with RNA and AA
    # NOTE this takes only atoms, that are in the structure. to use all atoms: read in sequence, generate vertices and edges with a force field
    # TODO edge types as parameter
    # TODO handle termini

    edges = []
    distances = []
    weights = []

    # TODO how to handle missing residues (numbering scheme)
    for residue in chain:
        resn = residue.get_resname()
        for edgeType in ['bondEdges', 'angleEdges', 'improperEdges']:
            for edge in chemicalDB[resn][edgeType]:
                residueIDs = [residue.get_id()[1], residue.get_id()[1]]
                atomNames = []
                edgeAtomSerialNumbers = []
                atomCoordinates = []
                for i in range(len(edge)):
                    if '+' in edge[i]:
                        residueIDs[i] += 1
                    elif '-' in edge[i]:
                        residueIDs[i] -= 1
                        pass
                    atomNames.append(edge[i].strip('+-'))
                    pass
                if residueIDs[0] in chain and residueIDs[1] in chain and chain[residueIDs[0]].has_id(atomNames[0]) and chain[residueIDs[1]].has_id(atomNames[1]):
                    edgeAtomSerialNumbers = [chain[residueIDs[0]][atomNames[0]].get_serial_number(), chain[residueIDs[1]][atomNames[1]].get_serial_number()]
                    atomCoordinates = [chain[residueIDs[0]][atomNames[0]].get_coord(), chain[residueIDs[1]][atomNames[1]].get_coord()]
                    pass
                else:
                    # TODO warn?
                    # TODO for this, proper handling of terminals is necessary
                    pass
                if len(edgeAtomSerialNumbers) == 2:
                    edges.append(tuple(sorted(edgeAtomSerialNumbers)))
                    if useStructureDistances:
                        distances.append(np.sqrt(np.dot(atomCoordinates[0] - atomCoordinates[1], atomCoordinates[0] - atomCoordinates[1])))
                    else:
                        distances.append(chemicalDB[resn][edgeType][edge])
                        pass

                    weights.append(1.)
                    pass
                pass
            pass
        pass

    # TODO warn if atoms that are in a database building block are missing in the structure
    return edges, distances, weights


def get_secondary_sequence_protein(fileName, topologyDB):
    # TODO the model seems unnecessary
    result = {}
    # model = Bio.PDB.PDBParser().get_structure('sec', fileName)[0]
    model = read_PDB('sec', fileName, topologyDB)[0]
    dssp = Bio.PDB.DSSP(model, fileName)

    for chain in model:
        result[chain.get_id()] = []
        pass

    for key in dssp.keys():
        chainID = key[0]
        resID = key[1][1]
        letter = dssp[key][2]

        if letter == 'H':
            result[chainID].append((letter, resID))
        elif letter == 'E':
            result[chainID].append((letter, resID))
        else:
            result[chainID].append(('-', resID))
            pass
        pass

    #################
    for chainID in result:
        print([x[0] for x in result[chainID]])
        print([x[1] for x in result[chainID]])
        pass
    #################

    return result


def get_secondary_edges_protein(model, secondaryStructureSequence, topologyDB, useStructureDistances=True):
    # TODO this works only for protein right now
    # TODO load edges to be used from MOBi.data
    # TODO add distance from database

    edges = []
    distances = []
    weights = []

    for chain in model:
        chainID = chain.get_id()
        if chainID in secondaryStructureSequence:
            # TODO check if first or last residue
            for i, r in enumerate(secondaryStructureSequence[chainID]):
                print('###')
                print(r)
                resID = r[1]
                letter = r[0]
                resn = chain[resID].get_resname()
                # prevresn = None
                # nextresn = None
                # if resID - 1 in chain:
                #     prevresn = chain[resID - 1].get_resname()
                #     pass
                # if resID + 1 in chain:
                #     prevresn = chain[resID - 1].get_resname()
                #     pass
                edge = None
                distance = None
                # omega
                if resn != 'PRO':  # TODO alphabet
                    edge, distance = get_dihedral_edge(model, chainID, resID, 'CA')
                    if edge:
                        print("omega:")
                        # ############## should be 3.8 angstrom for trans
                        A = chain[resID - 1]['CA'].get_vector()
                        B = chain[resID - 1]['C'].get_vector()
                        C = chain[resID]['N'].get_vector()
                        D = chain[resID]['CA'].get_vector()
                        dihedral = PDB.calc_dihedral(A, B, C, D)
                        # ##############
                        print("measured dihedral: ", dihedral*180./ np.pi)
                        print("measured distance: ", distance)
                        if not useStructureDistances:
                            prevresn = chain[resID - 1].get_resname()
                            dAB = topologyDB[prevresn]['bondEdges'][('C', 'CA')]
                            dBC = topologyDB[resn]['bondEdges'][('-C', 'N')]
                            dCD = topologyDB[resn]['bondEdges'][('CA', 'N')]
                            dAC = topologyDB[prevresn]['angleEdges'][('+N', 'CA')]
                            dBD = topologyDB[resn]['angleEdges'][('-C', 'CA')]
                            distance = math.dist_from_dihedral(180., dAB, dBC, dCD, dAC, dBD)
                            print("guessed angle: ", 180.)
                            print("guessed distance: ", distance)
                            pass
                        edges.append(edge)
                        distances.append(distance)
                        weights.append(1.)
                        print('omega', edge, resID)
                        pass
                    pass
                edge = None
                distance = None
                if letter == 'H':
                    # phi
                    edge, distance = get_dihedral_edge(model, chainID, resID, 'C')
                    if edge:
                        print("helix")
                        print("phi:")
                        # ##############
                        A = chain[resID - 1]['C'].get_vector()
                        B = chain[resID]['N'].get_vector()
                        C = chain[resID]['CA'].get_vector()
                        D = chain[resID]['C'].get_vector()
                        dihedral = PDB.calc_dihedral(A, B, C, D)
                        # ##############
                        print("measured distance: ", distance)
                        print("measured dihedral: ", dihedral * 180. / np.pi)
                        if not useStructureDistances:
                            dAB = topologyDB[resn]['bondEdges'][('-C', 'N')]
                            dBC = topologyDB[resn]['bondEdges'][('CA', 'N')]
                            dCD = topologyDB[resn]['bondEdges'][('C', 'CA')]
                            dAC = topologyDB[resn]['angleEdges'][('-C', 'CA')]
                            dBD = topologyDB[resn]['angleEdges'][('C', 'N')]
                            distance = math.dist_from_dihedral(data.alpha_phi, dAB, dBC, dCD, dAC, dBD)
                            print("guessed distance: ", distance)
                            print("guessed dihedral: ", data.alpha_phi)
                            pass
                        edges.append(edge)
                        distances.append(distance)
                        weights.append(1.)
                        print('phi', edge, resID)
                        pass
                    # psi
                    edge, distance = get_dihedral_edge(model, chainID, resID + 1, 'N')
                    if edge:
                        print("helix")
                        print("psi:")
                        # ##############
                        A = chain[resID]['N'].get_vector()
                        B = chain[resID]['CA'].get_vector()
                        C = chain[resID]['C'].get_vector()
                        D = chain[resID + 1]['N'].get_vector()
                        dihedral = PDB.calc_dihedral(A, B, C, D)
                        # ##############
                        print("measured distance: ", distance)
                        print("measured dihedral: ", dihedral * 180. / np.pi)
                        if not useStructureDistances:
                            nextresn = chain[resID + 1].get_resname()
                            dAB = topologyDB[resn]['bondEdges'][('CA', 'N')]
                            dBC = topologyDB[resn]['bondEdges'][('C', 'CA')]
                            dCD = topologyDB[nextresn]['bondEdges'][('-C', 'N')]
                            dAC = topologyDB[resn]['angleEdges'][('-C', 'CA')]
                            dBD = topologyDB[resn]['angleEdges'][('+N', 'CA')]
                            distance = math.dist_from_dihedral(data.alpha_psi, dAB, dBC, dCD, dAC, dBD)
                            print("guessed distance: ", distance)
                            print("guessed dihedral: ", data.alpha_psi)
                            pass
                        edges.append(edge)
                        distances.append(distance)
                        weights.append(1.)
                        print('psi', edge, resID)
                        pass
                    # h bond
                    # TODO include H
                    # TODO test this
                    if chain.has_id(resID - 4):
                        print('helix')
                        print('H bond')
                        if secondaryStructureSequence[chainID][i - 4][0] == 'H':
                            if secondaryStructureSequence[chainID][i - 4][1] == resID - 4:
                                atom1 = model[chainID][resID]['N']
                                atom2 = model[chainID][resID - 4]['C']
                                edge, distance = get_edge(atom1, atom2)
                                print("measured distance: ", distance)
                                if edge:
                                    # TODO maybe add edge between acceptor and hydrogen
                                    if not useStructureDistances:
                                        dAB = topologyDB[resn]['bondEdges'][('H', 'N')]
                                        dBC = data.hbond
                                        angle = data.hangle
                                        distance = math.dist_from_angle(angle, dAB, dBC)
                                        print("guessed distance: ", distance)
                                        pass
                                    edges.append(edge)
                                    distances.append(distance)
                                    weights.append(1.)
                                    # print('H', edge, resID)
                                    # print(distance)
                                    pass
                                pass
                            pass
                        pass
                elif letter == 'E':
                    # phi
                    edge, distance = get_dihedral_edge(model, chainID, resID, 'C')
                    if edge:
                        if not useStructureDistances:
                            dAB = topologyDB[resn]['bondEdges'][('-C', 'N')]
                            dBC = topologyDB[resn]['bondEdges'][('CA', 'N')]
                            dCD = topologyDB[resn]['bondEdges'][('C', 'CA')]
                            dAC = topologyDB[resn]['angleEdges'][('-C', 'CA')]
                            dBD = topologyDB[resn]['angleEdges'][('C', 'N')]
                            distance = math.dist_from_dihedral(data.beta_phi, dAB, dBC, dCD, dAC, dBD)
                            print("guessed distance: ", distance)
                            pass
                        edges.append(edge)
                        distances.append(distance)
                        weights.append(1.)
                        pass
                    # psi
                    edge, distance = get_dihedral_edge(model, chainID, resID + 1, 'N')
                    if edge:
                        if not useStructureDistances:
                            nextresn = chain[resID + 1].get_resname()
                            dAB = topologyDB[resn]['bondEdges'][('CA', 'N')]
                            dBC = topologyDB[resn]['bondEdges'][('C', 'CA')]
                            dCD = topologyDB[nextresn]['bondEdges'][('-C', 'N')]
                            dAC = topologyDB[resn]['angleEdges'][('-C', 'CA')]
                            dBD = topologyDB[resn]['angleEdges'][('+N', 'CA')]
                            distance = math.dist_from_dihedral(data.beta_psi, dAB, dBC, dCD, dAC, dBD)
                            print("guessed distance: ", distance)
                        edges.append(edge)
                        distances.append(distance)
                        weights.append(1.)
                        pass
                elif letter == '-':
                    pass
                else:
                    raise
                pass
            pass
        pass

    # dssp = PDB.DSSP(model, fileName)

    # for key in dssp.keys():
    #     chainID = key[0]
    #     resID = key[1][1]

    #     # omega
    #     edge, distance = get_dihedral_edge(model, chainID, resID, 'CA')
    #     if edge:
    #         edges.append(edge)
    #         if useStructureDistances:
    #             distances.append(distance)
    #         else:
    #             raise NotImplementedError()
    #         # weights.append(5.)
    #         weights.append(1.)
    #         pass
    #     # if chainID in model:
    #     #     if resID - 1 in model[chainID] and resID in model[chainID]:
    #     #         if model[chainID][resID - 1].has_id('CA') and model[chainID][resID].has_id('CA'):
    #     #             atom1 = model[chainID][resID - 1]['CA']
    #     #             atom2 = model[chainID][resID]['CA']
    #     #             edge, distance = get_edge(atom1, atom2)
    #     #             edges.append(edge)
    #     #             if useStructureDistances:
    #     #                 distances.append(distance)
    #     #             else:
    #     #                 raise NotImplementedError()
    #     #                 pass
    #     #             weights.append(5.)
    #     #             pass
    #     #         pass
    #     #     pass

    #     if dssp[key][2] == 'H':  # alpha Helix
    #         # phi
    #         edge, distance = get_dihedral_edge(model, chainID, resID - 1, 'C')
    #         if edge:
    #             edges.append(edge)
    #             if useStructureDistances:
    #                 distances.append(distance)
    #             else:
    #                 raise NotImplementedError()
    #             # weights.append(5.1)
    #             weights.append(1.)
    #             pass

    #         # psi
    #         edge, distance = get_dihedral_edge(model, chainID, resID, 'N')
    #         if edge:
    #             edges.append(edge)
    #             if useStructureDistances:
    #                 distances.append(distance)
    #             else:
    #                 raise NotImplementedError()
    #             # weights.append(5.2)
    #             weights.append(1.)
    #             pass

    #         # H bond
    #         if chainID in model:
    #             if resID in model[chainID] and resID + 4 in model[chainID]:
    #                 if model[chainID][resID].has_id('N') and model[chainID][resID + 4].has_id('C'):
    #                     # TODO take H atom into account
    #                     atom1 = model[chainID][resID]['N']
    #                     atom2 = model[chainID][resID]['C']
    #                     edge, distance = get_edge(atom1, atom2)
    #                     if edge:
    #                         edges.append(edge)
    #                         if useStructureDistances:
    #                             distances.append(distance)
    #                         else:
    #                             raise NotImplementedError()
    #                         # weights.append(5.3)
    #                         weights.append(1.)
    #                         pass
    #                 pass
    #             pass
    #         pass
    #     elif dssp[key][2] == 'E':  # beta strand

    #         # phi
    #         edge, distance = get_dihedral_edge(model, chainID, resID - 1, 'C')
    #         if edge:
    #             edges.append(edge)
    #             if useStructureDistances:
    #                 distances.append(distance)
    #             else:
    #                 raise NotImplementedError()
    #             # weights.append(6.1)
    #             weights.append(1.)
    #             pass

    #         # psi
    #         edge, distance = get_dihedral_edge(model, chainID, resID, 'N')
    #         if edge:
    #             edges.append(edge)
    #             if useStructureDistances:
    #                 distances.append(distance)
    #             else:
    #                 raise NotImplementedError()
    #             # weights.append(6.2)
    #             weights.append(1.)
    #             pass
    #         pass
    #     else:
    #         pass
    #     pass

    return edges, distances, weights


# NOTE default cutOff in units of the structure (assumed angstrom)
# TODO only return contacts? use fullID?
def get_tertiary_edges(model, cutOff=5., minSeqDist=5, getContacts=False):
    # NOTE these would include beta sheet hbonds

    tertiaryEdges = []
    distances = []
    weights = []
    contacts = []

    ns = PDB.NeighborSearch(list(model.get_atoms()))
    foundPairs = ns.search_all(cutOff)

    for pair in foundPairs:
        edge = tuple(sorted((pair[0].get_serial_number(), pair[1].get_serial_number())))
        if abs(pair[0].get_parent().get_id()[1] - pair[1].get_parent().get_id()[1]) >= minSeqDist:
            tertiaryEdges.append(edge)
            chainID1 = pair[0].get_parent().get_parent().get_id()
            chainID2 = pair[1].get_parent().get_parent().get_id()
            resID1 = pair[0].get_parent().get_id()[1]
            resID2 = pair[1].get_parent().get_id()[1]
            atom1 = pair[0].get_id()
            atom2 = pair[1].get_id()
            contacts.append(((chainID1, resID1, atom1), (chainID2, resID2, atom2)))
            distances.append(np.sqrt(np.dot(pair[0].get_coord() - pair[1].get_coord(), pair[0].get_coord() - pair[1].get_coord())))
            weights.append(1.)
            pass
        pass

    if getContacts:
        return tertiaryEdges, distances, weights, contacts
        pass
    return tertiaryEdges, distances, weights


# TODO better naming
def translate_to_edges(contacts, model):
    edges = []
    for contact in contacts:
        chainID1 = contact[0][0]
        chainID2 = contact[1][0]
        resID1 = contact[0][1]
        resID2 = contact[1][1]
        atom1 = contact[0][2]
        atom2 = contact[1][2]

        # print(model[chainID1][resID1].get_resname())
        # print(list(model[chainID1][resID1].child_dict.keys()))
        # for atom in model[chainID1][resID1]:
        #     print(atom.get_id())
        #     pass
        # print(model[chainID2][resID2].get_resname())
        # print(list(model[chainID2][resID2].child_dict.keys()))
        # for atom in model[chainID2][resID2]:
        #     print(atom.get_id())
        #     pass

        index1 = model[chainID1][resID1][atom1].get_serial_number()
        index2 = model[chainID2][resID2][atom2].get_serial_number()
        edges.append(tuple(sorted((index1, index2))))
        pass

    return edges


# TODO better naming
# def translate_secondary_structure(model, sequences, topologyDB):
#     edges = []
#     distances = []
#     weights = []
#
#     return edges, distances, weights


# TODO add parameters
# TODO work on model instead of structure?
def generate_graph(structure, fileName, topologyDB, cutOff=3., minSeqDist=5):

    atoms = list(structure[0].get_atoms())

    edges = []
    distances = []
    weights = []

    for chain in structure[0]:
        chainEdges, chainDistances, chainWeights = get_primary_edges(chain, topologyDB)
        edges += chainEdges
        distances += chainDistances
        weights += chainWeights
        pass

    print(str(len(atoms)) + " vertices found")
    print(str(len(edges)) + " primary edges found")
    # print(sorted(edges, key=lambda tup: tup[0]))

    # secondaryEdges, secondaryDistances, secondaryWeights = get_secondary_edges(structure[0], fileName)
    # print(str(len(secondaryEdges)) + " secondary edges found")

    # edges += secondaryEdges
    # distances += secondaryDistances
    # weights += secondaryWeights

    tertiaryEdges, tertiaryDistances, tertiaryWeights = get_tertiary_edges(structure[0], edges, cutOff, minSeqDist)
    print(str(len(tertiaryEdges)) + " tertiary edges found")

    # print(sorted(tertiaryEdges, key=lambda tup: tup[0]))

    edges += tertiaryEdges
    distances += tertiaryDistances
    weights += tertiaryWeights

    return atoms, edges, distances, weights


# TODO put this somewhere else
# TODO list of resIDs instead of simple offset
def build_structure(id, sequences, topologyDB, offsets=[]):
    structure = Bio.PDB.Structure.Structure(id)
    model = Bio.PDB.Model.Model(0, None)
    structure.add(model)

    if not offsets:
        offsets = [0] * len(sequences)
        pass

    chainID = 'A'
    atomCounter = 0
    for sequence, offset in zip(sequences, offsets):
        chain = PDB.Chain.Chain(chainID)
        chainID = chr(ord(chainID) + 1)
        model.add(chain)
        for i, res in enumerate(sequence):  # TODO put in proper residue ids
            res_ID = (" ", i + offset, " ")
            # NOTE assume, sequence has the right alphabet, convert before
            # resName = protein_letters_1to3[res].upper()
            resName = res
            segID = "   "
            residue = Bio.PDB.Residue.Residue(res_ID, resName, segID)
            chain.add(residue)

            vertices = set(topologyDB[resName]['vertices'])
            backbone = ['N', 'CA', 'C', 'O']
            # NOTE first four atoms are backbone, rest is alphabetically ordered
            for vertex in backbone:
                atomName = vertex
                coord = np.array((0., 0., 0.), "f")
                bfactor = 0.
                occupancy = 1,
                altloc = " "
                fullName = vertex
                serialNumber = atomCounter
                element = atomName[0]  # TODO valid for C, O, N, S, H, P? in proteins? definitely not universal!!!
                # element = None  # TODO fix this, giving None prints tons of warnings

                atom = Bio.PDB.Atom.Atom(atomName, coord, bfactor, occupancy, altloc, fullName, serialNumber, element)
                residue.add(atom)
                atomCounter += 1
                pass
            vertices -= set(backbone)
            # for vertex in topologyDB[resName]['vertices']:
            for vertex in sorted(vertices):
                atomName = vertex
                coord = np.array((0., 0., 0.), "f")
                bfactor = 0.
                occupancy = 1,
                altloc = " "
                fullName = vertex
                serialNumber = atomCounter
                element = atomName[0]  # TODO valid for C, O, N, S, H, P? in proteins? definitely not universal!!!
                # element = None  # TODO fix this, giving None prints tons of warnings

                atom = Bio.PDB.Atom.Atom(atomName, coord, bfactor, occupancy, altloc, fullName, serialNumber, element)
                residue.add(atom)
                atomCounter += 1
                pass
            pass
        pass
    print(len(list(structure.get_atoms())), ' atoms generated')
    return structure
