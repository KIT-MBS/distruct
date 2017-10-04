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
# Last Modified : Wed 04 Oct 2017 09:49:05 AM CEST
#
#####################################

import os

import numpy as np

from Bio import PDB
# from Bio import SeqUtils
from Bio import Seq
from Bio import SeqIO
# from Bio import Alphabet

from MOBi import data


def read_PDB(PDBCode, fileName, assign_serial_numbers=True):
    # TODO maybe handle hetero stuff if they are present in the force field

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
                    flaggedForRemoval.append((chain.id, residue.id))
                    assign_serial_numbers = True
                    pass
                # TODO remove this and handle termini properly
                if residue.has_id('OXT'):
                    residue.detach_child('OXT')
                pass
            pass
        for residue in flaggedForRemoval:
            model[residue[0]].detach_child(residue[1])
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

    # print structure
    # for model in structure:
    #     for chain in model:
    #         for residue in chain:
    #             print(str(residue.get_id()) + " " + str(residue.get_resname()))
    #             for atom in residue:
    #                 print(str(atom.get_serial_number()) + " " + str(atom.get_id()))
    #                 pass
    #             pass
    #         pass
    #     pass

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


def get_dihedral_edge(model, chainID, resID, atomID):
    if chainID in model:
        if resID in model[chainID] and resID + 1 in model[chainID]:
            if model[chainID][resID].has_id(atomID) and model[chainID][resID + 1].has_id(atomID):
                atom1 = model[chainID][resID][atomID]
                atom2 = model[chainID][resID + 1][atomID]
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
    # TODO add missing residues from a sequence and use distances from force field at the read pdb step?

    # TODO use get_edge(atom1, atom2) function
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
                    # edge[i].strip('+-')
                    atomNames.append(edge[i].strip('+-'))
                    pass
                # TODO this should have worked as well?
                # if resID in chain and chain[resID].has_id(atomName):
                #     atomIDs.append(chain[resID][atomName].get_serial_number())
                #     atomCoordinates.append(chain[resID][atomName].get_coord())
                if residueIDs[0] in chain and residueIDs[1] in chain and chain[residueIDs[0]].has_id(atomNames[0]) and chain[residueIDs[1]].has_id(atomNames[1]):
                    edgeAtomSerialNumbers = [chain[residueIDs[0]][atomNames[0]].get_serial_number(), chain[residueIDs[1]][atomNames[1]].get_serial_number()]
                    atomCoordinates = [chain[residueIDs[0]][atomNames[0]].get_coord(), chain[residueIDs[1]][atomNames[1]].get_coord()]
                    pass
                else:
                    # TODO warn?
                    pass
                if len(edgeAtomSerialNumbers) == 2:
                    edges.append(tuple(sorted(edgeAtomSerialNumbers)))
                    if useStructureDistances:
                        distances.append(np.sqrt(np.dot(atomCoordinates[0] - atomCoordinates[1], atomCoordinates[0] - atomCoordinates[1])))
                    else:
                        distances.append(chemicalDB[resn][edgeType][edge])
                        pass
                    # TODO weights from database
                    weights.append(10.)
                    pass
                pass
            pass
        pass

    # TODO warn if atoms that are in a database building block are missing in the structure
    return edges, distances, weights


# TODO this should be able to process RNA and proteins!
def get_secondary_edges(model, fileName, useStructureDistances=True):
    # TODO this works only for protein right now
    # TODO load edges to be used from MOBi.data
    # TODO add distance from database

    edges = []
    distances = []
    weights = []

    dssp = PDB.DSSP(model, fileName)

    for key in dssp.keys():
        chainID = key[0]
        resID = key[1][1]

        # omega
        edge, distance = get_dihedral_edge(model, chainID, resID, 'CA')
        if edge:
            edges.append(edge)
            if useStructureDistances:
                distances.append(distance)
            else:
                raise NotImplementedError()
            weights.append(5.)
            pass
        # if chainID in model:
        #     if resID - 1 in model[chainID] and resID in model[chainID]:
        #         if model[chainID][resID - 1].has_id('CA') and model[chainID][resID].has_id('CA'):
        #             atom1 = model[chainID][resID - 1]['CA']
        #             atom2 = model[chainID][resID]['CA']
        #             edge, distance = get_edge(atom1, atom2)
        #             edges.append(edge)
        #             if useStructureDistances:
        #                 distances.append(distance)
        #             else:
        #                 raise NotImplementedError()
        #                 pass
        #             weights.append(5.)
        #             pass
        #         pass
        #     pass

        if dssp[key][2] == 'H':  # alpha Helix
            # phi
            edge, distance = get_dihedral_edge(model, chainID, resID - 1, 'C')
            if edge:
                edges.append(edge)
                if useStructureDistances:
                    distances.append(distance)
                else:
                    raise NotImplementedError()
                weights.append(5.)
                pass

            # psi
            edge, distance = get_dihedral_edge(model, chainID, resID, 'N')
            if edge:
                edges.append(edge)
                if useStructureDistances:
                    distances.append(distance)
                else:
                    raise NotImplementedError()
                weights.append(5.)
                pass

            # H bond
            if chainID in model:
                if resID in model[chainID] and resID + 4 in model[chainID]:
                    if model[chainID][resID].has_id('N') and model[chainID][resID + 4].has_id('C'):
                        # TODO take H atom into account
                        atom1 = model[chainID][resID]['N']
                        atom2 = model[chainID][resID]['C']
                        edge, distance = get_edge(atom1, atom2)
                        if edge:
                            edges.append(edge)
                            if useStructureDistances:
                                distances.append(distance)
                            else:
                                raise NotImplementedError()
                            weights.append(5.)
                            pass
                    pass
                pass
            pass
        elif dssp[key][2] == 'E':  # beta strand

            # phi
            edge, distance = get_dihedral_edge(model, chainID, resID - 1, 'C')
            if edge:
                edges.append(edge)
                if useStructureDistances:
                    distances.append(distance)
                else:
                    raise NotImplementedError()
                weights.append(5.)
                pass

            # psi
            edge, distance = get_dihedral_edge(model, chainID, resID, 'N')
            if edge:
                edges.append(edge)
                if useStructureDistances:
                    distances.append(distance)
                else:
                    raise NotImplementedError()
                weights.append(5.)
                pass
            pass
        else:
            pass
        pass

    return edges, distances, weights


# NOTE default cutOff in units of the structure (assumed angstrom)
# TODO add sequence cutoff
def get_tertiary_edges(model, edges, cutOff=3., minSeqDist=5):
    # NOTE these would include beta sheet hbonds
    # TODO distinguish inter and intra chain contacts?
    # TODO add noise?

    tertiaryEdges = []
    distances = []
    weights = []

    ns = PDB.NeighborSearch(list(model.get_atoms()))
    foundPairs = ns.search_all(cutOff)

    for pair in foundPairs:
        edge = tuple(sorted((pair[0].get_serial_number(), pair[1].get_serial_number())))
        # TODO learn more about neighoring residue interactions
        if edge not in edges and abs(pair[0].get_parent().get_id()[1] - pair[1].get_parent().get_id()[1]) >= minSeqDist:
            tertiaryEdges.append(edge)
            distances.append(np.sqrt(np.dot(pair[0].get_coord() - pair[1].get_coord(), pair[0].get_coord() - pair[1].get_coord())))
            weights.append(1.)
            pass
        pass

    # TODO get this stuff in terms of chainID, residueID, atomID
    return tertiaryEdges, distances, weights


# TODO add parameters
# TODO work on model instead of structure?
def generateGraph(structure, fileName, topologyDB, cutOff=3., minSeqDist=5):

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

    secondaryEdges, secondaryDistances, secondaryWeights = get_secondary_edges(structure[0], fileName)

    tertiaryEdges, tertiaryDistances, tertiaryWeights = get_tertiary_edges(structure[0], edges, cutOff, minSeqDist)

    print(str(len(tertiaryEdges)) + " tertiary edges found")

    # print(sorted(tertiaryEdges, key=lambda tup: tup[0]))

    edges += tertiaryEdges
    distances += tertiaryDistances
    weights += tertiaryWeights

    return atoms, edges, distances, weights
