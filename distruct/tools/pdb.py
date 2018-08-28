#!/usr/bin/env python3
#####################################
#
# Filename : pdb.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Thu 11 May 2017 16:35:51 CEST
#
# Last Modified : Tue 28 Aug 2018 11:42:56 AM CEST
#
#####################################

from Bio import PDB
from Bio import Seq
from Bio.Alphabet import IUPAC


from distruct import data
# from distruct.tools import math
from . import math


def cull_atoms(atoms, structure):
    """
    Remove all atoms from a list, that are not in a structure.

    In most reference structures, there are missing atoms and residues and often there are additional ligands.
    To compute the RMSD between two structures, atom lists of identical length are needed.
    """
    for a in atoms:
        fullID = a.get_full_id()
        mID = fullID[1]
        cID = fullID[2]
        rID = fullID[3]
        aID = a.get_id()
        if mID in structure:
            if cID in structure[mID]:
                if rID in structure[mID][cID]:
                    if aID in structure[mID][cID][rID]:
                        yield a
                        pass
                    pass
                pass
            pass
        pass
    pass


# TODO redo this
def read_sequences(file):
    sequences = list()

    cifdict = PDB.MMCIF2Dict.MMCIF2Dict(file)

    typeKey = '_entity_poly.type'
    seqKey = '_entity_poly.pdbx_seq_one_letter_code_can'  # NOTE in this sequence, modified residues are identified with a one letter code

    if not isinstance(cifdict[typeKey], list):
        cifdict[typeKey] = [cifdict[typeKey]]
        cifdict[seqKey] = [cifdict[seqKey]]
        assert len(cifdict[typeKey]) == len(cifdict[seqKey])
        pass

    for s, t in zip(cifdict[seqKey], cifdict[typeKey]):
        # TODO make this safer
        a = None
        if 'peptide' in t:
            a = IUPAC.protein
        elif 'ribonucleotide' in t:
            a = IUPAC.unambiguous_rna
        else:
            raise

        sequences.append(Seq.Seq(s, a))
        pass

    return sequences


# # NOTE added chemicalDB to check for differing atom naming convention
# # TODO add missing residues from a sequence and use distances from force field at the read pdb step?
# def read_PDB(PDBCode, fileName=None, chemicalDB=data.defaultTopologyDB, assign_serial_numbers=True):
#
#     # TODO remove this and pass file instead of filename
#     import os
#     # TODO maybe handle hetero stuff if they are present in the force field
#
#     if fileName is None:
#         fileName = PDBCode + '.cif.gz'  # TODO replace this with .cif since it's the new default
#         pass
#
#     # TODO Glycin Hs are not read correctly!!! (others probably too)
#
#     structure = None
#
#     root, ext = os.path.splitext(fileName)
#
#     def _read_structure_file(structureName, structureFile, mode='.cif'):
#         nonlocal assign_serial_numbers
#         parser = None
#         if mode == '.cif':
#             assign_serial_numbers = True
#             parser = PDB.MMCIFParser()
#         else:
#             parser = PDB.PDBParser()
#             pass
#         return parser.get_structure(structureName, structureFile)
#
#     if ext == '.gz':
#         with gzip.open(fileName, 'rt') as f:
#             structure = _read_structure_file(PDBCode, f, os.path.splitext(root)[1])
#             pass
#         pass
#     elif ext == '.bz2':
#         with bzip.open(fileName, 'rt') as f:
#             structure = _read_structure_file(PDBCode, f, os.path.splitext(root)[1])
#     else:
#         structure = _read_structure_file(PDBCode, fileName, ext)
#         pass
#
#     # NOTE this is a bit clunky but it is apparently a bad idea to remove bits of an iterable while iterating over it
#     flaggedForRemoval = []
#     for model in structure:
#         for chain in model:
#             for residue in chain:
#                 if residue.id[0] != ' ':
#                     # NOTE this removes e.g. hetero residues
#                     flaggedForRemoval.append((chain.id, residue.id))
#                     assign_serial_numbers = True
#                 else:
#                     if residue.get_resname() not in chemicalDB:
#                         flaggedForRemoval.append((chain.id, residue.id))
#                         assign_serial_numbers = True
#                         pass
#                     # TODO remove this and handle termini properly
#                     if residue.has_id('OXT'):
#                         residue.detach_child('OXT')
#                         pass
#                     if chemicalDB:
#                         if residue.get_resname() in chemicalDB:
#                             for atom in residue:
#                                 if atom.get_name() not in chemicalDB[residue.get_resname()]['vertices']:
#                                     if not residue.has_id(atom.get_name()[:-1]) and atom.get_name()[:-1] in chemicalDB[residue.get_resname()]['vertices']:
#                                         # NOTE sometimes atoms in pdbfiles are called e.g. CD1 instead of CD
#                                         # print(str(residue), atom.get_name(), " not in ", chemicalDB[residue.get_resname()]['vertices'], " of " + residue.get_resname())
#                                         oldID = atom.get_id()
#                                         newID = oldID[:-1]
#                                         atom.name = newID
#                                         atom.id = newID
#                                         atom.parent.child_dict[newID] = atom.parent.child_dict[oldID]  # NOTE has_id checks in the parents child_dict. which still has the old id as key. since atoms are not entities for some reason, can't use builtin facilities to manage this
#                                         del atom.parent.child_dict[oldID]  # TODO maybe this should be done outside the loop
#                                         # print("renaming to ", atom.get_name())
#                                         # print(chemicalDB[residue.get_resname()]['bondEdges'])
#
#                                     else:
#                                         print(str(residue), atom.get_name(), " not in ", chemicalDB[residue.get_resname()]['vertices'], " of " + residue.get_resname())
#                                         # TODO handle this properly
#                                         assert False
#                                         pass
#                                     pass
#                                 pass
#                             pass
#                         else:
#                             print("unkown residue: ", str(residue), "in chain ", str(chain), " flagged for removal")
#                             flaggedForRemoval.append((chain.get_id(), residue.get_id()))
#                             assign_serial_numbers = True
#                             pass
#                     pass
#                 pass
#             pass
#         for cr in flaggedForRemoval:
#             residue = model[cr[0]][cr[1]]
#             model[cr[0]].detach_child(cr[1])
#             pass
#         flaggedForRemoval = []
#         for chain in model:
#             if len(chain) == 0:
#                 flaggedForRemoval.append(chain.id)
#                 pass
#             pass
#         for chain in flaggedForRemoval:
#             model.detach_child(chain)
#         pass
#
#     # NOTE this is necessary after removing some atoms or when loading a .cif, since the structure builder sets None for all serial numbers in that case
#     # TODO should happen when serial numbers start at 1 too!!
#     if assign_serial_numbers:
#         i = 0
#         for model in structure:
#             for chain in model:
#                 for residue in chain:
#                     for atom in residue:
#                         atom.set_serial_number(i)
#                         i += 1
#                         pass
#                     pass
#                 pass
#             pass
#         pass
#
#     print(len(list(structure.get_atoms())), ' atoms parsed')
#
#     return structure


from Bio.PDB import NeighborSearch

def get_contacts(model, cutOff=5., minSeqDist=5):
    contacts = list()

    ns = NeighborSearch(list(model.get_atoms()))
    foundPairs = ns.search_all(cutOff)

    for pair in foundPairs:
        fullAtomIDs = [
                pair[0].get_full_id(),
                pair[1].get_full_id()]
        distance = pair[1] - pair[0]
        weight = 1.

        contacts.append((fullAtomIDs, distance, weight))
        pass

    return contacts


# def get_edge(atom1, atom2):
#     edge = tuple(sorted([atom1.get_serial_number(), atom2.get_serial_number()]))
#     distance = np.sqrt(np.dot(atom1.get_coord() - atom2.get_coord(), atom1.get_coord() - atom2.get_coord()))
#     return edge, distance
#
#
# # TODO this is only for backbone dihedrals, rename
# def get_dihedral_edge(model, chainID, resID, atomID):
#     if chainID in model:
#         if resID in model[chainID] and resID - 1 in model[chainID]:
#             if model[chainID][resID].has_id(atomID) and model[chainID][resID - 1].has_id(atomID):
#                 atom1 = model[chainID][resID][atomID]
#                 atom2 = model[chainID][resID - 1][atomID]
#                 return get_edge(atom1, atom2)
#                 pass
#             pass
#         pass
#     return None, None


# TODO make this integrate better, it should take a primary sequence and a structure...
# def get_secondary_sequence_protein(fileName, topologyDB):
#     # TODO the model seems unnecessary
#     result = {}
#     # model = Bio.PDB.PDBParser().get_structure('sec', fileName)[0]
#     model = read_PDB('sec', fileName, topologyDB)[0]
#     dssp = Bio.PDB.DSSP(model, fileName)
#
#     for chain in model:
#         result[chain.get_id()] = []
#         pass
#
#     for key in dssp.keys():
#         chainID = key[0]
#         resID = key[1][1]
#         letter = dssp[key][2]
#
#         if letter == 'H':
#             result[chainID].append((letter, resID))
#         elif letter == 'E':
#             result[chainID].append((letter, resID))
#         else:
#             result[chainID].append(('-', resID))
#             pass
#         pass
#
#     #################
#     for chainID in result:
#         print([x[0] for x in result[chainID]])
#         print([x[1] for x in result[chainID]])
#         pass
#     #################
#
#     return result


# TODO integrate into Distructure
# def get_secondary_edges_protein(model, secondaryStructureSequence, topologyDB, useStructureDistances=True):
#     # TODO this works only for protein right now
#     # TODO load edges to be used from distruct.data
#     # TODO add distance from database
#
#     edges = []
#     distances = []
#     weights = []
#
#     for chain in model:
#         chainID = chain.get_id()
#         if chainID in secondaryStructureSequence:
#             # TODO check if first or last residue
#             for i, r in enumerate(secondaryStructureSequence[chainID]):
#                 print('###')
#                 print(r)
#                 resID = r[1]
#                 letter = r[0]
#                 resn = chain[resID].get_resname()
#                 # prevresn = None
#                 # nextresn = None
#                 # if resID - 1 in chain:
#                 #     prevresn = chain[resID - 1].get_resname()
#                 #     pass
#                 # if resID + 1 in chain:
#                 #     prevresn = chain[resID - 1].get_resname()
#                 #     pass
#                 edge = None
#                 distance = None
#                 # omega
#                 if resn != 'PRO':  # TODO alphabet
#                     edge, distance = get_dihedral_edge(model, chainID, resID, 'CA')
#                     if edge:
#                         print("omega:")
#                         # ############## should be 3.8 angstrom for trans
#                         A = chain[resID - 1]['CA'].get_vector()
#                         B = chain[resID - 1]['C'].get_vector()
#                         C = chain[resID]['N'].get_vector()
#                         D = chain[resID]['CA'].get_vector()
#                         dihedral = PDB.calc_dihedral(A, B, C, D)
#                         # ##############
#                         print("measured dihedral: ", dihedral*180./ np.pi)
#                         print("measured distance: ", distance)
#                         if not useStructureDistances:
#                             prevresn = chain[resID - 1].get_resname()
#                             dAB = topologyDB[prevresn]['bondEdges'][('C', 'CA')]
#                             dBC = topologyDB[resn]['bondEdges'][('-C', 'N')]
#                             dCD = topologyDB[resn]['bondEdges'][('CA', 'N')]
#                             dAC = topologyDB[prevresn]['angleEdges'][('+N', 'CA')]
#                             dBD = topologyDB[resn]['angleEdges'][('-C', 'CA')]
#                             distance = math.dist_from_dihedral(180., dAB, dBC, dCD, dAC, dBD)
#                             print("guessed dihedral: ", 180.)
#                             print("guessed distance: ", distance)
#                             pass
#                         edges.append(edge)
#                         distances.append(distance)
#                         weights.append(1.)
#                         print('omega', edge, resID)
#                         pass
#                     pass
#                 edge = None
#                 distance = None
#                 if letter == 'H':
#                     # phi
#                     edge, distance = get_dihedral_edge(model, chainID, resID, 'C')
#                     if edge:
#                         print("helix")
#                         print("phi:")
#                         # ##############
#                         A = chain[resID - 1]['C'].get_vector()
#                         B = chain[resID]['N'].get_vector()
#                         C = chain[resID]['CA'].get_vector()
#                         D = chain[resID]['C'].get_vector()
#                         dihedral = PDB.calc_dihedral(A, B, C, D)
#                         # ##############
#                         print("measured distance: ", distance)
#                         print("measured dihedral: ", dihedral * 180. / np.pi)
#                         if not useStructureDistances:
#                             dAB = topologyDB[resn]['bondEdges'][('-C', 'N')]
#                             dBC = topologyDB[resn]['bondEdges'][('CA', 'N')]
#                             dCD = topologyDB[resn]['bondEdges'][('C', 'CA')]
#                             dAC = topologyDB[resn]['angleEdges'][('-C', 'CA')]
#                             dBD = topologyDB[resn]['angleEdges'][('C', 'N')]
#                             distance = math.dist_from_dihedral(data.alpha_phi, dAB, dBC, dCD, dAC, dBD)
#                             print("guessed distance: ", distance)
#                             print("guessed dihedral: ", data.alpha_phi * 180. / np.pi)
#                             pass
#                         edges.append(edge)
#                         distances.append(distance)
#                         weights.append(1.)
#                         print('phi', edge, resID)
#                         pass
#                     # psi
#                     edge, distance = get_dihedral_edge(model, chainID, resID + 1, 'N')
#                     if edge:
#                         print("helix")
#                         print("psi:")
#                         # ##############
#                         A = chain[resID]['N'].get_vector()
#                         B = chain[resID]['CA'].get_vector()
#                         C = chain[resID]['C'].get_vector()
#                         D = chain[resID + 1]['N'].get_vector()
#                         dihedral = PDB.calc_dihedral(A, B, C, D)
#                         # ##############
#                         print("measured distance: ", distance)
#                         print("measured dihedral: ", dihedral * 180. / np.pi)
#                         if not useStructureDistances:
#                             nextresn = chain[resID + 1].get_resname()
#                             dAB = topologyDB[resn]['bondEdges'][('CA', 'N')]
#                             dBC = topologyDB[resn]['bondEdges'][('C', 'CA')]
#                             dCD = topologyDB[nextresn]['bondEdges'][('-C', 'N')]
#                             dAC = topologyDB[resn]['angleEdges'][('-C', 'CA')]
#                             dBD = topologyDB[resn]['angleEdges'][('+N', 'CA')]
#                             distance = math.dist_from_dihedral(data.alpha_psi, dAB, dBC, dCD, dAC, dBD)
#                             print("guessed distance: ", distance)
#                             print("guessed dihedral: ", data.alpha_psi * 180. / np.pi)
#                             pass
#                         edges.append(edge)
#                         distances.append(distance)
#                         weights.append(1.)
#                         print('psi', edge, resID)
#                         pass
#                     # h bond
#                     # TODO include H
#                     # TODO test this
#                     if chain.has_id(resID - 4):
#                         print('helix')
#                         print('H bond')
#                         if secondaryStructureSequence[chainID][i - 4][0] == 'H':
#                             if secondaryStructureSequence[chainID][i - 4][1] == resID - 4:
#                                 atom1 = model[chainID][resID]['N']
#                                 atom2 = model[chainID][resID - 4]['C']
#                                 edge, distance = get_edge(atom1, atom2)
#                                 print("measured distance: ", distance)
#                                 if edge:
#                                     # TODO maybe add edge between acceptor and hydrogen
#                                     if not useStructureDistances:
#                                         dAB = topologyDB[resn]['bondEdges'][('H', 'N')]
#                                         dBC = data.hbond
#                                         angle = data.hangle
#                                         distance = math.dist_from_angle(angle, dAB, dBC)
#                                         print("guessed distance: ", distance)
#                                         pass
#                                     edges.append(edge)
#                                     distances.append(distance)
#                                     weights.append(1.)
#                                     # print('H', edge, resID)
#                                     # print(distance)
#                                     pass
#                                 pass
#                             pass
#                         pass
#                 elif letter == 'E':
#                     # phi
#                     edge, distance = get_dihedral_edge(model, chainID, resID, 'C')
#                     print("measured distance: ", distance)
#                     if edge:
#                         if not useStructureDistances:
#                             dAB = topologyDB[resn]['bondEdges'][('-C', 'N')]
#                             dBC = topologyDB[resn]['bondEdges'][('CA', 'N')]
#                             dCD = topologyDB[resn]['bondEdges'][('C', 'CA')]
#                             dAC = topologyDB[resn]['angleEdges'][('-C', 'CA')]
#                             dBD = topologyDB[resn]['angleEdges'][('C', 'N')]
#                             distance = math.dist_from_dihedral(data.beta_phi, dAB, dBC, dCD, dAC, dBD)
#                             print("guessed distance: ", distance)
#                             pass
#                         edges.append(edge)
#                         distances.append(distance)
#                         weights.append(1.)
#                         pass
#                     # psi
#                     edge, distance = get_dihedral_edge(model, chainID, resID + 1, 'N')
#                     print("measured distance: ", distance)
#                     if edge:
#                         if not useStructureDistances:
#                             nextresn = chain[resID + 1].get_resname()
#                             dAB = topologyDB[resn]['bondEdges'][('CA', 'N')]
#                             dBC = topologyDB[resn]['bondEdges'][('C', 'CA')]
#                             dCD = topologyDB[nextresn]['bondEdges'][('-C', 'N')]
#                             dAC = topologyDB[resn]['angleEdges'][('-C', 'CA')]
#                             dBD = topologyDB[resn]['angleEdges'][('+N', 'CA')]
#                             distance = math.dist_from_dihedral(data.beta_psi, dAB, dBC, dCD, dAC, dBD)
#                             print("guessed distance: ", distance)
#                         edges.append(edge)
#                         distances.append(distance)
#                         weights.append(1.)
#                         pass
#                 elif letter == '-':
#                     pass
#                 else:
#                     raise
#                 pass
#             pass
#         pass
#
#     # dssp = PDB.DSSP(model, fileName)
#
#     # for key in dssp.keys():
#     #     chainID = key[0]
#     #     resID = key[1][1]
#
#     #     # omega
#     #     edge, distance = get_dihedral_edge(model, chainID, resID, 'CA')
#     #     if edge:
#     #         edges.append(edge)
#     #         if useStructureDistances:
#     #             distances.append(distance)
#     #         else:
#     #             raise NotImplementedError()
#     #         # weights.append(5.)
#     #         weights.append(1.)
#     #         pass
#     #     # if chainID in model:
#     #     #     if resID - 1 in model[chainID] and resID in model[chainID]:
#     #     #         if model[chainID][resID - 1].has_id('CA') and model[chainID][resID].has_id('CA'):
#     #     #             atom1 = model[chainID][resID - 1]['CA']
#     #     #             atom2 = model[chainID][resID]['CA']
#     #     #             edge, distance = get_edge(atom1, atom2)
#     #     #             edges.append(edge)
#     #     #             if useStructureDistances:
#     #     #                 distances.append(distance)
#     #     #             else:
#     #     #                 raise NotImplementedError()
#     #     #                 pass
#     #     #             weights.append(5.)
#     #     #             pass
#     #     #         pass
#     #     #     pass
#
#     #     if dssp[key][2] == 'H':  # alpha Helix
#     #         # phi
#     #         edge, distance = get_dihedral_edge(model, chainID, resID - 1, 'C')
#     #         if edge:
#     #             edges.append(edge)
#     #             if useStructureDistances:
#     #                 distances.append(distance)
#     #             else:
#     #                 raise NotImplementedError()
#     #             # weights.append(5.1)
#     #             weights.append(1.)
#     #             pass
#
#     #         # psi
#     #         edge, distance = get_dihedral_edge(model, chainID, resID, 'N')
#     #         if edge:
#     #             edges.append(edge)
#     #             if useStructureDistances:
#     #                 distances.append(distance)
#     #             else:
#     #                 raise NotImplementedError()
#     #             # weights.append(5.2)
#     #             weights.append(1.)
#     #             pass
#
#     #         # H bond
#     #         if chainID in model:
#     #             if resID in model[chainID] and resID + 4 in model[chainID]:
#     #                 if model[chainID][resID].has_id('N') and model[chainID][resID + 4].has_id('C'):
#     #                     # TODO take H atom into account
#     #                     atom1 = model[chainID][resID]['N']
#     #                     atom2 = model[chainID][resID]['C']
#     #                     edge, distance = get_edge(atom1, atom2)
#     #                     if edge:
#     #                         edges.append(edge)
#     #                         if useStructureDistances:
#     #                             distances.append(distance)
#     #                         else:
#     #                             raise NotImplementedError()
#     #                         # weights.append(5.3)
#     #                         weights.append(1.)
#     #                         pass
#     #                 pass
#     #             pass
#     #         pass
#     #     elif dssp[key][2] == 'E':  # beta strand
#
#     #         # phi
#     #         edge, distance = get_dihedral_edge(model, chainID, resID - 1, 'C')
#     #         if edge:
#     #             edges.append(edge)
#     #             if useStructureDistances:
#     #                 distances.append(distance)
#     #             else:
#     #                 raise NotImplementedError()
#     #             # weights.append(6.1)
#     #             weights.append(1.)
#     #             pass
#
#     #         # psi
#     #         edge, distance = get_dihedral_edge(model, chainID, resID, 'N')
#     #         if edge:
#     #             edges.append(edge)
#     #             if useStructureDistances:
#     #                 distances.append(distance)
#     #             else:
#     #                 raise NotImplementedError()
#     #             # weights.append(6.2)
#     #             weights.append(1.)
#     #             pass
#     #         pass
#     #     else:
#     #         pass
#     #     pass
#
#     return edges, distances, weights


# TODO redo this part and return ((fullid1, fullid2), distance, weight)
# # NOTE default cutOff in units of the structure (assumed angstrom)
# def get_tertiary_edges(model, cutOff=5., minSeqDist=5, getContacts=False):
#     # NOTE these would include beta sheet hbonds
#
#     tertiaryEdges = []
#     distances = []
#     weights = []
#     contacts = []
#
#     ns = PDB.NeighborSearch(list(model.get_atoms()))
#     foundPairs = ns.search_all(cutOff)
#
#     for pair in foundPairs:
#         edge = tuple(sorted((pair[0].get_serial_number(), pair[1].get_serial_number())))
#         if abs(pair[0].get_parent().get_id()[1] - pair[1].get_parent().get_id()[1]) >= minSeqDist:
#             tertiaryEdges.append(edge)
#             chainID1 = pair[0].get_parent().get_parent().get_id()
#             chainID2 = pair[1].get_parent().get_parent().get_id()
#             resID1 = pair[0].get_parent().get_id()[1]
#             resID2 = pair[1].get_parent().get_id()[1]
#             atom1 = pair[0].get_id()
#             atom2 = pair[1].get_id()
#             contacts.append(((chainID1, resID1, atom1), (chainID2, resID2, atom2)))
#             distances.append(np.sqrt(np.dot(pair[0].get_coord() - pair[1].get_coord(), pair[0].get_coord() - pair[1].get_coord())))
#             weights.append(1.)
#             pass
#         pass
#
#     if getContacts:
#         return tertiaryEdges, distances, weights, contacts
#         pass
#     return tertiaryEdges, distances, weights
