# distutils: language = c++
from libc.stdint cimport uint64_t
from libc.stdint cimport uint32_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair

import numpy as np

# TODO maybe replace this? NetworKit::Point was deprecated and then undeprecated again.
cdef extern from "viz/Point.h" namespace "NetworKit":
    cdef cppclass Point[T]:
        Point()
        Point(uint64_t dim)
        T& operator[](const uint64_t i) except +
        uint64_t getDimensions() except +
        pass
    pass


# TODO remove this; needed temprarily for thesis

# the wrapperwrapper since the other stuff did not work

cdef extern from "DuckingWrapper.h":
    cdef vector[Point[double]] runMaxent(
            uint64_t numNodes,
            double alpha,
            double q,
            uint64_t solves,
            vector[pair[uint64_t, uint64_t]] edges,
            vector[double] distances,
            vector[double] probabilites
            ) except +
    pass

# #cdef doublyWrappedMaxent(uint64_t numNodes, double alpha=1., double q=0., uint64_t solves=300, vector[pair[uint64_t, uint64_t]] edges=[], vector[double] distances=[], vector[double] probabilites=[]):
# def doublyWrappedMaxent(uint64_t numNodes, double alpha=1., double q=0., uint64_t solves=300, vector[pair[uint64_t, uint64_t]] edges=[], vector[double] distances=[], vector[double] probabilites=[]):
# 
#     cdef vector[Point[double]] coord = runMaxent(numNodes, alpha, q, solves, edges, distances, probabilites)
# 
#     result = []
#     for vertex in coord:
#         l = []
#         for i in range(vertex.getDimensions()):
#             l.append(vertex[i])
#             pass
#         result.append(l)
#         pass
#     return result

# cdef extern from "../src/BioMaxentStressOld.h" namespace "NetworKit":
#     cdef cppclass _BioMaxentStress "NetworKit::BioMaxentStress":
#         # TODO find a way to use actual NetworKit interface
#         # _BioMaxentStress(_Graph G, const vector[Point[double]] initialCoordinates, vector[double] distances, const uint64_t dim, double alpha) except +
#         _BioMaxentStress(
#                 uint64_t numNodes,
#                 vector[pair[uint64_t, uint64_t]] edges,
#                 vector[double] weights,
#                 vector[double] distances,
#                 double alpha,
#                 #uint64_t solvesPerAlpha,
#                 double q,
#                 uint32_t loggingFrequency) except +
#         void run(uint64_t) except +
#         # TODO diagnostics, get result, set params
#         vector[Point[double]] getCoordinates() except +
#         pass
#     pass
#
#
# cdef class BioMaxentStress:
#     cdef _BioMaxentStress *_this
#     # TODO find a way to use the (python) NetworKit Graph object and the _NetworKit _Graph object
#
#     def __cinit__(self, uint64_t numNodes, vector[pair[uint64_t, uint64_t]] edges=[], vector[double] weights=[], vector[double] distances=[], double alpha=1., double q=2., uint32_t loggingFrequency=0):
#         cdef uint64_t dim = 3
#         cdef Point[double] p = Point[double](dim)
#
#         cdef uint64_t d
#         self._this = new _BioMaxentStress(numNodes, edges, weights, distances, alpha, q, loggingFrequency)
#         return
#
#
#     def __dealloc__(self):
#         del self._this
#         return
#
#
#     def run(self, maxSolves):
#         (<_BioMaxentStress*>(self._this)).run(maxSolves)
#         return self
#
#     # TODO rename and make faster
#     def getCoordinates(self):
#         vertexCoordinates = self._this.getCoordinates()
#         result = []
#         for vertex in vertexCoordinates:
#             l = []
#             for i in range(vertex.getDimensions()):
#                 l.append(vertex[i])
#                 pass
#             result.append(l)
#             pass
#
#         return result
#
#
#     pass

# cdef extern from "../src/BioMaxentStress.h" namespace "MOBi":
#     cdef cppclass _BioMaxentStress "MOBi::BioMaxentStress":
#         # TODO find a way to use actual NetworKit interface
#         # _BioMaxentStress(_Graph G, const vector[Point[double]] initialCoordinates, vector[double] distances, const uint64_t dim, double alpha) except +
#         _BioMaxentStress(
#                 uint64_t numNodes,
#                 vector[pair[uint64_t, uint64_t]] edges,
#                 vector[double] weights,
#                 vector[double] distances,
#                 vector[Point[double]] initialCoordinates,
#                 double alpha,
#                 double q,
#                 uint32_t loggingFrequency) except +
#         void run(uint64_t maxSolves) except +
#         # TODO diagnostics, get result, set params
#         vector[Point[double]] getCoordinates() except +
#         pass
#     pass
#
#
# cdef class BioMaxentStress:
#     cdef _BioMaxentStress *_this
#     # TODO find a way to use the (python) NetworKit Graph object and the _NetworKit _Graph object
#
#     def __cinit__(self, uint64_t numNodes, vector[pair[uint64_t, uint64_t]] edges=[], vector[double] weights=[], vector[double] distances=[], initialCoordinateList=[], double alpha=1., double q=2., uint32_t loggingFrequency=0):
#         cdef uint64_t dim = 3
#         cdef Point[double] p = Point[double](dim)
#         cdef vector[Point[double]] initialCoordinates = vector[Point[double]]()
#
#         cdef uint64_t d
#         for coordinates in initialCoordinateList:
#             d = 0
#             while d < dim:
#                 p[d] = coordinates[d]
#                 d += 1
#                 pass
#             initialCoordinates.push_back(p)
#             pass
#         self._this = new _BioMaxentStress(numNodes, edges, weights, distances, initialCoordinates, alpha, q, loggingFrequency)
#         return
#
#
#     def __dealloc__(self):
#         del self._this
#         return
#
#
#     def run(self, maxSolves):
#         (<_BioMaxentStress*>(self._this)).run(maxSolves)
#         return self
#
#     # TODO rename and make faster
#     def getCoordinates(self):
#         vertexCoordinates = self._this.getCoordinates()
#         result = []
#         for vertex in vertexCoordinates:
#             l = []
#             for i in range(vertex.getDimensions()):
#                 l.append(vertex[i])
#                 pass
#             result.append(l)
#             pass
#
#         return result
#
#
#     pass

# TODO separate the Distructure class

# cdef cppclass _edge:
#     cdef tuple fullAtomID1
#     cdef tuple fullAtomID2
#     cdef double weight
#     cdef double distance
#     pass

from Bio.PDB.Structure import Structure
from Bio.PDB.Model import Model
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom

from networkit import Graph
from distruct import data

from itertools import zip_longest

class Distructure(Structure):
    """
    Interface between hierarchical Bio.PDB.Structure representation of a molecule and a graph.

    Contains Structure, Model, Chain and Residue objects.
    Also maintains a list of contacts and corresponding edges between vertices (atoms).
    """

    def __init__(self, id, sequences = [], resIDLists = [], SSsequences = None, topDB=data.defaultTopologyDB):
        """
        Initialize Distructure object.

        Optionally takes primary structure (residue sequence),
        protein secondary structure (3 state)
        """
        # TODO better (markup) doc

        # TODO redo entire initialization!!!!!!!
        # TODO use args and kwargs here for more flexibility
        # TODO separate chains and sequences
        # TODO implement terminals
        # TODO implement gaps

        # TODO implement nucleic acid SS
        # TODO implement 8 state protein SS

        self.topDB = topDB
        Structure.__init__(self, id)

        self._serial_set = False
        atomCounter = 0  # NOTE atomCount starts at 0 since graph vertices do
        if sequences:
            model = Model(0, None)
            self.add(model)

            chainCounter = 1
            for sequence, resIDs in zip_longest(sequences, resIDLists, fillvalue=list()):

                chainID = self._chain_count2ID(chainCounter)
                resCounter = 1
                chainCounter += 1

                chain = Chain(chainID)
                model.add(chain)

                for resID, letter in zip_longest(resIDs, sequence):
                    assert letter is not None
                    polymerType = data.polymer_type(sequence.alphabet)
                    if letter not in self.topDB['alphabets'][polymerType]:
                        print(letter + " is an unknown residue.")
                        print("Add Vertices manually and make sure to add enough edges to connect them all to the rest of the graph.")
                        continue
                    resName = self.topDB['alphabets'][polymerType][letter]
                    segID = "   "  # TODO check this
                    if resID is None:
                        hetField = " "
                        iCode = " "
                        resID = (hetField, resCounter, iCode)
                        pass
                    residue = Residue(resID, resName, segID)
                    chain.add(residue)
                    resCounter += 1

                    for atomName, element in self.topDB[resName]['vertices']:
                        # atomName = vertex
                        coord = np.full(3, np.nan)
                        bFactor = 0.
                        occupancy = 1.
                        altloc = " "
                        fullName = atomName  # TODO get the correct full name from somewhere
                        serialNumber = atomCounter
                        # element = atomName = [0]  # TODO get from top db
                        atom = Atom(atomName, coord, bFactor, occupancy, altloc, fullName, serialNumber, element)
                        residue.add(atom)
                        atomCounter += 1
                        pass
                    pass
                pass
            self._serial_set = True
            pass
        self.graph = Graph(atomCounter, True, False)
        self.graph.setName(id)
        self._primaryContacts = list()
        self._secondaryContacts = list()
        self._tertiaryContacts = list()
        self._edgesSet = False
        return

    def _chain_count2ID(self, count):
        ID = ''

        base26 = []
        while count:
            base26.append(count % 26)
            count = count // 26
            pass

        for i in range(len(base26) - 1):
            if base26[i] <= 0:
                base26[i] += 26
                base26[i+1] -= 1
                pass
            pass
        for x in base26:
            if x > 0:
                ID += chr(ord('A') - 1 + x)
                pass
            pass
        return ID

    # TODO use biopython unfold entities instead
    def _get_entity(self, fullID):
        """
        Get a child entity from the Structure from its full ID.

        Primarily intended to convert an Atoms full_id to the integer ID of the corresponding vertex.
        """

        e = self
        for ID in fullID[1:]:
            if e.get_level() == 'R':
                e = e[ID[0]]  # NOTE for biopython weirdness reasons, they don't return the atom_id in the full_id, but (atom_name, altloc)
            else:
                e = e[ID]
                pass
            pass
        return e

    def _generate_chain_primary_contacts(self, chain, useStructureDistances=False):
        """
        Generate the primary contacts for a single chain.
        """

        # TODO maybe this should be a dict
        contacts = list()

        # TODO this is to handle insertion codes, test it works correctly.
        chainResidues = list(chain.get_residues())
        for idx, r in enumerate(chainResidues):
            resn = r.get_resname()
            if resn not in self.topDB:
                print(resn + ' ' + str(r.get_id()) + " in " + str(chain.get_id()) + " is an unknown residue.")
                print("Add Vertices manually and make sure to add enough edges to connect them all to the rest of the graph.")
                continue

            for edgeType in ["bondEdges", "angleEdges", "improperEdges"]:
                for edge in self.topDB[resn][edgeType]:
                    # NOTE edges between neighboring residues in sequence are constructed
                    # NOTE with sequential sequence ids. May cause problems with non-
                    # NOTE blank insertion codes and the like (hetero field should be fine).

                    # TODO clean up, functionify

                    # TODO Correctly handle insertion code and gaps
                    # TODO this is to handle insertion codes, test it works correctly.
                    # resIDs = [r.get_id()[1], r.get_id()[1]]
                    resIndices = [idx, idx]

                    for i in range(len(edge)):
                        if '+' in edge[i]:
                            resIndices[i] += 1
                        elif '-' in edge[i]:
                            resIndices[i] -= 1
                            pass
                        pass
                    # NOTE check indices are not out of bounds and residues are adjacent in sequence
                    if resIndices[0] < 0:
                        continue
                    if resIndices[0] >= len(chainResidues):
                        continue
                    if resIndices[1] < 0:
                        continue
                    if resIndices[1] >= len(chainResidues):
                        continue

                    resIDs = [chainResidues[resIndices[0]].get_id(), chainResidues[resIndices[1]].get_id()]
                    if abs(resIDs[0][1] - resIDs[1][1]) > 1:
                        continue

                    # for i in range(len(edge)):
                    #     if '+' in edge[i]:
                    #         resIDs[i] += 1
                    #     elif '-' in edge[i]:
                    #         resIDs[i] -= 1
                    #         pass
                    #     atomIDs.append(edge[i].strip('+-'))
                    #     pass
                    atomIDs = list()
                    atomCoords = list()
                    atomIDs = [x.strip('+-') for x in edge]
                    if resIDs[0] in chain and resIDs[1] in chain:
                        if chain[resIDs[0]].has_id(atomIDs[0]) and chain[resIDs[1]].has_id(atomIDs[1]):
                            edgeFullAtomIDs = [
                                    chain[resIDs[0]][atomIDs[0]].get_full_id(),
                                    chain[resIDs[1]][atomIDs[1]].get_full_id()]
                            atoms = [
                                    chain[resIDs[0]][atomIDs[0]],
                                    chain[resIDs[1]][atomIDs[1]]]

                            distance = None
                            # TODO check how coordinates are initialized, when unknown
                            if useStructureDistances and (atoms[1] - atoms[0] > 0.1):
                                distance = atoms[1] - atoms[0]
                            else:
                                distance = self.topDB[resn][edgeType][edge]
                                pass
                            weight = 1.
                            contacts.append((edgeFullAtomIDs, distance, weight))
                        else:
                            # NOTE atoms that are present in topology are not in the structure
                            # this may happen when there are missing atoms in a parsed structure
                            # or when hydrogens are deliberately left out
                            # print("one of the atoms is missing in the structure")
                            # print(chain.get_id())
                            # print(resIDs)
                            # print(atomIDs)
                            pass
                    else:
                        # NOTE the next or previous residue of r is not present
                        # this happens, when there are missing residues in the structure
                        # or when the end of the chain is reached
                        # TODO handle termini
                        # print("one of the residues is missing in the structure")
                        # print("TODO this also warns at terminals, which is dumb")
                        # print(chain.get_id())
                        # print(resIDs)
                        pass
                    pass
                pass
            pass

        return contacts

    def generate_primary_contacts(self):
        """
        Generate the contacts for bonds, angles and improper/fixed dihedrals.
        """

        contacts = list()
        for c in self[0]:
            contacts += self._generate_chain_primary_contacts(c)
            pass

        self._primaryContacts = contacts
        self._edgesSet = False
        return

    def generate_secondary_contacts(self):
        """
        Generate the contacts for backbone dihedrals and helix hydrogen bridges from SS sequence.

        Note that long range hydrogen bonds like between beta strands are tertiary contacts in this
        context.
        """

        # TODO check for protein / rna
        # TODO implement other SS elements
        self._edgesSet = False
        raise NotImplementedError

    def set_tertiary_contacts(self, contacts):
        """
        Check supplied edges for missing atoms and warn, then add only the useful ones.
        """
        self._tertiaryContacts = list()
        for contact in contacts:
            fullID1 = contact[0][0]
            fullID2 = contact[0][1]
            distance = contact[1]
            weight = contact[2]

            if fullID1[-3] in self[0] and fullID2[-3] in self[0]:
                chain1 = self[0][fullID1[-3]]
                chain2 = self[0][fullID2[-3]]
                if fullID1[-2] in chain1 and fullID2[-2] in chain2:
                    residue1 = chain1[fullID1[-2]]
                    residue2 = chain2[fullID2[-2]]

                    atomID1 = fullID1[-1]
                    atomID2 = fullID2[-1]
                    if isinstance(atomID1, str):
                        atomID1 = (atomID1, ' ')
                    if isinstance(atomID2, str):
                        atomID2 = (atomID2, ' ')

                    if atomID1[0] in residue1 and atomID2[0] in residue2:
                        fullAtomIDs = [
                                residue1[atomID1[0]].get_full_id(),
                                residue2[atomID2[0]].get_full_id()]
                        self._tertiaryContacts.append((fullAtomIDs, distance, weight))
                    else:
                        # TODO implement verbosity
                        # print("at least one atom in the edge between")
                        # print(chain1.get_id() + '(' + str(residue1.get_id()[1]) + ' ' + residue1.get_resname() + ')' + fullID1[4][0])
                        # print("and")
                        # print(chain2.get_id() + '(' + str(residue2.get_id()[1]) + ' ' + residue2.get_resname() + ')' + fullID2[4][0])
                        # print("is missing in the structure (dist = " + str(contact[1]) + ')')
                        # print([a.get_id() for a in residue1.get_atoms()])
                        # print([a.get_id() for a in residue2.get_atoms()])
                        # print('')
                        pass
                else:
                    print("at least one residue in", contact, "is missing in the structure")
                    pass
            else:
                print("at least one chain in", contact, "is missing in the structure")
                pass
            pass
        self._edgesSet = False
        return

    # def add_tertiary_contacts(self, contacts):
    #     self._tertiaryContacts += contacts

    #     self._edgesSet = False
    #     return

    def generate_edges(self):
        # TODO improve
        # TODO add redundancy checks, setWeights just generates the edge, if it does not exist yet
        # and otherwise just sets the new weight.
        # because of this it iterates over tertiary contacts first and overwrites them later with
        # primary contact weights. this should probably not be in the shipping version!!!!!

        # TODO better way to do the distances
        print("generating edges...")

        # TODO set this to false when adding new atoms
        if not self._serial_set:
            atomCounter = 0
            for atom in self.get_atoms():
                atom.set_serial_number(atomCounter)
                atomCounter += 1
                pass
            self._serial_set = True
            self.graph = Graph(atomCounter, True, False)
            pass

        self.distDict = dict()
        for contact in self._tertiaryContacts + self._primaryContacts + self._secondaryContacts:
            atom1, atom2 = contact[0]
            distance = contact[1]
            weight = contact[2]

            vertex1 = self._get_entity(atom1).get_serial_number()
            vertex2 = self._get_entity(atom2).get_serial_number()

            edge = tuple(sorted((vertex1, vertex2), reverse=True))

            self.distDict[edge] = distance
            self.graph.setWeight(vertex1, vertex2, weight)
            pass
        self._edgesSet = True
        self.graph.indexEdges()

        return

    def run(self, double alpha = 1., double q = 0., uint64_t solves = 300):
        """
        Generate atomic coordinates from the supplied edges, running MaxEnt-Stress graph drawing.
        """

        if not self._edgesSet:
            self.generate_edges()
            self._edgesSet = True
            pass

        # TODO move all the checks from the wrapper here
        # TODO work on graph directly
        edges = self.graph.edges()
        distances = [self.distDict[e] for e in edges]
        weights = [self.graph.weight(u, v) for (u, v) in edges]

        # NOTE connectivity check

        print("running MaxEnt stress")
        cdef vector[Point[double]] coord = runMaxent(self.graph.numberOfNodes(), alpha, q, solves, edges, distances, weights)

        for atom in self.get_atoms():
            atomCoord = coord[atom.get_serial_number()]
            atom.set_coord(np.array([atomCoord[0], atomCoord[1], atomCoord[2]]))
            pass

        self.error2()

        for atom in self.get_atoms():
            atom.set_bfactor(self.nodeErrors2[atom.get_serial_number()])
            pass
        return

    # TODO put these in cpp
    # TODO write tests
    def error(self):
        atoms = list(self[0].get_atoms())
        self.edgeErrors = np.zeros(self.graph.numberOfEdges())
        self.nodeErrors = np.zeros(self.graph.numberOfNodes())

        def edgeError(u, v, weight, id):
            assert u == atoms[u].get_serial_number()
            assert v == atoms[v].get_serial_number()
            edge = tuple(sorted((u, v), reverse=True))
            d = self.distDict[edge]
            self.edgeErrors[id] = np.abs(d - (atoms[u] - atoms[v]))
            return

        self.graph.forEdges(edgeError)

        def nodeError(u):
            neighborErrors = np.zeros(self.graph.degree(u))
            i = 0
            def neighborError(u, v, weight, id):
                edge = tuple(sorted((u, v), reverse=True))
                d = self.distDict[edge]
                neighborErrors[i] = np.abs(d - (atoms[u] - atoms[v]))
                return
            self.graph.forEdgesOf(u, neighborError)
            self.nodeErrors[u] = np.sum(neighborErrors) / len(neighborErrors)
            return
        self.graph.forNodes(nodeError)

        return np.sum(self.edgeErrors) / len(self.edgeErrors)

    # TODO put these in cpp
    def error2(self):
        atoms = list(self[0].get_atoms())
        self.edgeErrors2 = np.zeros(self.graph.numberOfEdges())
        self.nodeErrors2 = np.zeros(self.graph.numberOfNodes())

        def edgeError2(u, v, weight, id):
            assert u == atoms[u].get_serial_number()
            assert v == atoms[v].get_serial_number()
            edge = tuple(sorted((u, v), reverse=True))
            d = self.distDict[edge]
            self.edgeErrors2[id] = (d - (atoms[u] - atoms[v]))*(d - (atoms[u] - atoms[v]))
            return

        self.graph.forEdges(edgeError2)

        def nodeError2(u):
            neighborErrors2 = np.zeros(self.graph.degree(u))
            i = 0
            def neighborError2(u, v, weight, id):
                edge = tuple(sorted((u, v), reverse=True))
                d = self.distDict[edge]
                neighborErrors2[i] = (d - (atoms[u] - atoms[v]))*(d - (atoms[u] - atoms[v]))
                return
            self.graph.forEdgesOf(u, neighborError2)
            self.nodeErrors2[u] = np.sqrt(np.sum(neighborErrors2) / len(neighborErrors2))
            return
        self.graph.forNodes(nodeError2)

        return np.sqrt(np.sum(self.edgeErrors2) / len(self.edgeErrors2))


    # TODO put these in cpp
    def stress(self):
        raise NotImplementedError

    pass
