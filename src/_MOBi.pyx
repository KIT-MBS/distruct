from libc.stdint cimport uint64_t
from libc.stdint cimport uint32_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair

import numpy as np

# TODO maybe replace this? NetworKit::Point was deprecated and then undeprecated again.
cdef extern from "cpp/viz/Point.h" namespace "NetworKit":
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

# TODO add IDGPOptimizerOld

# NOTE actual MOBi interface follows


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

# TODO a more efficient MOBi::edge class may be useful at some point, but the atom IDs will have to be python objects
# cdef cppclass _edge:
#     cdef tuple fullAtomID1
#     cdef tuple fullAtomID2
#     cdef double weight
#     cdef double distance
#     pass

from Bio.PDB.Structure import Structure
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Atom import Atom

from networkit import Graph
from MOBi import data

from itertools import zip_longest

class Distructure(Structure):
    """
    Interface between hierarchical Bio.PDB.Structure representation of a molecule and a graph.

    Contains Structure, Model, Chain and Residue objects.
    Also maintains a list of contacts and corresponding edges between vertices (atoms).
    """

    def __init__(self, id, sequences = [], resIDLists = [], SSsequences = None, topologyDB=data.defaultTopologyDB):
        """
        Initialize Distructure object.

        Optionally takes primary structure (residue sequence),
        protein secondary structure (3 state)
        """
        # TODO better (markup) doc
        # TODO implement nucleic acid SS
        # TODO implement 8 state protein SS

        Structure.__init__(self, id)

        if sequences:
            chainCounter = 1
            atomCounter = 0  # NOTE atomCount starts at 0 since graph vertices do
            for sequence, resIDs in zip_longest(sequences, resIDLists, fillvalue=list()):

                chainID = self._chain_count2ID(chainCounter)
                resCounter = 1
                chainCounter += 1

                chain = Chain(chainID)

                for resID, letter in zip_longest(resIDs, sequence):
                    assert letter is not None
                    polymerType = data.polymer_type(alphabet)
                    resName = topDB['alphabets'][polymerType][letter]
                    segID = "   "  # TODO check this
                    if resID is None:
                        hetField = " "
                        iCode = " "
                        resID = [hetField, resCounter, iCode]
                        pass
                    residue = Residue(resID, resName, segID)
                    chain.add(residue)
                    resCounter += 1

                    for atomName, element in topologyDB[resName]['vertices']:
                        # atomName = vertex
                        coord = np.full(3, np.nan)
                        bFactor = 0.
                        occupancy = 1.
                        altloc = " "
                        fullName = vertex  # TODO get the correct full name from somewhere
                        serialNumber = atomCounter
                        # element = atomName = [0]  # TODO get from top db
                        atom = Atom(atomName, coord, bFactor, occupancy, altloc, fullName, serialNumber, element)
                        residue.add(atom)
                        atomCounter += 1
                        pass
                    pass
                pass
            pass
        self.graph = Graph(atomCounter, True, False)
        self.graph.setName(id)
        self._primaryContacts = list()
        self._secondaryContacts = list()
        self._tertiaryContacts = list()
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

    def _get_entity(self, fullID):
        """
        Get a child entity from the Structure from its full ID.

        Primarily intended to convert an Atoms full_id to the integer ID of the corresponding vertex.
        """

        assert self.id == fullID[0]

        e = self
        for ID in fullID[1:]:
            e = e[ID]
            pass
        return e

    def _generate_chain_primary_contacts(self, chain, useStructureDistances=False):
        """
        Generate the primary contacts for a single chain.
        """

        # TODO maybe this should be a dict
        contacts = list()

        for r in chain:
            resn = r.get_resname()
            # TODO maybe rename keys in dict
            for edgeType in ["bondEdges", "angleEdges", "improperEdges"]:
                for edge in self.topologyDB[resn][edgeType]:
                    # NOTE edges between neighboring residues in sequence are constructed
                    # NOTE with sequential sequence ids. May cause problems with non-
                    # NOTE blank insertion codes and the like (hetero field should be fine).

                    # TODO Correctly handle insertion code
                    resIDs = [r.get_id()[1], r.get_id()[1]]
                    atomIDs = list()
                    edgeAtomIDs = list()
                    atomCoords = list()
                    for i in range(len(edge)):
                        if '+' in edge[i]:
                            resIDs[i] += 1
                        elif '-' in edge[i]:
                            resIDs[i] += 1
                            pass
                        atomIDs.append(edge[i].strip('+-'))
                        pass
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
                                distance = self.topologyDB[resn][edgeType][edge]
                                pass
                            weight = 1.
                            contacts.append((edgeFullAtomIDs, distance, weight))
                        else:
                            # NOTE atoms that are present in topology are not in the structure
                            # this may happen when there are missing atoms in a parsed structure
                            # or when hydrogens are deliberately left out
                            pass
                    else:
                        # NOTE the next or previous residue of r is not present
                        # this happens, when there are missing residues in the structure
                        # or when the end of the chain is reached
                        # TODO handle termini
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
        return

    def generate_secondary_contacts(self):
        """
        Generate the contacts for backbone dihedrals and helix hydrogen bridges from SS sequence.

        Note that long range hydrogen bonds like between beta strands are tertiary contacts in this
        context.
        """

        # TODO check for protein / rna
        # TODO implement other SS elements
        raise NotImplementedError

    def set_tertiary_contacts(self, contacts):
        self._tertiaryContacts = contacts
        return

    def add_tertiary_contacts(self, contacts):
        self._tertiaryContacts += contacts
        return

    def generate_edges(self):
        # TODO improve
        # TODO add redundancy checks, setWeights just generates the edge, if it does not exist yet
        # and otherwise just sets the new weight.
        # because of this it iterates over tertiary contacts first and overwrites them later with
        # primary contact weights. this should probably not be in the shipping version!!!!!
        for contact in self._tertiaryContacts + self._primaryContacts + self._secondaryContacts:
            self.graph.setWeight(contact[0][0], contact[0][1], contact[1])
            pass

        return

    def run(self, double alpha = 1., double q = 0., uint64_t solves = 300):
        """
        Generate atomic coordinates from the supplied edges, running MaxEnt-Stress graph drawing.
        """

        # TODO work on graph directly
        edges = self.graph.edges()
        # TODO redundancy check befor this. at this point there should be no overlap between the different sets of contacts.
        distDict = {(u, v): d for ((u, v), w, d) in self._tertiaryContacts + self._primaryContacts + self._secondaryContacts}
        distances = [distDict[(u, v)] for (u, v) in edges]
        weights = [self.graph.weight(u, v) for (u, v) in edges]

        # TODO move checks here
        cdef vector[Point[double]] coord = runMaxent(self.graph.numberOfNodes(), alpha, q, solves, edges, distances, weights)

        for atom in self.get_atoms():
            atomCoord = coord[atom.get_serial_number()]
            atom.set_coord(np.array([atomCoord[0], atomCoord[1], atomCoord[2]]))
            pass
        return

    def update_serial_numbers(self):
        """
        Renumber all atoms in the structure.

        This should be called after adding new atoms to the structure.
        """
        for i, a in enumerate(self.get_atoms()):
            a.set_serial_number(i)
            pass
        return

    pass
