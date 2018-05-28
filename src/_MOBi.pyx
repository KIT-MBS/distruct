from libc.stdint cimport uint64_t
from libc.stdint cimport uint32_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair

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

#cdef doublyWrappedMaxent(uint64_t numNodes, double alpha=1., double q=0., uint64_t solves=300, vector[pair[uint64_t, uint64_t]] edges=[], vector[double] distances=[], vector[double] probabilites=[]):
def doublyWrappedMaxent(uint64_t numNodes, double alpha=1., double q=0., uint64_t solves=300, vector[pair[uint64_t, uint64_t]] edges=[], vector[double] distances=[], vector[double] probabilites=[]):

    cdef vector[Point[double]] coord = runMaxent(numNodes, alpha, q, solves, edges, distances, probabilites)

    result = []
    for vertex in coord:
        l = []
        for i in range(vertex.getDimensions()):
            l.append(vertex[i])
            pass
        result.append(l)
        pass
    return result

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

# TODO separate the Distruct class
from Bio.PDB.Structure import Structure
from networkit import Graph

class Distruct(Structure):
    """
    Interface between hierarchical Bio.Structure representation of a molecule and a graph.

    Contains Structure, Model, Chain and Residue objects.
    Also contains a list of edges between vertices (atoms).
    """

    def __init__(self, id, sequences = None, SSsequences = None):
        """
        Initialize Distruct object.

        Optionally takes primary structure (residue sequence),
        protein secondary structure (3 state)
        """
        # TODO better (markup) doc
        # TODO implement nucleic acid SS
        # TODO implement 8 state protein SS

        Structure.__init__(self, id)
        self.graph = Graph(0, True, False)
        if sequences is not None:
            # TODO
            raise NotImplementedError
        return

    # TODO Distructure builder

    def generate_primary_edges():
        """
        Generate the edges for bonds, angles and improper/fixed dihedrals.
        """
        raise NotImplementedError

    def generate_secondary_edges():
        """
        Generate the edges for backbone dihedrals and helix hydrogen bridges from SS.
        """
        raise NotImplementedError

    def generate_tertiary_edges():
        """
        Generate tertiary edges from supplied contacts
        """
        raise NotImplementedError

    def run():
        """
        Generate atomic coordinates from the supplied edges, running MaxEnt-Stress graph drawing.
        """
        raise NotImplementedError

    pass
