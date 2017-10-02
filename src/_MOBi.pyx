import networkit as nwk

from libc.stdint cimport uint64_t
from libcpp.vector cimport vector
from libcpp.utility cimport pair

# TODO maybe replace this? NetworKit::Point was deprecated and then undeprecated again.
cdef extern from "NetworKit/viz/Point.h" namespace "NetworKit":
    cdef cppclass Point[T]:
        Point()
        Point(uint64_t dim)
        T& operator[](const uint64_t i) except +
        uint64_t getDimensions() except +
        pass
    pass


# cdef extern from "NetworKit/viz/GraphLayoutAlgorithm.h":
#     cdef cppclass _GraphLayoutAlgorithm "NetworKit::GraphLayoutAlgorithm"[T]:
#         _GraphLayoutAlgorithm(_Graph, uint64_t) except +
#         vector[Point[double]] getCoordinates() except +
#         pass
#     pass
#
#
# cdef class GraphLayoutAlgorithm:
#     """Abstract base class for graph embedding algorithms"""
#     cdef _GraphLayoutAlgorithm[double] *_this
#     cdef object _G
#
#     def __init__(self, *args, **kwargs):
#         if type(self) == GraphLayoutAlgorithm:
#             raise RuntimeError("Error: don't use GraphLayoutAlgorithm directly")
#         return
#
#     def __dealloc__(self):
#         self._G = None
#         return
#
#     def getCoordinates(self):
#         """returns vertexCoordinates as list of dim dimensional vectors"""
#         if self._this == NULL:
#             raise RuntimeError("Error: GraphLayoutAlgorithm object not properly initialized")
#         pointCoord = self._this.getCoordinates()
#         coordList = []
#         for coord in pointCoord:
#             l = []
#             dim = coord.getDimensions()
#             d = 0
#             while d < dim:
#                 l.append(coord[d])
#                 d += 1
#                 pass
#             coordList.append(l)
#             pass
#         return
#
#     pass


# NOTE actual MOBi interface follows


cdef extern from "../src/BioMaxentStress.h" namespace "MOBi":
    cdef cppclass _BioMaxentStress "MOBi::BioMaxentStress":
        # TODO find a way to use actual NetworKit interface
        # _BioMaxentStress(_Graph G, const vector[Point[double]] initialCoordinates, vector[double] distances, const uint64_t dim, double alpha) except +
        _BioMaxentStress(
                uint64_t numNodes,
                vector[pair[uint64_t, uint64_t]] edges,
                vector[double] weights,
                vector[double] distances,
                vector[Point[double]] initialCoordinates) except +
        void run(uint64_t maxSolves) except +
        # TODO diagnostics, get result, set params
        vector[Point[double]] getCoordinates() except +
        pass
    pass


cdef class BioMaxentStress:
    cdef _BioMaxentStress *_this
    # TODO find a way to use the (python) NetworKit Graph object and the _NetworKit _Graph object

    def __cinit__(self, uint64_t numNodes, vector[pair[uint64_t, uint64_t]] edges=[], vector[double] weights=[], vector[double] distances=[], initialCoordinateList=[]):
        # TODO
        cdef uint64_t dim = 3
        cdef Point[double] p = Point[double](dim)
        cdef vector[Point[double]] initialCoordinates = vector[Point[double]]()

        cdef uint64_t d
        for coordinates in initialCoordinateList:
            d = 0
            while d < dim:
                p[d] = coordinates[d]
                d += 1
                pass
            initialCoordinates.push_back(p)
            pass
        self._this = new _BioMaxentStress(numNodes, edges, weights, distances, initialCoordinates)
        return


    def __dealloc__(self):
        del self._this
        return


    def run(self, maxSolves):
        (<_BioMaxentStress*>(self._this)).run(maxSolves)
        return self

    # TODO rename and make faster
    def getCoordinates(self):
        vertexCoordinates = self._this.getCoordinates()
        result = []
        for vertex in vertexCoordinates:
            l = []
            for i in range(vertex.getDimensions()):
                l.append(vertex[i])
                pass
            result.append(l)
            pass
        
        return result


    pass
