
/*************************************
 *
 *  Filename : BioMaxentStress.h
 *
 *  Projectname : MOBi
 *
 *  Author : Oskar Taubert, Michael Wegner
 *
 *  Creation Date : Wed 23 May 2018 05:27:01 PM CEST
 *
 *  Last Modified : Wed 23 May 2018 05:28:46 PM CEST
 *
 * *************************************/

#ifndef BIOMAXENTSTRESS_H
#define BIOMAXENTSTRESS_H



#include <cpp/viz/GraphLayoutAlgorithm.h>
#include <cpp/numerics/LinearSolver.h>
#include <cpp/algebraic/CSRMatrix.h>
#include <cpp/viz/Octree.h>

//#include <NetworKit/graph/BFS.h>
#include <cpp/distance/Dijkstra.h>
#include <cpp/distance/AlgebraicDistance.h>

#include <cpp/viz/FruchtermanReingold.h>
#include <cpp/io/LineFileReader.h>
#include <cpp/auxiliary/StringTools.h>

#include <sys/time.h>

#include <memory>

#include "aux/sig.h"

namespace NetworKit {

typedef std::vector<Vector> CoordinateVector; // more meaningful and shorter name for coordinates stored by dimension

/**
 * @ingroup viz
 *
 * Implementation of MaxentStress by Gansner et al. using a Laplacian system solver for the application of protein structure determination.
 *
 * @see Gansner, Emden R., Yifan Hu, and Steve North. "A maxent-stress model for graph layout." Visualization and Computer Graphics, IEEE Transactions on 19, no. 6 (2013): 927-940.
 */
class BioMaxentStress : public GraphLayoutAlgorithm<double> {
public:
	/**
	 * Creates a BioMaxentStress object for graph @a G. The algorithm embeds the graph in @a dim dimensional Euclidean
	 * space using the Laplacian solver specified by @a solver. If @a fastComputation is set to true, the algorithm
	 * will converge faster but the quality will usually be worse. By default, @a fastComputation is false.
	 * @param G The graph.
	 * @param dim Number of dimensions (usually we want 3 dimensions)
	 * @param solver A linear solver capable of solving Laplacian linear systems.
	 * @param fastComputation
	 */
	BioMaxentStress(const Graph& G, const count dim, LinearSolver<CSRMatrix>& solver, bool fastComputation=false);

	/**
	 * Creates a BioMaxentStress object for graph @a G. The algorithm embeds the graph in @a dim dimensional Euclidean
	 * space using the Laplacian solver specified by @a solver. If @a fastComputation is set to true, the algorithm
	 * will converge faster but the quality will usually be worse. By default, @a fastComputation is false.
	 * The vector @a probability gives a probability for each edge in @a G that states how probable it is that the
	 * actual distance lies in the distance interval.
	 * @param G The graph.
	 * @param dim Number of dimensions (usually we want 3 dimensions)
	 * @param solver A linear solver capable of solving Laplacian linear systems.
	 * @param probability Probability for each edge (by edgeId)
	 * @param fastComputation
	 */
	BioMaxentStress(const Graph& G, const count dim, LinearSolver<CSRMatrix>& solver, std::vector<double>& probability, bool fastComputation=false);

	/**
	 * Creates a BioMaxentStress object for graph @a G. The algorithm embeds the graph in @a dim dimensional Euclidean
	 * space using the Laplacian solver specified by @a solver. The @a coordinates vector stores initial coordinates
	 * for each vertex.  If @a fastComputation is set to true, the algorithm will converge faster but the quality will
	 * usually be worse. By default, @a fastComputation is false.
	 * @param G The graph.
	 * @param dim Number of dimensions (usually we want 3 dimensions)
	 * @param coordinates Initial coordinates of the vertices of @a G.
	 * @param solver A linear solver capable of solving Laplacian linear systems.
	 * @param fastComputation
	 */
	BioMaxentStress(const Graph& G, const count dim, const std::vector<Point<double>>& coordinates, LinearSolver<CSRMatrix> &solver, bool fastComputation=false);

	/**
	 * Creates a BioMaxentStress object for graph @a G. The algorithm embeds the graph in @a dim dimensional Euclidean
	 * space using the Laplacian solver specified by @a solver. The @a coordinates vector stores initial coordinates
	 * for each vertex.  If @a fastComputation is set to true, the algorithm will converge faster but the quality will
	 * usually be worse. By default, @a fastComputation is false.
	 * The vector @a probability gives a probability for each edge in @a G that states how probable it is that the
	 * actual distance lies in the distance interval.
	 * @param G The graph.
	 * @param dim Number of dimensions (usually we want 3 dimensions)
	 * @param coordinates Initial coordinates of the vertices of @a G.
	 * @param solver A linear solver capable of solving Laplacian linear systems.
	 * @param probability Probability for each edge (by edgeId)
	 * @param fastComputation
	 */
	BioMaxentStress(const Graph& G, const count dim, const std::vector<Point<double>>& coordinates, LinearSolver<CSRMatrix> &solver, std::vector<double>& probability, bool fastComputation=false);


	/** Default destructor. */
	virtual ~BioMaxentStress() = default;

	/** Computes a graph drawing according to the Maxent-Stress model. */
	virtual void run();

	/**
	 * Set parameter @a q.
	 * @param q
	 */
	void setQ(double q) {
		this->q = q;
	}

	/**
	 * Set parameter @a alpha.
	 * @param alpha
	 */
	void setAlpha(double alpha) {
		this->alpha = alpha;
	}

	/**
	 * Set parameter @a alphaReduction.
	 * @param alphaReduction
	 */
	void setAlphaReduction(double alphaReduction) {
		this->alphaReduction = alphaReduction;
	}

	/**
	 * Set parameter @a finalAlpha.
	 * @param finalAlpha
	 */
	void setFinalAlpha(double finalAlpha) {
		this->finalAlpha = finalAlpha;
	}

	/**
	 * Set convergence threshold used by the maxent-stress algorithm.
	 * @param convThreshold
	 */
	void setConvergenceThreshold(double convThreshold) {
		this->convThreshold = convThreshold * convThreshold;
	}

	/**
	 * Set the maximum number of solves per alpha.
	 * @param maxSolvesPerAlpha
	 */
	void setMaxSolvesPerAlpha(count maxSolvesPerAlpha) {
		this->maxSolvesPerAlpha = maxSolvesPerAlpha;
	}

private:
	/**
	 * Reference to the linear solver to use during the maxent-stress algorithm.
	 */
	LinearSolver<CSRMatrix>& solver;

	/** Parameters of the MaxentStress model **/
	double q, alpha, alphaReduction, finalAlpha, convThreshold;

	/** Specifies whether initial coordinates have been provided in the constructor */
	bool coordinatesProvided;

	/** Specifies whether probabilities (for weighting the distances) are provided in the constructor */
	bool probabilityProvided;

	/** Defines whether the algorithm stops when converged on a higher than the lowest level. This saves some time but
	 *  usually leads to slightly worse results.
	 */
	bool fastComputation;

	/** Maximum number of solves for the same value of alpha **/
	count maxSolvesPerAlpha;

	/** Probability for edge having the given distance **/
	std::vector<double> probability;

	/** points of vertices are in R^{dim} */
	count dim;

	/**
	 * Determines whether the run() method has already been called.
	 */
	bool hasRun;

	/**
	 * Checks whether the MaxentStress algorithm converged, i.e. ||newCoords - oldCoords|| / ||oldCoords|| < convThreshold.
	 * @param newCoords The new coordinates computed in the current round of the algorithm.
	 * @param oldCoords The coordinates from the previous round of the algorithm.
	 * @return @code True when converged, otherwise @code false.
	 */
	bool isConverged(const CoordinateVector& newCoords, const CoordinateVector& oldCoords);

	/**
	 * Create a weighted Laplacian matrix from G and setup the solver for this matrix.
	 */
	void setupWeightedLaplacianMatrix();

	/**
	 * Computes the vector L_{w,d}*x where x is stored in @a coordinates and stores the result in @a rhs.
	 * @param coordinates The coordinate vector (x in the thesis)
	 * @param rhs The right-hand side that stores the result of the matrix-vector multiplication.
	 */
	void computeCoordinateLaplacianTerm(const CoordinateVector& coordinates, CoordinateVector& rhs);

	/**
	 * Computes the repulsive forces according to Equation (8) in Gansner et al. exactly (i.e. using no approximation).
	 * @param coordinates The current coordinates of the vertices.
	 * @param b Repulsive force vector to compute
	 */
	CoordinateVector computeRepulsiveForces(const CoordinateVector& coordinates, CoordinateVector& b) const;

	/**
	 * Approximates the repulsive forcse by means of an octree (Barnes and Hut).
	 * @param coordinates The current coordinates of the vertices.
	 * @param octree Octree for Barnes-Hut approximation.
	 * @param theta Parameter for Barnes-Hut cell-opening criterion.
	 * @param b Repulsive force vector to compute.
	 */
	void approxRepulsiveForces(const CoordinateVector& coordinates, const Octree<double>& octree, const double theta, CoordinateVector& b) const;

	/**
	 * Initializes the @a coordinates corresponding to vertices in the Graph to a random point in d-dimensional space 50^d pixel.
	 * @param coordinates The coordinates of the graph vertices.
	 */
	void randomInitCoordinates(CoordinateVector& coordinates) const;

	/**
	 * Initialized the @a coordinates corresponding to vertices in the Graph acoording to the random sphere placement algorithm.
	 * @param coordinates The coordinates of the graph vertices.
	 */
	void randomSphereCoordinates(CoordinateVector& coordinates) const;

	/**
	 * Computes the squared distance ||c_i - c_j||^2 between @a coordinates @a i and @a j
	 * @param coordinates
	 * @param i
	 * @param j
	 */
	double squaredDistance(const CoordinateVector& coordinates, const index i, const index j) const;


	/**
	 * Computes the distance ||c_i - c_j|| between @a coordinates @a i and @a j.
	 * @param coordinates
	 * @param i
	 * @param j
	 */
	inline double distance(const CoordinateVector& coordinates, const index i, const index j) const {
		return sqrt(squaredDistance(coordinates, i, j));
	}

	/**
	 * Computes the squared distance ||c1_i - c2_j||^2 between coordinate @a i from @a coordinates1 and coordinate @a j from @a coordinates 2.
	 * @param coordinates1
	 * @param coordinates2
	 * @param i
	 * @param j
	 */
	double squaredDistance(const CoordinateVector& coordinates1, const CoordinateVector& coordinates2, const index i, const index j) const;

	/**
	 * Computes the squared distance ||c1_i - c2_j|| between coordinate @a i from @a coordinates1 and coordinate @a j from @a coordinates 2.
	 * @param coordinates1
	 * @param coordinates2
	 * @param i
	 * @param j
	 */
	inline double distance(const CoordinateVector& coordinates1, const CoordinateVector& coordinates2, const index i, const index j) const {
		return sqrt(squaredDistance(coordinates1, coordinates2, i, j));
	}

	/**
	 * Computes the squared length ||c_i||^2 of coordinate @a i in @a coordinates.
	 * @param coordinates
	 * @param i
	 */
	double squaredLength(const CoordinateVector& coordinates, const index i) const;

	/**
	 * Computes the length ||c_i|| of coordinate @a i in @a coordinates.
	 * @param coordinates
	 * @param i
	 */
	inline double length(const CoordinateVector& coordinates, const index i) const {
		return sqrt(squaredLength(coordinates, i));
	}



	/**
	 * Weighting factor, Gansner et al. propose 1/(edgeWeight^2). If no probabilities have been provided we use uniform
	 * weights. Otherwise we use the weight as a penalty function s.t. high probabilities lead to higher weights.
	 * @param edgeWeight The graph-theoretic distance between two vertices (i.e. the edge weight).
	 * @param edgeId The id of the corresponding edge in the graph.
	 * @return The corresponding weighting factor.
	 */
	inline double weightingFactor(double edgeWeight, index edgeId) const {
        return probabilityProvided? probability[edgeId] : 1.0;
		//return probabilityProvided? 1+5*exp(-5*(1-probability[edgeId])) : 1.0;
	}

	/**
	 * Returns the sign of @a value. The sign of 0.0 is defined to be positive.
	 * @param value
	 */
	inline double sign(const double value) const {
		return (value >= 0.0) - (value < 0.0);
	}

	/**
	 * Returns the point at index @a i stored in @a coordinates.
	 * @param coordinates
	 * @param i
	 */
	inline Point<double> getPoint(const CoordinateVector& coordinates, index i) const {
		Point<double> p(coordinates.size());
		for (index d = 0; d < p.getDimensions(); ++d) {
			assert(coordinates[d].getDimension() > i);
			p[d] = coordinates[d][i];
		}

		return p;
	}
};









} /* namespace NetworKit */
#endif /* BIOMAXENTSTRESS_H */
