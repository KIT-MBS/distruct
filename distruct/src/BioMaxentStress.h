
/*************************************
 *
 *  Filename : BioMaxentStress.h
 *
 *  Projectname : diSTruct
 *
 *  Author : Oskar Taubert, Michael Wegner
 *
 *  Creation Date : Wed 23 May 2018 05:27:01 PM CEST
 *
 *  Last Modified : Fri 29 Mar 2019 04:07:03 PM CET
 *
 * *************************************/

#ifndef BIOMAXENTSTRESS_H
#define BIOMAXENTSTRESS_H

#include <viz/GraphLayoutAlgorithm.h>
#include <numerics/LinearSolver.h>
#include <algebraic/CSRMatrix.h>
#include <viz/Octree.h>

//provides keyboard interrupt handling
#include "aux/sig.h"

namespace diSTruct {

    /**
     * Implementation of MaxentStress by Gansner et al. using a Laplacian system solver for the application of protein structure determination.
     *
     * @see Gansner, Emden R., Yifan Hu, and Steve North. "A maxent-stress model for graph layout." Visualization and Computer Graphics, IEEE Transactions on 19, no. 6 (2013): 927-940.
     */
    class BioMaxentStress : public NetworKit::GraphLayoutAlgorithm<double> {
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
            BioMaxentStress(const NetworKit::Graph& G, const uint64_t dim, NetworKit::LinearSolver<NetworKit::CSRMatrix>& solver, bool fastComputation=false);

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
            BioMaxentStress(const NetworKit::Graph& G, const uint64_t dim, NetworKit::LinearSolver<NetworKit::CSRMatrix>& solver, std::vector<double>& probability, bool fastComputation=false);

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
            BioMaxentStress(const NetworKit::Graph& G, const uint64_t dim, const std::vector<NetworKit::Point<double>>& coordinates, NetworKit::LinearSolver<NetworKit::CSRMatrix> &solver, bool fastComputation=false);

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
            BioMaxentStress(const NetworKit::Graph& G, const uint64_t dim, const std::vector<NetworKit::Point<double>>& coordinates, NetworKit::LinearSolver<NetworKit::CSRMatrix> &solver, std::vector<double>& probability, bool fastComputation=false);


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
            void setMaxSolvesPerAlpha(uint64_t maxSolvesPerAlpha) {
                this->maxSolvesPerAlpha = maxSolvesPerAlpha;
            }

        private:
            /**
             * Reference to the linear solver to use during the maxent-stress algorithm.
             */
            NetworKit::LinearSolver<NetworKit::CSRMatrix>& solver;

            /** Parameters of the MaxentStress model **/
            double q, alpha, alphaReduction, finalAlpha, convThreshold;

            /** Defines whether the algorithm stops when converged on a higher than the lowest level. This saves some time but
             *  usually leads to slightly worse results.
             */
            bool fastComputation;

            /** Maximum number of solves for the same value of alpha **/
            uint64_t maxSolvesPerAlpha;

            /** Probability for edge having the given distance **/
            std::vector<double> probability;

            /** points of vertices are in R^{dim} */
            uint64_t dim;

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
            bool isConverged(const std::vector<NetworKit::Vector>& newCoords, const std::vector<NetworKit::Vector>& oldCoords);

            /**
             * Create a weighted Laplacian matrix from G and setup the solver for this matrix.
             */
            void setupWeightedLaplacianMatrix();

            /**
             * Computes the vector L_{w,d}*x where x is stored in @a coordinates and stores the result in @a rhs.
             * @param coordinates The coordinate vector (x in the thesis)
             * @param rhs The right-hand side that stores the result of the matrix-vector multiplication.
             */
            void computeCoordinateLaplacianTerm(const std::vector<NetworKit::Vector>& coordinates, std::vector<NetworKit::Vector>& rhs);

            /**
             * Computes the repulsive forces according to Equation (8) in Gansner et al. exactly (i.e. using no approximation).
             * @param coordinates The current coordinates of the vertices.
             * @param b Repulsive force vector to compute
             */
            std::vector<NetworKit::Vector> computeRepulsiveForces(const std::vector<NetworKit::Vector>& coordinates, std::vector<NetworKit::Vector>& b) const;

            /**
             * Approximates the repulsive forcse by means of an octree (Barnes and Hut).
             * @param coordinates The current coordinates of the vertices.
             * @param octree Octree for Barnes-Hut approximation.
             * @param theta Parameter for Barnes-Hut cell-opening criterion.
             * @param b Repulsive force vector to compute.
             */
            void approxRepulsiveForces(const std::vector<NetworKit::Vector>& coordinates, const NetworKit::Octree<double>& octree, const double theta, std::vector<NetworKit::Vector>& b) const;

            /**
             * Computes the squared distance ||c_i - c_j||^2 between @a coordinates @a i and @a j
             * @param coordinates
             * @param i
             * @param j
             */
            double squaredDistance(const std::vector<NetworKit::Vector>& coordinates, const uint64_t i, const uint64_t j) const;


            /**
             * Computes the distance ||c_i - c_j|| between @a coordinates @a i and @a j.
             * @param coordinates
             * @param i
             * @param j
             */
            inline double distance(const std::vector<NetworKit::Vector>& coordinates, const uint64_t i, const uint64_t j) const {
                return sqrt(squaredDistance(coordinates, i, j));
            }

            /**
             * Computes the squared distance ||c1_i - c2_j||^2 between coordinate @a i from @a coordinates1 and coordinate @a j from @a coordinates 2.
             * @param coordinates1
             * @param coordinates2
             * @param i
             * @param j
             */
            double squaredDistance(const std::vector<NetworKit::Vector>& coordinates1, const std::vector<NetworKit::Vector>& coordinates2, const uint64_t i, const uint64_t j) const;

            /**
             * Computes the squared distance ||c1_i - c2_j|| between coordinate @a i from @a coordinates1 and coordinate @a j from @a coordinates 2.
             * @param coordinates1
             * @param coordinates2
             * @param i
             * @param j
             */
            inline double distance(const std::vector<NetworKit::Vector>& coordinates1, const std::vector<NetworKit::Vector>& coordinates2, const uint64_t i, const uint64_t j) const {
                return sqrt(squaredDistance(coordinates1, coordinates2, i, j));
            }

            /**
             * Computes the squared length ||c_i||^2 of coordinate @a i in @a coordinates.
             * @param coordinates
             * @param i
             */
            double squaredLength(const std::vector<NetworKit::Vector>& coordinates, const uint64_t i) const;

            /**
             * Computes the length ||c_i|| of coordinate @a i in @a coordinates.
             * @param coordinates
             * @param i
             */
            inline double length(const std::vector<NetworKit::Vector>& coordinates, const uint64_t i) const {
                return sqrt(squaredLength(coordinates, i));
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
            inline NetworKit::Point<double> getPoint(const std::vector<NetworKit::Vector>& coordinates, uint64_t i) const {
                NetworKit::Point<double> p(coordinates.size());
                for (uint64_t d = 0; d < p.getDimensions(); ++d) {
                    assert(coordinates[d].getDimension() > i);
                    p[d] = coordinates[d][i];
                }

                return p;
            }
    };

}
#endif /* BIOMAXENTSTRESS_H */
