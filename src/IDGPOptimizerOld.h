/*
 * IDGPOptimizer.h
 *
 *  Created on: Jul 4, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#ifndef NETWORKIT_CPP_VIZ_IDGPOPTIMIZER_H_
#define NETWORKIT_CPP_VIZ_IDGPOPTIMIZER_H_

#include <NetworKit/base/Algorithm.h>
#include <NetworKit/auxiliary/Random.h>
#include <NetworKit/graph/Graph.h>
#include <NetworKit/viz/Point.h>

#include <cmath>
#include <algorithm>
#include <vector>


namespace NetworKit {

/**
 * Stores a distance interval modeled by [lowerBound, upperBound] and optionally also a probability for the interval.
 */
struct DistanceInfo {
	DistanceInfo() : lowerBound(0), upperBound(0), probability(1.0) {}
	DistanceInfo(double lowerBound, double upperBound, double probability = 1.0) : lowerBound(lowerBound), upperBound(upperBound), probability(probability) {}

	double lowerBound;
	double upperBound;
	double probability;
};

/**
 * @ingroup viz
 *
 * Optimizer for interval distance geometry problems.
 * The class incorporates two local optimization algorithms: A simple local optimizer and a simulated annealing local optimizer
 */
class IDGPOptimizer : public Algorithm {
public:
	/**
	 * Determines whether the optimization algorithms should consider distance intervals or not (EXACT vs. BOUNDS) and whether
	 * the optimizers should consider probabilities (BOUNDS_WITH_CONSTRAINTS) for the given distance intervals.
	 */
	enum OptimizationType {
		EXACT,
		BOUNDS,
		BOUNDS_WITH_CONSTRAINTS
	};


	/**
	 * Creates a IDGPOptimizer object for the given Graph @a graph. The coordinates of each vertex are stored in
	 * @a coordinates and are altered by this algorithm (reference). The @a distanceInfo vector stores distance
	 * intervals and optionally the probability for each edge. The @a optimizationType determines how to deal with
	 * the given distanceInfo and @a useSAOptimization determines whether the SA algorithm should be used (in
	 * addition to the simple local optimizer).
	 * @param graph
	 * @param coordinates
	 * @param distanceInfo
	 * @param optimizationType
	 * @param useSAOptimization
	 */
	IDGPOptimizer(const Graph& graph, std::vector<Point<double>>& coordinates, std::vector<DistanceInfo>& distanceInfo, OptimizationType optimizationType, bool useSAOptimization = true) : graph(graph),
																					coordinates(coordinates), distanceInfo(distanceInfo), optimizationType(optimizationType), useSAOptimization(useSAOptimization) {}
	/**
	 * Runs the optimizer algorithm for the graph and the settings given in the constructor.
	 */
	void run() override;

	/**
	 * Computes the largest distance mean error (ldme) for the graph and its coordinates given in the constructor.
	 */
	double ldme() const;


private:
	const Graph& graph;
	std::vector<Point<double>>& coordinates;
	std::vector<DistanceInfo>& distanceInfo;
	OptimizationType optimizationType;

	bool useSAOptimization;

	/**
	 * Edge object that stores incident vertices u and v as well as its weight, the edgeId in the graph
	 * object and the current error.
	 */
	struct Edge {
		node u;
		node v;
		edgeweight weight;
		index edgeId;
		double error;
	};

	/**
	 * Runs the simple local optimizer algorithm.
	 */
	void localOptimization();

	/**
	 * Runs the simulated annealing local optimizer algorithm.
	 */
	void saOptimization();

	/**
	 * Computes the local ldme error of the @a edge where the incident vertices are placed at @a uPos and
	 * @a vPos respectively.
	 * @param edge
	 * @param uPos
	 * @param vPos
	 * @return
	 */
	double computeLocalEdgeError(Edge& edge, Point<double>& uPos, Point<double>& vPos) const;

	/**
	 * Runs a local force optimization for the edge (u,v)
	 * @param u
	 * @param v
	 * @param newU
	 * @param newV
	 */
	void localForceOptimization(node u, node v, Point<double>& newU, Point<double>& newV) const;

	/**
	 * Runs a simple local optimization for the edge (u,v)
	 * @param u
	 * @param v
	 * @param edgeId
	 * @param newU
	 * @param newV
	 */
	void simpleLocalOptimization(node u, node v, index edgeId, Point<double>& newU, Point<double>& newV) const;

	/**
	 * Weighting factor to use for the given @a probability.
	 * @param probability
	 * @return
	 */
	inline double weightingFactor(double probability) const {
		return exp(-5*(1-probability));
	}

	/**
	 * Computes the ldme using the probabilities and the weights to model the penalty function.
	 * @param distInfo
	 * @param dist
	 */
	double computeConstraintError(const DistanceInfo& distInfo, double dist) const;
};

} /* namespace NetworKit */

#endif /* NETWORKIT_CPP_VIZ_IDGPOPTIMIZER_H_ */
