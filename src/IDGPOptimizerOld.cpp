/*
 * IDGPOptimizer.cpp
 *
 *  Created on: Jul 4, 2016
 *      Author: Michael Wegner (michael.wegner@student.kit.edu)
 */

#include "IDGPOptimizer.h"
#include <NetworKit/algebraic/Vector.h>

namespace NetworKit {

void IDGPOptimizer::run() {
	if (!graph.hasEdgeIds()) throw std::runtime_error("IDGPOptimizer::run(): The graph does not have edge ids. Index edges first.");

	double startLdme = ldme();
	INFO("ldme = ", startLdme);

	if (useSAOptimization) {
		std::vector<Point<double>> temp = coordinates;

		INFO("Running SA optimization");
		saOptimization();

		double saLdme = ldme();
		INFO("Bounds ldme after SA: ", saLdme);

		if (saLdme > startLdme) {
			coordinates = temp;
		}
	}



	INFO("Running local optimization");
	localOptimization();

	INFO("Final ldme: ", ldme());
}

double IDGPOptimizer::ldme() const {
	double ldme = 0.0;
	count numViolatingEdges = 0;
	graph.forEdges([&](node u, node v, double val, index edgeId) {
		double dist = coordinates[u].distance(coordinates[v]);
		double error = std::max(distanceInfo[edgeId].lowerBound - dist, std::max(dist - distanceInfo[edgeId].upperBound, 0.0));
		if (error >= 1e-9) numViolatingEdges++;
		ldme += error * error;
	});
	ldme /= graph.upperEdgeIdBound();
	INFO("numViolatingEdges = ", numViolatingEdges);
	return sqrt(ldme);
}

void IDGPOptimizer::localOptimization() {
	std::vector<bool> lock(graph.upperEdgeIdBound(), false);
	bool improvement = true;

	count numRounds = 0;
	while (improvement && numRounds < 50) {
		lock = std::vector<bool>(graph.upperEdgeIdBound(), false);
		improvement = false;

		std::vector<Edge> violatingEdges;
		graph.forEdges([&](node u, node v, double weight, index edgeId) {
			double dist = coordinates[u].distance(coordinates[v]);
			double error = 0.0;

			switch (optimizationType) {
			case BOUNDS:
				error = std::max(distanceInfo[edgeId].lowerBound - dist, std::max(dist - distanceInfo[edgeId].upperBound, 0.0));
				break;
			case BOUNDS_WITH_CONSTRAINTS:
				error = computeConstraintError(distanceInfo[edgeId], dist);
				break;
			default:
				error = fabs(weight - dist);
				break;
			}

			if (error > 1e-9) {
				violatingEdges.push_back({u,v, weight, edgeId, error * error});
			}
		});


		if (optimizationType == BOUNDS_WITH_CONSTRAINTS) { // sort by descending certaintyFactor and within the same certaintyFactor by descending error
			std::sort(violatingEdges.begin(), violatingEdges.end(), [&](const Edge& e1, const Edge& e2) {
				return distanceInfo[e1.edgeId].probability >= distanceInfo[e2.edgeId].probability && e1.error > e2.error;
			});
		} else {// sort by descending error
			std::sort(violatingEdges.begin(), violatingEdges.end(), [&](const Edge& e1, const Edge& e2) {return e1.error > e2.error;});
		}

		for (Edge& edge : violatingEdges) {
			index edgeId = edge.edgeId;
			if (lock[edgeId]) continue;
			node u = edge.u;
			node v = edge.v;
			double dist = coordinates[u].distance(coordinates[v]);
			double error = 0.0;
			switch (optimizationType) {
			case BOUNDS:
				error = std::max(distanceInfo[edgeId].lowerBound - dist, std::max(dist - distanceInfo[edgeId].upperBound, 0.0));
				break;
			case BOUNDS_WITH_CONSTRAINTS:
				error = computeConstraintError(distanceInfo[edgeId], dist);
				break;
			default:
				error = fabs(edge.weight - dist);
				break;
			}

			if (error <= 1e-9) continue;

			// compute ldme impact of this edge and its neighbors
			double localError = computeLocalEdgeError(edge, coordinates[u], coordinates[v]);

			// translate the nodes u and v in the right direction
			Point<double> newU = coordinates[u];
			Point<double> newV = coordinates[v];
			simpleLocalOptimization(u, v, edgeId, newU, newV);

			double newError = computeLocalEdgeError(edge, newU, newV);
			if (localError > newError) { // improvement
				coordinates[u] = newU;
				coordinates[v] = newV;
				improvement = true;

				lock[edgeId] = true;

				graph.forNeighborsOf(u, [&](node u, node w, double uwWeight, index uwEdgeId) {
					if (w != v) {
						lock[uwEdgeId] = true;
					}
				});

				graph.forNeighborsOf(v, [&](node v, node w, double vwWeight, index vwEdgeId) {
					if (w != u) {
						lock[vwEdgeId] = true;
					}
				});
			}
		}

		numRounds++;
	}
}

void IDGPOptimizer::saOptimization() {
	count numStepsWithNoImprovement = 0;
	count numModifications = 0;
	count numRounds = 0;

	std::srand(Aux::Random::getSeed());
	std::vector<Edge> edges(graph.upperEdgeIdBound());
	graph.forEdges([&](node u, node v, edgeweight weight, index edgeId) {
		double dist = coordinates[u].distance(coordinates[v]);
		double error = 0.0;
		switch (optimizationType) {
		case BOUNDS:
			error = std::max(distanceInfo[edgeId].lowerBound - dist, std::max(dist - distanceInfo[edgeId].upperBound, 0.0));
			break;
		case BOUNDS_WITH_CONSTRAINTS:
			error = computeConstraintError(distanceInfo[edgeId], dist);
			break;
		default:
			error = fabs(weight - dist);
			break;
		}
		edges[edgeId] = {u,v,weight,edgeId, error};
	});
	count m = graph.numberOfEdges();

	double t = 0.3;
	while (numStepsWithNoImprovement < m && t > 1e-7) {
		INFO("Running with t = ", t, " numRoundsWithNoImprovement = ", numStepsWithNoImprovement);
		bool change = true;
		while (numRounds < 2*m && numModifications < 0.5*m && change) {
			change = false;
#pragma omp parallel for reduction(+:numStepsWithNoImprovement, numRounds, numModifications)
			for (index i = 0; i < edges.size(); ++i) {
				Edge& edge = edges[i];

				double localError = computeLocalEdgeError(edge, coordinates[edge.u], coordinates[edge.v]);

				Point<double> newU, newV;
				localForceOptimization(edge.u, edge.v, newU, newV);
				double newLocalError = computeLocalEdgeError(edge, newU, newV);

				bool acceptModification = newLocalError < localError;
				if (!acceptModification) {
					double boltzman = std::exp((localError - newLocalError)/t);
					double randNum = Aux::Random::real();
					acceptModification =  randNum < boltzman;
				}

				if (acceptModification) {
					coordinates[edge.u] = newU;
					coordinates[edge.v] = newV;
					edge.error = newLocalError;
					numModifications++;
					change = true;
				} else {
					numStepsWithNoImprovement++;
				}

				numRounds++;
			}
		}

		numRounds = 0;
		numModifications = 0;
		t *= 0.1;
	}
}

double IDGPOptimizer::computeLocalEdgeError(Edge& edge, Point<double>& uPos, Point<double>& vPos) const {
	double dist = uPos.distance(vPos);
	double error = 0.0;
	switch (optimizationType) {
	case BOUNDS:
		error = std::max(distanceInfo[edge.edgeId].lowerBound - dist, std::max(dist - distanceInfo[edge.edgeId].upperBound, 0.0));
		break;
	case BOUNDS_WITH_CONSTRAINTS:
		error = computeConstraintError(distanceInfo[edge.edgeId], dist);
		break;
	default:
		error = fabs(edge.weight - dist);
		break;
	}

	double localError = error * error;
	graph.forNeighborsOf(edge.u, [&](node u, node w, double uwWeight, index uwEdgeId) {
		if (w != edge.v) {
			double dist = uPos.distance(coordinates[w]);
			double uwError = 0.0;
			switch (optimizationType) {
			case BOUNDS:
				uwError = std::max(distanceInfo[uwEdgeId].lowerBound - dist, std::max(dist - distanceInfo[uwEdgeId].upperBound, 0.0));
				break;
			case BOUNDS_WITH_CONSTRAINTS:
				uwError = computeConstraintError(distanceInfo[uwEdgeId], dist);
				break;
			default:
				uwError = fabs(uwWeight - dist);
				break;
			}
			localError += uwError * uwError;
		}
	});

	graph.forNeighborsOf(edge.v, [&](node v, node w, double vwWeight, index vwEdgeId) {
		if (w != edge. u) {
			double dist = vPos.distance(coordinates[w]);
			double vwError = 0.0;
			switch (optimizationType) {
			case BOUNDS:
				vwError = std::max(distanceInfo[vwEdgeId].lowerBound - dist, std::max(dist - distanceInfo[vwEdgeId].upperBound, 0.0));
				break;
			case BOUNDS_WITH_CONSTRAINTS:
				vwError = computeConstraintError(distanceInfo[vwEdgeId], dist);
				break;
			default:
				vwError = fabs(vwWeight - dist);
				break;
			}
			localError += vwError * vwError;
		}
	});

	return localError;
}

void IDGPOptimizer::localForceOptimization(node u, node v, Point<double>& newU, Point<double>& newV) const {
	newU = coordinates[u];
	newV = coordinates[v];
	count iter = 0;
	double alpha = 0.01;
	while (iter < 20) {
		Point<double> fU(coordinates[u].getDimensions());
		Point<double> fV(coordinates[v].getDimensions());

		graph.forNeighborsOf(u, [&](node u, node w, double uwWeight, index uwEdgeId) {
			if (w != v) {
				double dist = std::max(newU.distance(coordinates[w]), 1e-9);
				if (dist < distanceInfo[uwEdgeId].lowerBound) { // distance smaller than lower bound => add repulsive force
					double lowerBound = distanceInfo[uwEdgeId].lowerBound;
					fU += (newU - coordinates[w]) * (lowerBound * lowerBound / (dist * dist));
				} else if (dist > distanceInfo[uwEdgeId].upperBound) { // distance larger than upper bound => add attractive force
					double upperBound = distanceInfo[uwEdgeId].upperBound;
					fU += (coordinates[w] - newU) * (upperBound * upperBound / (dist * dist));
				}

			} else {
				double dist = std::max(newU.distance(newV), 1e-9);
				if (dist < distanceInfo[uwEdgeId].lowerBound) { // distance smaller than lower bound => add repulsive force
					double lowerBound = distanceInfo[uwEdgeId].lowerBound;
					fU += (newU - newV) * (lowerBound * lowerBound / (dist * dist));
				} else if (dist > distanceInfo[uwEdgeId].upperBound) { // distance larger than upper bound => add attractive force
					double upperBound = distanceInfo[uwEdgeId].upperBound;
					fU += (newV - newU) * (upperBound * upperBound / (dist * dist));
				}
			}
		});

		graph.forNeighborsOf(v, [&](node v, node w, double vwWeight, index vwEdgeId) {
			if (w != u) {
				double dist = std::max(newV.distance(coordinates[w]), 1e-9);
				if (dist < distanceInfo[vwEdgeId].lowerBound) { // distance smaller than lower bound => add repulsive force
					double lowerBound = distanceInfo[vwEdgeId].lowerBound;
					fV += (newV - coordinates[w]) * (lowerBound * lowerBound / (dist * dist));
				} else if (dist > distanceInfo[vwEdgeId].upperBound) { // distance larger than upper bound => add attractive force
					double upperBound = distanceInfo[vwEdgeId].upperBound;
					fV += (coordinates[w] - newV) * (upperBound * upperBound / (dist * dist));
				}
			} else {
				double dist = std::max(newV.distance(newU), 1e-9);
				if (dist < distanceInfo[vwEdgeId].lowerBound) { // distance smaller than lower bound => add repulsive force
					double lowerBound = distanceInfo[vwEdgeId].lowerBound;
					fV += (newV - newU) * (lowerBound * lowerBound / (dist * dist));
				} else if (dist > distanceInfo[vwEdgeId].upperBound) { // distance larger than upper bound => add attractive force
					double upperBound = distanceInfo[vwEdgeId].upperBound;
					fV += (newU - newV) * (upperBound * upperBound / (dist * dist));
				}
			}
		});

		newU += fU * alpha;
		newV += fV * alpha;

		alpha *= 0.7;
		iter++;
	}
}

void IDGPOptimizer::simpleLocalOptimization(node u, node v, index edgeId, Point<double>& newU, Point<double>& newV) const {
	newU = coordinates[u];
	newV = coordinates[v];
	Point<double> vu = newU - newV;
	double dist = coordinates[u].distance(coordinates[v]);
	if (dist > distanceInfo[edgeId].upperBound) {
		double s = (1 - distanceInfo[edgeId].upperBound / newU.distance(newV)) / 2;
		newU -= vu * s;
		newV += vu * s;
	} else {
		double s = (distanceInfo[edgeId].lowerBound / newU.distance(newV) - 1) / 2;
		newU += vu * s;
		newV -= vu * s;
	}
}

double IDGPOptimizer::computeConstraintError(const DistanceInfo& distInfo, double dist) const {
	return std::max(distInfo.lowerBound - dist, std::max(dist - distInfo.upperBound, 0.0)) * ( 1+5*weightingFactor(distInfo.probability));
}

} /* namespace NetworKit */
