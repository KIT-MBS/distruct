/*
 * BioMaxentStress.cpp
 *
 *  Created on: 06.07.2016
 *      Author: Michael Wegner
 */

#include "BioMaxentStressOldOld.h"

//#include "SplitTree.h"

#include <NetworKit/auxiliary/Log.h>
#include <NetworKit/auxiliary/PrioQueue.h>

//#include <NetworKit/numerics/LAMG/Lamg.h>
#include <NetworKit/numerics/LinearSolver.h>

#include <queue>

namespace NetworKit {

BioMaxentStress::BioMaxentStress(const Graph& G, const count dim, LinearSolver<CSRMatrix>& solver, bool fastComputation) : GraphLayoutAlgorithm<double>(G, dim), solver(solver), q(0.0), alpha(1.0), alphaReduction(0.3), finalAlpha(0.008), convThreshold(0.001*0.001), coordinatesProvided(false), probabilityProvided(false), fastComputation(fastComputation), maxSolvesPerAlpha(50), dim(dim), hasRun(false) {
    signal(SIGINT, handle_signals);
}

BioMaxentStress::BioMaxentStress(const Graph& G, const count dim, LinearSolver<CSRMatrix>& solver, std::vector<double>& probability, bool fastComputation) : GraphLayoutAlgorithm<double>(G, dim), solver(solver), q(0.0), alpha(1.0), alphaReduction(0.3), finalAlpha(0.008), convThreshold(0.001*0.001), coordinatesProvided(false), probabilityProvided(true), fastComputation(fastComputation), maxSolvesPerAlpha(50), probability(probability), dim(dim), hasRun(false) {
    signal(SIGINT, handle_signals);
}

BioMaxentStress::BioMaxentStress(const Graph& G, const count dim, const std::vector<Point<double>>& coordinates, LinearSolver<CSRMatrix>& solver, bool fastComputation) : GraphLayoutAlgorithm<double>(G, dim), solver(solver), q(0.0), alpha(1.0), alphaReduction(0.3), finalAlpha(0.008), convThreshold(0.001*0.001), coordinatesProvided(true), probabilityProvided(false), fastComputation(fastComputation), maxSolvesPerAlpha(50), dim(dim), hasRun(false) {
	vertexCoordinates = coordinates;
    signal(SIGINT, handle_signals);
}

BioMaxentStress::BioMaxentStress(const Graph& G, const count dim, const std::vector<Point<double>>& coordinates, LinearSolver<CSRMatrix>& solver, std::vector<double>& probability, bool fastComputation) : GraphLayoutAlgorithm<double>(G, dim), solver(solver), q(0.0), alpha(1.0), alphaReduction(0.3), finalAlpha(0.008), convThreshold(0.001*0.001), coordinatesProvided(true), probabilityProvided(true), fastComputation(fastComputation), maxSolvesPerAlpha(50), probability(probability), dim(dim), hasRun(false) {
	vertexCoordinates = coordinates;
    signal(SIGINT, handle_signals);
}


void BioMaxentStress::run() {
	// Initialization of timers for performance benchmarks
	Aux::Timer t;
	double solveTime = 0;
	double rhsTime = 0;
	double approxTime = 0.0;

    /////////////////////////
    std::cout << "############################" << std::endl;
    std::cout << "setup laplacian and solver" << std::endl;
    /////////////////////////

	t.start();
	setupWeightedLaplacianMatrix(); // create weighted Laplacian matrix and setup the solver
	t.stop();
	solveTime += t.elapsedMicroseconds();


    /////////////////////////
    std::cout << "initializing coordinates" << std::endl;
    /////////////////////////
	CoordinateVector oldCoordinates(dim, Vector(this->G.upperNodeIdBound()));

	if (!coordinatesProvided && !hasRun) { // no coordinates have been provided => use coordinates from random sphere placement algorithm
		//randomSphereCoordinates(oldCoordinates);
        randomInitCoordinates(oldCoordinates);
	} else { // use the initial coordinates provided in the constructor
#pragma omp parallel for
		for (index i = 0; i < vertexCoordinates.size(); ++i) {
			for (index d = 0; d < dim; ++d) {
				oldCoordinates[d][i] = vertexCoordinates[i][d];
			}
		}
	}

	CoordinateVector newCoordinates = oldCoordinates;

	Aux::Timer timer; // timer that measures how long the whole algorithm took
	timer.start();

	// Initialization of the currentAlpha <- alpha and the converged parameter
	double currentAlpha = alpha;
	bool converged = false;

    /////////////////////////
    //std::cout << "edges:" << std::endl;
    //G.forEdges([&](uint64_t u, uint64_t v, double weight, uint64_t ){
    //        //std::cout << u << " " << v << " " << weight << " " << probability[edgeID] << std::endl;
    //        std::cout << u << " " << v << " " << weight << " " << std::endl;
    //        });
    /////////////////////////
    
	CoordinateVector repulsiveForces(dim, Vector(this->G.numberOfNodes(), 0)); // Vector that stores the repulsive forces (entropy term)
	count currentLowerBound = 0;
	count newLowerBound = 0;
	while (!converged) { // Run until converged (usually when currentAlpha == finalAlpha)
		INFO("Running with alpha = ", currentAlpha);

		for (count numSolves = 0; numSolves < maxSolvesPerAlpha; ++numSolves) { // solve up to maxSolvesPerAlpha linear systems
			oldCoordinates = newCoordinates;

            for(uint64_t d=0; d<dim; ++d)
            {
                for(uint64_t i=0; i<G.numberOfNodes(); ++i)
                {
                    if(std::isnan(oldCoordinates[d][i]))  // this needs C++11
                    {
                        throw std::range_error("ERROR: NaN encountered");
                    }
                }
            }
            /////////////////////////
            //std::cout << "numSolves: "<< numSolves << std::endl;
            //std::cout << currentAlpha << std::endl;
            //for(uint64_t i=0; i<G.numberOfNodes(); ++i)
            //{
            //    for(uint64_t d=0; d<dim; ++d)
            //    {
            //        std::cout << oldCoordinates[d][i] << " ";
            //    }
            //    std::cout << std::endl;
            //}
            /////////////////////////
            
			t.start();
			newLowerBound = floor(5 * std::log(numSolves));
			if (newLowerBound != currentLowerBound) { // Lazy approximation of entropy terms, if bounds are different we trigger a recomputation
				repulsiveForces = CoordinateVector(dim, Vector(G.numberOfNodes(), 0));
				Octree<double> octree(oldCoordinates); // initialize the octree
				approxRepulsiveForces(oldCoordinates, octree, 0.6, repulsiveForces); // Barnes-Hut-Approximation using the octree
				currentLowerBound = newLowerBound;
			}
			t.stop();
			approxTime += t.elapsedMicroseconds();

			t.start();
			CoordinateVector rhs(dim, Vector(this->G.numberOfNodes()));
			computeCoordinateLaplacianTerm(oldCoordinates, rhs); // compute L_{w,d}*x and store it in rhs vector
			t.stop();
			rhsTime += t.elapsedMicroseconds();

			t.start();
			for (index d = 0; d < dim; ++d) {
				if (numSolves < maxSolvesPerAlpha/5) { // normalize the rhs for the first 20% of solves
					rhs[d] /= rhs[d].length();
				}
				rhs[d] += currentAlpha * repulsiveForces[d]; // add alpha * b(x) to the rhs
			}
			t.stop();
			rhsTime += t.elapsedMicroseconds();

			// correcting rhs to be zero-sum (since it is an approximation)
			Point<double> sum(dim);
			for (index d = 0; d < dim; ++d) {
				for (index i = 0; i < G.numberOfNodes(); ++i) {
					sum[d] += rhs[d][i];
				}
				sum[d] /= G.numberOfNodes();
			}

#pragma omp parallel for
			for (index i = 0; i < G.numberOfNodes(); ++i) {
				for (index d = 0; d < dim; ++d) {
					rhs[d][i] -= sum[d];
				}
			}

			t.start();
			solver.parallelSolve(rhs, newCoordinates, 1000, 75); // solve the Laplacian linear system for each dimension
			t.stop();
			solveTime += t.elapsedMicroseconds();

			converged = isConverged(newCoordinates, oldCoordinates); // checks if the algorithm converged in this iteration
			if (converged) {
				INFO("Converged after ", numSolves+1, " solves");
				if (!fastComputation) { // if fastComputation is turned off (false) then we continue with the next smaller alpha
					converged = false;
				}
				break;
			}
		}

		currentAlpha *= alphaReduction; // Cooling: reduce alpha for next round

		converged = converged || currentAlpha < finalAlpha;
	}

	// write coordinates to vertexCoordinates vector
#pragma omp parallel for
	for (index i = 0; i < this->G.upperNodeIdBound(); ++i) {
		for (index d = 0; d < dim; ++d) {
			this->vertexCoordinates[i][d] = newCoordinates[d][i];
		}
	}
	timer.stop();

	hasRun = true;

	// output runtime statistics
	INFO("Time spent on rhs ", rhsTime/1000, " and solve time was ", solveTime/1000);
	INFO("approxTime: ", approxTime/1000);
	INFO("Graph drawn in ", timer.elapsedMilliseconds());
    /////////////////////////
    std::cout << "############################" << std::endl;
    /////////////////////////
}



bool BioMaxentStress::isConverged(const CoordinateVector& newCoords, const CoordinateVector& oldCoords) {
	assert(newCoords.size() == oldCoords.size());

	double relChange = 0.0;
	double oldCoordsSqLength = 0.0;
#pragma omp parallel for reduction(+:relChange,oldCoordsSqLength)
	for (index i = 0; i < newCoords[0].getDimension(); ++i) {
		relChange += squaredDistance(newCoords, oldCoords, i, i);
		oldCoordsSqLength += squaredLength(oldCoords, i);
	}

	return relChange / oldCoordsSqLength < convThreshold;
}


void BioMaxentStress::setupWeightedLaplacianMatrix() {
	count n = this->G.numberOfNodes();
	std::vector<index> rowIdx(n+1, 0);
	std::vector<index> columnIdx(n + G.upperEdgeIdBound()*2, 0); // currently only supports simple graphs
	std::vector<double> nonzeros(columnIdx.size());

	index idx = 0;
	G.forNodes([&](node u) {
		double weightedDegree = 0.0;
		G.forNeighborsOf(u, [&](node u, node v, edgeweight weight, index edgeId) {
			double weightFactor = weightingFactor(weight, edgeId);
			columnIdx[idx] = v;
			nonzeros[idx] = -weightFactor;
			weightedDegree += weightFactor;
			idx++;
		});

		// add diagonal element
		columnIdx[idx] = u;
		nonzeros[idx] = weightedDegree;
		idx++;

		rowIdx[u+1] = G.degree(u)+1; // +1 for diagonal
	});

	// compute correct rowIdx offsets
	for (index i = 1; i < rowIdx.size(); ++i) {
		rowIdx[i] += rowIdx[i-1];
	}

	CSRMatrix laplacian(n, n, rowIdx, columnIdx, nonzeros);
	solver.setupConnected(laplacian);
}


void BioMaxentStress::computeCoordinateLaplacianTerm(const CoordinateVector& coordinates, CoordinateVector& rhs) {
	G.parallelForNodes([&](node u) {
		double weightedDegree = 0.0;
		G.forNeighborsOf(u, [&](node u, node v, edgeweight weight, index edgeId) {
			double dist = std::max(distance(coordinates, u, v), 1e-5);
			double w = weightingFactor(weight, edgeId) * weight / dist; // w_{ij} * d_{i,j} / ||x_i - x_j|| NOTE: The last term is multiplied in the paper of Gansner et al. which is wrong!
			for (index d = 0; d < dim; ++d) {
				rhs[d][u] += -w * coordinates[d][v];
			}
			weightedDegree += w;
		});

		for (index d = 0; d < dim; ++d) {
			rhs[d][u] += weightedDegree * coordinates[d][u];
		}
	});
}

CoordinateVector BioMaxentStress::computeRepulsiveForces(const CoordinateVector& coordinates, CoordinateVector &b) const {
	count n = this->G.numberOfNodes();
	double qSign = sign(this->q);
	double q2 = (q+2)/2;


	G.parallelForNodes([&](node u) {
		std::vector<bool> knownDist(n, false);
		G.forNeighborsOf(u, [&](node u, node v, edgeweight weight, index egdeId) {
			knownDist[v] = true;
		});

		G.forNodes([&](node v) {
			if (!knownDist[v] && u != v) {
				double sqDist = std::max(squaredDistance(coordinates, u, v), 1e-3); // ||x_i - x_j||
				double factor = qSign * 1.0/std::pow(sqDist, q2);
				for (index d = 0; d < dim; ++d) {
					b[d][u] += factor * (coordinates[d][u] - coordinates[d][v]); // sum_{\{i,k\} \in S} dist^{-q-2} * (x_{i,d} - x_{j,d})
				}
			}
		});
	});

	// normalize b
	for (index d = 0; d < dim; ++d) {
		b[d] /= b[d].length();
	}

	return b;
}

void BioMaxentStress::approxRepulsiveForces(const CoordinateVector& coordinates, const Octree<double>& octree, const double theta, CoordinateVector& b) const {
	double qSign = sign(q);
	double q2 = (q+2)/2;

	G.parallelForNodes([&](node u) {
		Point<double> posU = getPoint(coordinates, u);
		auto approximateNeighbor = [&](const count numNodes, const Point<double>& centerOfMass, const double sqDist) {
			if (sqDist < 1e-5) return;
			double factor = qSign * numNodes * 1.0/pow(sqDist, q2);
			for (index d = 0; d < dim; ++d) {
				b[d][u] += factor * (posU[d] - centerOfMass[d]);
			}
		};

		octree.approximateDistance(posU, theta, approximateNeighbor);

	});

	// normalize b
	for (index d = 0; d < dim; ++d) {
		b[d] /= b[d].length();
	}
}

void BioMaxentStress::randomInitCoordinates(CoordinateVector &coordinates) const {

#pragma omp parallel for
	for (index i = 0; i < coordinates[0].getDimension(); ++i) {
		for (index d = 0; d < dim; ++d) {
			coordinates[d][i] = Aux::Random::real() * 50; // 100 x 100 pixel
		}
	}
}

void BioMaxentStress::randomSphereCoordinates(CoordinateVector &coordinates) const {
	// find node with highest degree
	node maxDegNode = 0;
	count maxDeg = G.degree(maxDegNode);
	G.forNodes([&](node u) {
		if (G.degree(u) > maxDeg) {
			maxDegNode = u;
			maxDeg = G.degree(u);
		}
	});


	// set coordinate of node 0 to (0,0,...,0)
	for (index d = 0; d < dim; ++d) {
		coordinates[d][maxDegNode] = 0.0;
	}

	std::vector<bool> coordinateSet(G.upperNodeIdBound(), false);
	coordinateSet[maxDegNode] = true;
	count numSet = 1;

	while (numSet < G.numberOfNodes()) {
		node start = 0;
		G.forNodes([&](node u) {
			if (coordinateSet[u]) {
				start = u;
				return;
			}
		});

		// perform BFS from start
		std::queue<node> Q;
		Q.push(start);
		while (!Q.empty()) {
			node u = Q.front(); Q.pop();
			G.forNeighborsOf(u, [&](node u, node v, edgeweight w, index) {
				if (!coordinateSet[v]) {
					Vector p(dim);
					for (index d = 0; d < dim; ++d) {
						p[d] = 2*Aux::Random::real() - 1;
					}
					p *= w / p.length();
					for (index d = 0; d < dim; ++d) {
						coordinates[d][v] = coordinates[d][u] + p[d];
					}
					coordinateSet[v] = true;
					numSet++;
					Q.push(v);
				}
			});
		}
	}
}

double BioMaxentStress::squaredDistance(const CoordinateVector& coordinates, const index i, const index j) const {
	double dist = 0.0;
	for (index d = 0; d < dim; ++d) {
		double diff = coordinates[d][i] - coordinates[d][j];
		dist += diff * diff;
	}
	return dist;
}

double BioMaxentStress::squaredDistance(const CoordinateVector& coordinates1, const CoordinateVector& coordinates2, const index i, const index j) const {
	double dist = 0.0;
	for (index d = 0; d < dim; ++d) {
		double diff = coordinates1[d][i] - coordinates2[d][j];
		dist += diff * diff;
	}

	return dist;
}

double BioMaxentStress::squaredLength(const CoordinateVector& coordinates, const index i) const {
	double length = 0.0;
	for (index d = 0; d < dim; ++d) {
		length += coordinates[d][i] * coordinates[d][i];
	}

	return length;
}



} /* namespace NetworKit */
