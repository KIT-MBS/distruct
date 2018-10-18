
/*************************************
 *
 *  Filename : BioMaxentStress.cpp
 *
 *  Projectname : diSTruct
 *
 *  Author : Michael Wegner
 *
 *  Creation Date : Wed 23 May 2018 05:28:02 PM CEST
 *
 *  Last Modified : Thu 18 Oct 2018 11:25:55 AM CEST
 *
 * *************************************/

// TODO remove and replace this entire thing

#include "BioMaxentStress.h"

namespace diSTruct {

    BioMaxentStress::BioMaxentStress(
            const NetworKit::Graph& G,
            const uint64_t dim,
            const std::vector<NetworKit::Point<double>>& coordinates,
            NetworKit::LinearSolver<NetworKit::CSRMatrix>& solver,
            std::vector<double>& probability,
            bool fastComputation) :
        GraphLayoutAlgorithm<double>(G, dim),
        solver(solver),
        q(0.0),
        alpha(1.0),
        alphaReduction(0.3),
        finalAlpha(0.008),
        convThreshold(0.001*0.001),
        fastComputation(fastComputation),
        maxSolvesPerAlpha(300),
        probability(probability),
        dim(dim),
        hasRun(false)
    {
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
        std::vector<NetworKit::Vector> oldCoordinates(dim, NetworKit::Vector(this->G.upperNodeIdBound()));

        // NOTE init coordinate vectors
#pragma omp parallel for
        for (uint64_t i = 0; i < vertexCoordinates.size(); ++i) {
            for (uint64_t d = 0; d < dim; ++d) {
                oldCoordinates[d][i] = vertexCoordinates[i][d];
            }
        }

        std::vector<NetworKit::Vector> newCoordinates = oldCoordinates;

        Aux::Timer timer; // timer that measures how long the whole algorithm took
        timer.start();

        // Initialization of the currentAlpha <- alpha and the converged parameter
        double currentAlpha = alpha;
        bool converged = false;

        std::vector<NetworKit::Vector> repulsiveForces(dim, NetworKit::Vector(this->G.numberOfNodes(), 0)); // Vector that stores the repulsive forces (entropy term)
        uint64_t currentLowerBound = 0;
        uint64_t newLowerBound = 0;
        while (!converged) { // Run until converged (usually when currentAlpha == finalAlpha)
            INFO("Running with alpha = ", currentAlpha);

            for (uint64_t numSolves = 0; numSolves < maxSolvesPerAlpha; ++numSolves) { // solve up to maxSolvesPerAlpha linear systems
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
                    repulsiveForces = std::vector<NetworKit::Vector>(dim, NetworKit::Vector(G.numberOfNodes(), 0));
                    NetworKit::Octree<double> octree(oldCoordinates); // initialize the octree
                    approxRepulsiveForces(oldCoordinates, octree, 0.6, repulsiveForces); // Barnes-Hut-Approximation using the octree
                    currentLowerBound = newLowerBound;
                }
                t.stop();
                approxTime += t.elapsedMicroseconds();

                t.start();
                std::vector<NetworKit::Vector> rhs(dim, NetworKit::Vector(this->G.numberOfNodes()));
                computeCoordinateLaplacianTerm(oldCoordinates, rhs); // compute L_{w,d}*x and store it in rhs vector
                t.stop();
                rhsTime += t.elapsedMicroseconds();

                t.start();
                for (uint64_t d = 0; d < dim; ++d) {
                    if (numSolves < maxSolvesPerAlpha/5) { // normalize the rhs for the first 20% of solves
                        rhs[d] /= rhs[d].length();
                    }
                    rhs[d] += currentAlpha * repulsiveForces[d]; // add alpha * b(x) to the rhs
                }
                t.stop();
                rhsTime += t.elapsedMicroseconds();

                // correcting rhs to be zero-sum (since it is an approximation)
                NetworKit::Point<double> sum(dim);
                for (uint64_t d = 0; d < dim; ++d) {
                    for (uint64_t i = 0; i < G.numberOfNodes(); ++i) {
                        sum[d] += rhs[d][i];
                    }
                    sum[d] /= G.numberOfNodes();
                }

#pragma omp parallel for
                for (uint64_t i = 0; i < G.numberOfNodes(); ++i) {
                    for (uint64_t d = 0; d < dim; ++d) {
                        rhs[d][i] -= sum[d];
                    }
                }

                t.start();
                // TODO maybe allow more time to solve in later iterations
                solver.parallelSolve(rhs, newCoordinates, maxSolvesPerAlpha*10/3, maxSolvesPerAlpha/3); // solve the Laplacian linear system for each dimension
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
        for (uint64_t i = 0; i < this->G.upperNodeIdBound(); ++i) {
            for (uint64_t d = 0; d < dim; ++d) {
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



    bool BioMaxentStress::isConverged(const std::vector<NetworKit::Vector>& newCoords, const std::vector<NetworKit::Vector>& oldCoords) {
        assert(newCoords.size() == oldCoords.size());

        double relChange = 0.0;
        double oldCoordsSqLength = 0.0;
#pragma omp parallel for reduction(+:relChange,oldCoordsSqLength)
        for (uint64_t i = 0; i < newCoords[0].getDimension(); ++i) {
            relChange += squaredDistance(newCoords, oldCoords, i, i);
            oldCoordsSqLength += squaredLength(oldCoords, i);
        }

        return relChange / oldCoordsSqLength < convThreshold;
    }


    void BioMaxentStress::setupWeightedLaplacianMatrix() {
        uint64_t n = this->G.numberOfNodes();
        std::vector<uint64_t> rowIdx(n+1, 0);
        std::vector<uint64_t> columnIdx(n + G.upperEdgeIdBound()*2, 0); // currently only supports simple graphs
        std::vector<double> nonzeros(columnIdx.size());

        uint64_t idx = 0;
        G.forNodes([&](uint64_t u) {
                double weightedDegree = 0.0;
                G.forNeighborsOf(u, [&](uint64_t u, uint64_t v, double weight, uint64_t edgeId) {
                        //double weightFactor = weightingFactor(weight, edgeId);
                        double weightFactor = probability[edgeId];
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
        for (uint64_t i = 1; i < rowIdx.size(); ++i) {
            rowIdx[i] += rowIdx[i-1];
        }

        NetworKit::CSRMatrix laplacian(n, n, rowIdx, columnIdx, nonzeros);
        solver.setupConnected(laplacian);
    }


    void BioMaxentStress::computeCoordinateLaplacianTerm(const std::vector<NetworKit::Vector>& coordinates, std::vector<NetworKit::Vector>& rhs) {
        G.parallelForNodes([&](uint64_t u) {
                double weightedDegree = 0.0;
                G.forNeighborsOf(u, [&](uint64_t u, uint64_t v, double weight, uint64_t edgeId) {
                        double dist = std::max(distance(coordinates, u, v), 1e-5);
                        double w = probability[edgeId] * weight / dist; // w_{ij} * d_{i,j} / ||x_i - x_j|| NOTE: The last term is multiplied in the paper of Gansner et al. which is wrong!
                        for (uint64_t d = 0; d < dim; ++d) {
                        rhs[d][u] += -w * coordinates[d][v];
                        }
                        weightedDegree += w;
                        });

                for (uint64_t d = 0; d < dim; ++d) {
                rhs[d][u] += weightedDegree * coordinates[d][u];
                }
                });
    }

    std::vector<NetworKit::Vector> BioMaxentStress::computeRepulsiveForces(const std::vector<NetworKit::Vector>& coordinates, std::vector<NetworKit::Vector> &b) const {
        uint64_t n = this->G.numberOfNodes();
        double qSign = sign(this->q);
        double q2 = (q+2)/2;


        G.parallelForNodes([&](uint64_t u) {
                std::vector<bool> knownDist(n, false);
                G.forNeighborsOf(u, [&](uint64_t u, uint64_t v, double weight, uint64_t egdeId) {
                        knownDist[v] = true;
                        });

                G.forNodes([&](uint64_t v) {
                        if (!knownDist[v] && u != v) {
                        double sqDist = std::max(squaredDistance(coordinates, u, v), 1e-3); // ||x_i - x_j||
                        double factor = qSign * 1.0/std::pow(sqDist, q2);
                        for (uint64_t d = 0; d < dim; ++d) {
                        b[d][u] += factor * (coordinates[d][u] - coordinates[d][v]); // sum_{\{i,k\} \in S} dist^{-q-2} * (x_{i,d} - x_{j,d})
                        }
                        }
                        });
                });

        // normalize b
        for (uint64_t d = 0; d < dim; ++d) {
            b[d] /= b[d].length();
        }

        return b;
    }

    void BioMaxentStress::approxRepulsiveForces(const std::vector<NetworKit::Vector>& coordinates, const NetworKit::Octree<double>& octree, const double theta, std::vector<NetworKit::Vector>& b) const {
        double qSign = sign(q);
        double q2 = (q+2)/2;

        G.parallelForNodes([&](uint64_t u) {
                NetworKit::Point<double> posU = getPoint(coordinates, u);
                auto approximateNeighbor = [&](const uint64_t numNodes, const NetworKit::Point<double>& centerOfMass, const double sqDist) {
                if (sqDist < 1e-5) return;
                double factor = qSign * numNodes * 1.0/pow(sqDist, q2);
                for (uint64_t d = 0; d < dim; ++d) {
                b[d][u] += factor * (posU[d] - centerOfMass[d]);
                }
                };

                octree.approximateDistance(posU, theta, approximateNeighbor);

                });

        // normalize b
        for (uint64_t d = 0; d < dim; ++d) {
            b[d] /= b[d].length();
        }
    }

    double BioMaxentStress::squaredDistance(const std::vector<NetworKit::Vector>& coordinates, const uint64_t i, const uint64_t j) const {
        double dist = 0.0;
        for (uint64_t d = 0; d < dim; ++d) {
            double diff = coordinates[d][i] - coordinates[d][j];
            dist += diff * diff;
        }
        return dist;
    }

    double BioMaxentStress::squaredDistance(const std::vector<NetworKit::Vector>& coordinates1, const std::vector<NetworKit::Vector>& coordinates2, const uint64_t i, const uint64_t j) const {
        double dist = 0.0;
        for (uint64_t d = 0; d < dim; ++d) {
            double diff = coordinates1[d][i] - coordinates2[d][j];
            dist += diff * diff;
        }

        return dist;
    }

    double BioMaxentStress::squaredLength(const std::vector<NetworKit::Vector>& coordinates, const uint64_t i) const {
        double length = 0.0;
        for (uint64_t d = 0; d < dim; ++d) {
            length += coordinates[d][i] * coordinates[d][i];
        }

        return length;
    }


} /* namespace NetworKit */
