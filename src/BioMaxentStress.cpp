
/*************************************
 *
 *  Filename : BioMaxentStress.cpp
 *
 *  Projectname : MOBi
 *
 *  Author : Michael Wegner, Oskar Taubert
 *
 *  Creation Date : Tue 27 Jun 2017 11:53:49 CEST
 *
 *  Last Modified : Tue 22 Aug 2017 03:19:12 PM CEST
 *
 * *************************************/

#include "BioMaxentStress.h"

#include <NetworKit/components/ConnectedComponents.h>


namespace MOBi
{
    /*
    BioMaxentStress::BioMaxentStress(
            const NetworKit::Graph& G,
            std::vector< NetworKit::Point<double> >& initialCoordinates,
            std::vector<double> distances,
            const uint64_t dim,
            double alpha
            ) : GraphLayoutAlgorithm<double>(G, dim), solver(NetworKit::Lamg<NetworKit::CSRMatrix>(1e-5)), distances(distances), dim(dim), alpha(alpha), q(0.), convergenceThreshold(0.001*0.001), theta(.6), loggingFrequency(10)
    {
        // TODO check if input valid here (without asserts)
        // NOTE distances have to be indexed correctly!!!
        vertexCoordinates = initialCoordinates;
        assert(vertexCoordinates.size() == G.numberOfNodes());
    }
    */

    BioMaxentStress::BioMaxentStress(
            uint64_t numNodes,
            std::vector<std::pair<uint64_t, uint64_t>>& edges,
            std::vector<double>& weights,
            std::vector<double>& distances,
            std::vector<NetworKit::Point<double>>& initialCoordinates
            ) : dim(3), alpha(1.), q(0.), convergenceThreshold(0.001*0.001), theta(.6), loggingFrequency(0)
    {
        // TODO maybe there are faster ways to generate graphs from edges?
        std::cout << "BioMaxentStress init: " << std::endl;
        G = NetworKit::Graph(numNodes, true, false);

        // TODO replace with exceptions, check valid graph
        uint64_t numEdges = edges.size();
        if(numEdges != weights.size())
        {
            std::cout << numEdges << " " << weights.size() << std::endl;
            throw std::invalid_argument("ERROR: Number of edges and weights does not match.");
        }

        for(uint64_t i=0; i<numEdges; ++i)
        {
            G.addEdge(edges[i].first, edges[i].second, weights[i]);
        }
        std::cout << "indexing..." << std::endl;
        G.indexEdges();

        NetworKit::ConnectedComponents cc(this->G);
        cc.run();
        if(cc.numberOfComponents() != 1)
        {
            std::cout << "Number of connected components: " << cc.numberOfComponents() << std::endl;
            std::cout << "Sizes of connected components: " << std::endl;
            std::map<uint64_t, uint64_t> ccsSizes = cc.getComponentSizes();
            for(auto iter : ccsSizes)
            {
                // NOTE dereferencing a map iterator yields a pair
                std::cout << iter.first << " " << iter.second << std::endl;
            }
            
            throw std::invalid_argument("ERROR: Input graph is not connected.");
        }

        if(distances.size() != G.numberOfEdges())
        {
            std::cout << distances.size() << " " << G.numberOfEdges() << std::endl;
            throw std::invalid_argument("ERROR: Number of edges and distances does not match.");
        }

        std::cout << "resorting edges..." << std::endl;
        std::vector<std::pair<uint64_t, uint64_t>> sortedEdges (numEdges, std::pair<uint64_t, uint64_t> (0, 0));
        std::vector<double> sortedWeights (numEdges, 0.);
        std::vector<double> sortedDistances (numEdges, 0.);
        for(uint64_t i=0; i<numEdges; ++i)
        {
            uint64_t edgeID = G.edgeId(edges[i].first, edges[i].second);
            sortedEdges[edgeID] = edges[i];
            sortedWeights[edgeID] = weights[i];
            sortedDistances[edgeID] = distances[i];
        }

        edges = sortedEdges;
        weights = sortedWeights;
        // TODO better naming
        BioMaxentStress::distances = sortedDistances;

        vertexCoordinates = initialCoordinates;

        //TODO this is not nice, replace with nwk handler?
        signal(SIGINT, handle_signals);
    }

    void BioMaxentStress::run(uint64_t maxSolves)
    {

        NetworKit::CSRMatrix weighted_laplacian = NetworKit::CSRMatrix::laplacianMatrix(G); // L_{w}
        solver.setupConnected(weighted_laplacian);

        // TODO use other vector data structure??
        // NOTE: list of coordinate vectors is 'transposed' to a list of dim lists of coordinates
        std::vector<NetworKit::Vector> newCoordinates(dim, NetworKit::Vector(G.numberOfNodes()));
        std::vector<NetworKit::Vector> oldCoordinates;

#pragma omp parallel for
        for(uint64_t i = 0; i < vertexCoordinates.size(); ++i)
        {
            for(uint64_t j = 0; j < dim; j++)
            {
                newCoordinates[j][i] = vertexCoordinates[i][j];
            }
        }

        std::cout << "BioMaxentStress running with:" << std::endl;
        std::cout << "alpha = " << alpha << std::endl;
        std::cout << "theta = " << theta << std::endl;

        converged = false;

        uint64_t iterationCount = 0;

        uint64_t currentStaggerThreshold = 0;
        uint64_t newStaggerThreshold = 0;

        std::vector<NetworKit::Vector> repulsiveForces(dim, std::vector<double>(this->G.numberOfNodes()));
        std::vector<NetworKit::Vector> rightHandSides(dim, NetworKit::Vector(this->G.numberOfNodes())); // =L_{w,d} * x
        while(!converged)
        {
            // TODO remove
            std::cout << "=============================" << std::endl;
            std::cout << "iteration: " << iterationCount << std::endl;
            std::cout << "alpha = " << alpha << std::endl;
            {
                if(loggingFrequency != 0 && iterationCount % loggingFrequency == 0)
                {
                    set_vertexCoordinates(newCoordinates);
                    // TODO compute intermediate embedding and save it
                    double stress = compute_stress();
                    double entropy = compute_entropy();
                    double ldme = compute_largest_distance_mean_error();

                    std::cout << "stress: " << stress << " entropy: " << entropy << " ldme: " << ldme <<std::endl;
                    // TODO add average relative error
                    // TODO add loggingFrequency input facilities
                    // TODO add logging
                    // TODO add units in output
                }
            }
            // TODO stagger entropy calculation (probably only useful for large number of iterations)

            oldCoordinates = newCoordinates;

            newStaggerThreshold = floor(5 * std::log(iterationCount));
            if(newStaggerThreshold != currentStaggerThreshold)
            {
                NetworKit::Octree<double> octree(oldCoordinates);
                approximate_repulsive_forces(oldCoordinates, octree, repulsiveForces); // TODO pass theta
                currentStaggerThreshold = newStaggerThreshold;
            }

            compute_distance_laplacian_term(oldCoordinates, rightHandSides);

            //TODO test how normalization affects solution quality

            for(uint64_t d=0; d<dim; ++d)
            {
                // add entropy term
                rightHandSides[d] += alpha * repulsiveForces[d];
            }

            //correct rightHandSides to be sum zero
            {
                std::vector<double> sums(dim, 0.);
#pragma omp parallel for
                for(uint64_t d=0; d<dim; ++d)
                {
                    for(uint64_t i=0; i<G.numberOfNodes(); ++i)
                    {
                        sums[d] += rightHandSides[d][i];
                    }
                    sums[d] /= G.numberOfNodes();
                }
#pragma omp parallel for
                for(uint64_t i=0; i<G.numberOfNodes(); ++i)
                {
                    for(uint64_t d=0; d<dim; ++d)
                    {
                        rightHandSides[d][i] -= sums[d];
                    }
                }
            }

            // TODO check these values
            uint64_t maxSolverConvergenceTime = 1000;
            uint64_t maxSolverIterations = 75;
            solver.parallelSolve(rightHandSides, newCoordinates, maxSolverConvergenceTime, maxSolverIterations);

            //TODO check on the ldme as well? (and max steps)
            converged = check_converged(newCoordinates, oldCoordinates);

            ++iterationCount;
            if(converged)
            {
                std::cout << "converged after " << iterationCount << " solves" << std::endl;
            }
            if(iterationCount >= maxSolves)
            {
                alpha *= .3;
                // TODO
                if(alpha < 0.008)
                {
                    converged = true;
                    std::cout << "converged (not really) after " << iterationCount << " solves" << std::endl;
                }
                else
                {
                    run(maxSolves);
                }
            }
        }

        set_vertexCoordinates(newCoordinates);
        // TODO save stats on run
        // TODO center
        // TODO invert if laevus
    }


    void BioMaxentStress::approximate_repulsive_forces(std::vector<NetworKit::Vector> coordinates, NetworKit::Octree<double>& octree, std::vector<NetworKit::Vector>& result) const
    {
        // TODO implement different qs
        if(fabs(q) > 0.001)
        {
            std::cout << "currently only q=0 supported" << std::endl;
            throw;
        }
        // TODO should take vdw radius into account
        G.parallelForNodes([&](uint64_t u) {

            NetworKit::Point<double> uCoordinates(dim);
            for(uint64_t d=0; d<dim; ++d)
            {
                uCoordinates[d] = coordinates[d][u];
            }

            auto approximate_force_entry = [&](const uint64_t weight, const NetworKit::Point<double> centerOfMass, double dist2)
            {
                if (dist2 < 1e-5) return;
                for(uint64_t d=0; d<dim; ++d)
                {
                    result[d][u] = weight * (uCoordinates[d] - centerOfMass[d]) / dist2;
                }
                return;
            };

            octree.approximateDistance(uCoordinates, theta, approximate_force_entry);
        });
    }


    // alternativly create matrix L_{w, d} and multiply by coordinate vector for every dimension
    void BioMaxentStress::compute_distance_laplacian_term(std::vector<NetworKit::Vector> coordinates, std::vector<NetworKit::Vector> result) const
    {
        // NOTE result assumed to be zero initialized
        G.parallelForNodes([&](uint64_t u) {
            double weightedDegree = 0.;
            G.forNeighborsOf(u, [&](uint64_t u, uint64_t v, double weight, uint64_t edgeID) {
                double targetDistance = distances[edgeID];
                double currentDistance = sqrt(dist2(coordinates, coordinates, u, v));
                for(uint64_t d=0; d<dim; ++d)
                {
                    result[d][u] += - weight * targetDistance / currentDistance * coordinates[d][v];
                    weightedDegree += weight * targetDistance / currentDistance;
                }
            });
            for(uint64_t d=0; d<dim; ++d)
            {
                result[d][u] += weightedDegree * coordinates[d][u];
            }
        });
        return;
    }


    bool BioMaxentStress::check_converged(const std::vector<NetworKit::Vector>& newCoordinates, const std::vector<NetworKit::Vector>& oldCoordinates) const
    {
        double distSum = 0.;
        double oldLengthSum = 0.;
#pragma omp parallel for reduction(+:distSum,oldLengthSum)
        for(uint64_t u=0; u<G.numberOfNodes(); ++u)
        {
            distSum += dist2(newCoordinates, oldCoordinates, u, u);
            for(uint64_t d=0; d<dim; ++d)
            {
                oldLengthSum += oldCoordinates[d][u] * oldCoordinates[d][u];
            }
        }
        std::cout << "distSum: " << distSum  << " oldLengthSum: " << oldLengthSum << std::endl;
        std::cout << "convergence criterion: " << distSum / oldLengthSum << std::endl;
        return distSum / (oldLengthSum) < convergenceThreshold;
    }

    // TODO rename this to internal energy?
    double BioMaxentStress::compute_stress() const
    {
        double result = 0.;

#pragma omp parallel for reduction(+:result)
        for(uint64_t u=0; u<G.numberOfNodes(); ++u)
        {
            G.forNeighborsOf(u, [&](uint64_t v, double w)
                    {
                        double currentDistance = std::max(vertexCoordinates[u].distance(vertexCoordinates[v]), 1e-5);
                        uint64_t edgeID = G.edgeId(u, v);

                        result += w * (distances[edgeID] - currentDistance) * (distances[edgeID] - currentDistance);
                    });
        }
        result /= 2;

        return result;
    }


    double BioMaxentStress::compute_entropy() const
    {
        double result = 0.;
        // NOTE this is not the entropy mentioned in the theory. this is summed over all pairs of atoms
#pragma omp parallel for reduction(+:result)
        for(uint64_t u=0; u<G.numberOfNodes(); ++u)
        {
            G.forNodes([&](uint64_t v)
                    {
                        if(u!=v)
                        {
                            double currentDistance = std::max(vertexCoordinates[u].distance(vertexCoordinates[v]), 1e-5);
                            result += (fabs(q) < 0.001) ? log(currentDistance) : pow(currentDistance, -q);
                        }
                    });
        }
        return alpha * result;
    }


    double BioMaxentStress::compute_largest_distance_mean_error() const
    {
        double result = 0.;

        G.forEdges([&](uint64_t u, uint64_t v, uint64_t edgeID){
                double currentDistance = vertexCoordinates[u].distance(vertexCoordinates[v]);
                result += (distances[edgeID] - currentDistance) * (distances[edgeID] - currentDistance);
                });
        result /= distances.size();
        return sqrt(result);
    }


    void BioMaxentStress::set_vertexCoordinates(std::vector<NetworKit::Vector> coordinates)
    {
        //NOTE this accepts both the list of geometric coordinate vectors and a list of algebraic vectors by dimension as they are used in maxent stress.
        if(dim == G.numberOfNodes())
        {
            // TODO warn here
        }

        if(coordinates.size() == dim && coordinates[0].getDimension() == G.numberOfNodes())
        {
#pragma omp parallel
            for(uint64_t i=0; i<G.numberOfNodes(); ++i)
            {
                for(uint64_t d=0; d<dim; ++d)
                {
                    vertexCoordinates[i][d] = coordinates[d][i];
                }
            }
        }
        else if(coordinates.size() == G.numberOfNodes() && coordinates[0].getDimension() == dim)
        {
#pragma omp parallel
            for(uint64_t i=0; i<G.numberOfNodes(); ++i)
            {
                for(uint64_t d=0; d<dim; ++d)
                {
                    vertexCoordinates[i][d] = coordinates[i][d];
                }
            }
        }
        else
        {
            throw std::invalid_argument("ERROR: Tried to set invalid set of coordinates.");
        }
    }


    inline double BioMaxentStress::dist2(const std::vector<NetworKit::Vector>& coordinateSet1, const std::vector<NetworKit::Vector>& coordinateSet2, uint64_t u, uint64_t v) const
    {
        assert(coordinateSet1.size() == coordinateSet2.size());
        double result = 0.;
        for(uint64_t d = 0; d<dim; ++d)
        {
            double diff = (coordinateSet1[d][u] - coordinateSet2[d][v]);
            result += diff * diff;
        }
        return result;
    }
}
