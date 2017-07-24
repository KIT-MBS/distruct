
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
 *  Last Modified : Mon 17 Jul 2017 18:34:49 CEST
 *
 * *************************************/

#include "BioMaxentStress.h"

#include <NetworKit/components/ConnectedComponents.h>

namespace MOBi
{
    BioMaxentStress::BioMaxentStress(
            const NetworKit::Graph& G,
            NetworKit::LinearSolver<NetworKit::CSRMatrix>& solver,
            std::vector< NetworKit::Point<double> >& initialCoordinates,
            std::vector<double> distances,
            const uint64_t dim,
            double alpha
            ) : GraphLayoutAlgorithm<double>(G, dim), solver(solver), distances(distances), dim(dim), alpha(alpha), q(0.), convergenceThreshold(0.001*0.001), theta(.6), loggingFrequency(10)
    {
        // TODO check if input valid here (without asserts)
        // NOTE distances have to be indexed correctly!!!
        vertexCoordinates = initialCoordinates;
        assert(vertexCoordinates.size() == G.numberOfNodes());
    }

    void BioMaxentStress::run()
    {

        NetworKit::ConnectedComponents cc(this->G);
        cc.run();
        if(cc.numberOfComponents() != 1)
        {
            throw std::invalid_argument("ERROR: Input graph is not connected.");
        }

        if(distances.size() != G.numberOfEdges())
        {
            throw std::invalid_argument("ERROR: Number of edges and distances does not match.");
        }

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

        bool converged = false;

        uint64_t iterationCount = 0;
        while(!converged)
        {
            //TODO
            {
                if(loggingFrequency != 0 && iterationCount % loggingFrequency == 0)
                {
                    // TODO compute stress terms, approximate entropy, intermediate embedding and save them
                    // TODO add output facilities
                    // TODO add loggingFrequency input facilities
                }
            }
            // TODO stagger entropy calculation

            oldCoordinates = newCoordinates;

            std::vector<NetworKit::Vector> repulsiveForces(dim, std::vector<double>(this->G.numberOfNodes()));
            NetworKit::Octree<double> octree(oldCoordinates);
            approximate_repulsive_forces(oldCoordinates, octree, repulsiveForces);

            std::vector<NetworKit::Vector> rightHandSides(dim, NetworKit::Vector(this->G.numberOfNodes())); // =L_{w,d} * x
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
        }

        set_vertexCoordinates(newCoordinates);
        // TODO save stats on run
    }


    void BioMaxentStress::approximate_repulsive_forces(std::vector<NetworKit::Vector> coordinates, NetworKit::Octree<double>& octree, std::vector<NetworKit::Vector>& result) const
    {
        // TODO implement different qs
        // TODO this also takes neighbours into acccount?
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
        return distSum / (oldLengthSum) < convergenceThreshold;
    }

    double BioMaxentStress::compute_stress() const
    {
        // TODO
        return 0.;
    }


    double BioMaxentStress::compute_entropy() const
    {
        // TODO
        return 0.;
    }


    double BioMaxentStress::compute_largest_distance_mean_error() const
    {
        // TODO
        return 0.;
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