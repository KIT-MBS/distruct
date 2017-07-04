
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
 *  Last Modified : Tue 04 Jul 2017 11:44:43 CEST
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
            std::vector<double>& weightingFactors,
            const uint64_t dim,
            double alpha
            ) : GraphLayoutAlgorithm<double>(G, dim), solver(solver), weightingFactors(weightingFactors), dim(dim), alpha(alpha), q(0.), convergenceThreshold(0.001*0.001)
    {
        vertexCoordinates = initialCoordinates;
    }

    void BioMaxentStress::run()
    {

        NetworKit::ConnectedComponents cc(this->G);
        cc.run();
        if(cc.numberOfComponents() != 1)
        {
            throw std::invalid_argument("ERROR: Input graph is not connected.");
        }

        if(weightingFactors.size() != G.numberOfEdges())
        {
            throw std::invalid_argument("ERROR: Number of edges and weighting factors does not match");
        }

        setup_weighted_laplacian_matrix(); // = L_{w}

        // NOTE: list of coordinate vectors is 'transposed' to a list of dim lists of coordinates
        // TODO NetworKit::Vector is probably an unnecessarily complex data structure for this. using it for now... 
        // std::vector< std::vector< double > > newCoordinates(dim, std::vector<double>(this->G.upperNodeIdBound()));
        // std::vector< std::vector< double > > oldCoordinates;
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
            // TODO compute diagnostics, do conditional logging of intermediate embeddings
            //TODO stagger entropy calculation

            std::vector<std::vector < double > > entropyTerm(dim, std::vector<double>(this->G.upperNodeIdBound()));

            oldCoordinates = newCoordinates;

            NetworKit::Octree<double> octree(oldCoordinates);
            // TODO clean this up / implement this
            approximate_repulsive_forces();

            std::vector<NetworKit::Vector> rhs(dim, NetworKit::Vector(this->G.numberOfNodes())); // =L_{w,d} * x
            compute_distance_laplacian_term(oldCoordinates, rhs);

            //TODO test how normalization affects solution quality
            //
            
            // TODO correct rhs to be sum zero
            
            uint64_t maxSolverConvergenceTime = 1000;
            uint64_t maxSolverIterations = 75;
            solver.parallelSolve(rhs, newCoordinates, maxSolverConvergenceTime, maxSolverIterations);

            converged = check_converged(newCoordinates, oldCoordinates);

            ++iterationCount;
        }
        // TODO save final result
        // TODO save stats on run
    }


    //void BioMaxentStress::setup_weighted_laplacian_matrix()
    //{
    //    //TODO
    //    return
    //}

    //void BioMaxentStress::approximate_repulsive_forces() const
    //{
    //    //TODO
    //    return

    //}
}
