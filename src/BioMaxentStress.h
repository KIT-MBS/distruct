
/*************************************
 *
 *  Filename : BioMaxentStress.h
 *
 *  Projectname : MOBi
 *
 *  Author : Michael Wegner, Oskar Taubert
 *
 *  Creation Date : Tue 27 Jun 2017 11:53:39 CEST
 *
 *  Last Modified : Wed 04 Oct 2017 03:34:18 PM CEST
 *
 * *************************************/

#ifndef BIOMAXENTSTRESS_H
#define BIOMAXENTSTRESS_H

#include <cstdint>
#include <math.h>

#include <NetworKit/viz/GraphLayoutAlgorithm.h>
//#include <NetworKit/numerics/LinearSolver.h>
#include <NetworKit/numerics/LAMG/Lamg.h>
#include <NetworKit/algebraic/CSRMatrix.h>

#include <NetworKit/viz/Octree.h>

//provides keyboard interrupt handling
#include "aux/sig.h"

namespace MOBi
{

class BioMaxentStress /*: public NetworKit::GraphLayoutAlgorithm<double>*/
    {
        public:
            // TODO choose alpha depending on system stats (sum/average of weighting factors?), add rest of parameters
            //constructors
            /*
            BioMaxentStress(
                    const NetworKit::Graph& G,
                    std::vector< NetworKit::Point<double> >& initialCoordinates,
                    std::vector<double> distances,
                    const uint64_t dim=3,
                    double alpha=1.
                    );
                    */

            BioMaxentStress(
                    uint64_t numNodes,
                    std::vector<std::pair<uint64_t, uint64_t>>& edges,
                    std::vector<double>& weights,
                    std::vector<double>& distances,
                    std::vector<NetworKit::Point<double>>& initialCoordinates,
                    double alpha, // TODO default args with cython -_-
                    double q,
                    uint32_t loggingFrequency
                    );
            //destructors
            virtual ~BioMaxentStress() = default;

            void run(uint64_t maxSolves = 100);

            std::vector<NetworKit::Point<double>> getCoordinates() const
            {
                return vertexCoordinates;
            }

        private:
            // TODO if there will ever be different solvers
            // NetworKit::LinearSolver<NetworKit::CSRMatrix>& solver;
            NetworKit::Lamg<NetworKit::CSRMatrix> solver;

            std::vector<double> distances;

            //MaxentStress parameters:

            uint64_t dim;
            //entropy weight / temperataure
            double alpha;
            //entropy shape parameter
            double q;
            // TODO entropy range parameter
            // TODO atom radius!!!

            //implementation parameters
            // TODO should depend on number of atoms!
            double convergenceThreshold;
            //entropy approximation threshold
            double theta;
            //logging frequency
            uint32_t loggingFrequency;
            // TODO track number of runs?

            bool converged;

            NetworKit::Graph G;
            std::vector<NetworKit::Point<double>> vertexCoordinates;
            //

            void approximate_repulsive_forces(std::vector<NetworKit::Vector> coordinates, NetworKit::Octree<double>& octree, std::vector<NetworKit::Vector>& forces) const;
            void compute_distance_laplacian_term(const std::vector<NetworKit::Vector>, std::vector<NetworKit::Vector>& result) const;

            bool check_converged(const std::vector<NetworKit::Vector>&, const std::vector<NetworKit::Vector>&) const;
            double compute_stress() const;
            double compute_entropy() const;
            double compute_distance_error() const;

            void set_vertexCoordinates(std::vector<NetworKit::Vector> coordinates);
            inline double dist2(const std::vector<NetworKit::Vector>&, const std::vector<NetworKit::Vector>&, uint64_t, uint64_t) const;
    };
}

#endif /* BIOMAXENTSTRESS_H */
