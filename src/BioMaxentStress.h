
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
 *  Last Modified : Fri 21 Jul 2017 15:47:25 CEST
 *
 * *************************************/

#ifndef BIOMAXENTSTRESS_H
#define BIOMAXENTSTRESS_H

#include <cstdint>
#include <math.h>

#include <NetworKit/viz/GraphLayoutAlgorithm.h>
#include <NetworKit/numerics/LinearSolver.h>
#include <NetworKit/algebraic/CSRMatrix.h>

#include <NetworKit/viz/Octree.h>

namespace MOBi
{

class BioMaxentStress : public NetworKit::GraphLayoutAlgorithm<double>
    {
        public:
            // TODO choose alpha depending on system stats (sum/average of weighting factors?), add rest of parameters
            // TODO would one ever want fast computation in this context?

            //constructors
            BioMaxentStress(
                    const NetworKit::Graph& G,
                    std::vector< NetworKit::Point<double> >& initialCoordinates,
                    std::vector<double> distances,
                    const uint64_t dim=3,
                    double alpha=1.
                    );
            // TODO a constructor from sequence and edge list
            //destructors
            virtual ~BioMaxentStress() = default;

            virtual void run();

        private:
            NetworKit::LinearSolver<NetworKit::CSRMatrix>& solver;

            std::vector<double> distances;

            //MaxentStress parameters:

            uint64_t dim;
            //entropy weight
            double alpha;
            //entropy shape parameter
            double q;
            // TODO entropy range parameter

            //implementation parameters
            double convergenceThreshold;
            //entropy approximation threshold
            double theta;
            //logging frequency
            uint32_t loggingFrequency;

            //

            void approximate_repulsive_forces(std::vector<NetworKit::Vector> coordinates, NetworKit::Octree<double>& octree, std::vector<NetworKit::Vector>& forces) const;
            void compute_distance_laplacian_term(const std::vector<NetworKit::Vector>, std::vector<NetworKit::Vector> result) const;

            bool check_converged(const std::vector<NetworKit::Vector>&, const std::vector<NetworKit::Vector>&) const;
            // TODO
            double compute_stress() const;
            double compute_entropy() const;
            double compute_largest_distance_mean_error() const;

            void set_vertexCoordinates(std::vector<NetworKit::Vector> coordinates);
            inline double dist2(const std::vector<NetworKit::Vector>&, const std::vector<NetworKit::Vector>&, uint64_t, uint64_t) const;
    };
}

#endif /* BIOMAXENTSTRESS_H */
