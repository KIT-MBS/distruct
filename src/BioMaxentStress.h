
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
 *  Last Modified : Fri 30 Jun 2017 16:23:39 CEST
 *
 * *************************************/

#ifndef BIOMAXENTSTRESS_H
#define BIOMAXENTSTRESS_H

#include<cstdint>

#include <NetworKit/viz/GraphLayoutAlgorithm.h>
#include <NetworKit/numerics/LinearSolver.h>
#include <NetworKit/algebraic/CSRMatrix.h>

#include <NetworKit/viz/Octree.h>

namespace MOBi
{

class BioMaxentStress : public NetworKit::GraphLayoutAlgorithm<double>
    {
        public:
            // TODO choose alpha depending on system stats (sum of weighting factors?), add rest of parameters
            // TODO would one ever want fast computation in this context?

            //constructors
            BioMaxentStress(
                    const NetworKit::Graph& G,
                    NetworKit::LinearSolver<NetworKit::CSRMatrix>& solver,
                    std::vector< NetworKit::Point<double> >& initialCoordinates,
                    std::vector<double>& weightingFactors,
                    const uint64_t dim=3,
                    double alpha=1.
                    );
            //destructors
            virtual ~BioMaxentStress() = default;

            virtual void run();

        private:
            NetworKit::LinearSolver<NetworKit::CSRMatrix>& solver;

            std::vector<double> weightingFactors;

            //MaxentStress parameters:

            uint64_t dim;
            //entropy weight
            double alpha;
            //entropy range parameter
            double q;

            double convergenceThreshold;

            //

            void setup_weighted_laplacian_matrix();
            void approximate_repulsive_forces() const;
            void compute_distance_laplacian_term(const std::vector<NetworKit::Vector>, std::vector<NetworKit::Vector> result) const;

            bool check_converged(const std::vector<NetworKit::Vector>, const std::vector<NetworKit::Vector>) const;
            // TODO
            double compute_maxent_stress() const;
            double compute_stress() const;
            double compute_entropy() const;
    };
}

#endif /* BIOMAXENTSTRESS_H */
