
/*************************************
 *
 *  Filename : DuckingWrapper.h
 *
 *  Projectname : diSTruct
 *
 *  Author : Oskar Taubert
 *
 *  Creation Date : Fri 27 Oct 2017 01:41:27 PM CEST
 *
 *  Last Modified : Mon 30 Jul 2018 03:03:35 PM CEST
 *
 * *************************************/

#ifndef DUCKINGWRAPPER_H
#define DUCKINGWRAPPER_H

#include <vector>

#include <cpp/graph/Graph.h>
#include <cpp/numerics/LAMG/Lamg.h>
#include <cpp/algebraic/CSRMatrix.h>

#include <cpp/components/ConnectedComponents.h>

#include "BioMaxentStress.h"

// TODO remove this and use the proper networkit cython interface

std::vector<NetworKit::Point<double>> runMaxent(uint64_t numNodes, double alpha, double q, uint64_t solves, std::vector<std::pair<uint64_t, uint64_t>> edges, std::vector<double> distances, std::vector<double> probabilities/*, std::vector<NetworKit::Point<double>> initCoords*/)
{
    NetworKit::Graph graph(numNodes, true, false);

    std::cout << "BioMaxentStress init: " << std::endl;

    uint64_t numEdges = edges.size();


    for(uint64_t i=0; i<numEdges; ++i)
    {
        if(graph.hasEdge(edges[i].first, edges[i].second))
        {
            throw std::invalid_argument("ERROR: Input graph is not simple.");
        }
        graph.addEdge(edges[i].first, edges[i].second, distances[i]);
    }
    std::cout << "vertices: " << graph.numberOfNodes() << " edges: " << graph.numberOfEdges() << std::endl;
    graph.indexEdges();

    std::vector<double> sortedProbabilities(numEdges);
    for(uint64_t i=0; i<numEdges; ++i)
    {
        uint64_t edgeID = graph.edgeId(edges[i].first, edges[i].second);
        sortedProbabilities[edgeID] = probabilities[i];
    }
    probabilities = sortedProbabilities;

    NetworKit::ConnectedComponents cc(graph);
    cc.run();
    if(cc.numberOfComponents() != 1)
    {
        std::cout << "Number of connected components: " << cc.numberOfComponents() << std::endl;
        std::cout << "Sizes of connected components: " << std::endl;
        //std::map<uint64_t, uint64_t> ccsSizes = cc.getComponentSizes();
        //for(auto iter : ccsSizes)
        //{
        //    // NOTE dereferencing a map iterator yields a pair
        //    std::cout << iter.first << " " << iter.second << std::endl;
        //}
        std::vector<std::vector<uint64_t>> concom = cc.getComponents();
        int i = 0;
        for(auto iter : concom)
        {
            std::cout << "component: " << ++i << " size: " << iter.size() << " root: " << iter[0] << std::endl;
        }

        throw std::invalid_argument("ERROR: Input graph is not connected.");
    }

    NetworKit::Lamg<NetworKit::CSRMatrix> lamg(1e-5);

    NetworKit::BioMaxentStress maxent(graph, 3, lamg, probabilities, false);
    //if(initCoords.size() != numNodes)
    //{
    //    NetworKit::BioMaxentStress maxent(graph, 3, lamg, probabilities, false);
    //}
    //else
    //{
    //    NetworKit::BioMaxentStress maxent(graph, 3, initCoords, lamg, probabilities, false);
    //}
    //
    // TODO do all the checks

    // defaults
    alpha = 1.;
    double alphaReduction = .3;
    double finalAlpha = .008;
    solves = 300;


    maxent.setAlpha(alpha);
    maxent.setAlphaReduction(alphaReduction);
    maxent.setFinalAlpha(finalAlpha);
    maxent.setMaxSolvesPerAlpha(solves);
    maxent.setConvergenceThreshold(1e-5);
    maxent.setQ(q);

    maxent.run();

    return maxent.getCoordinates();
}

#endif /* DUCKINGWRAPPER_H */