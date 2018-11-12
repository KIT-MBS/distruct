diSTruct
========

python package name: distruct
Implementation of Maxent-stress Optimization of Biomolecular systems.

*diSTruct* is in essence an implementation of the MaxEnt-Stress graph drawing algorithm (Gansner, E.; Hu, Y. and North, S. C.: "A Maxent-Stress Model for Graph Layout" in *IEE Trans. Vis. Comput. Graph.* 2013) for generating biomolecular structures from distance constraints.

## Requirements
The actual structure generation in diSTruct is build on the MaxEnt-stress graph drawing implementation in [NetworKit](https://networkit.iti.kit.edu).

It aims to conveniently combine this with the PDB module in [Biopython](https://biopython.org).

It also depends on *lxml* and requires *cython*.

## Installation Instructions

First install NetworKit following the instructions in their [homepage](https://networkit.iti.kit.edu).
Make sure you can
import networkit
without error. Sometimes there are issues with automatically installing all dependenicies.
Simply install missing packages manually.

Install cython, Biopython and lxml.
Install *diSTruct* with

pip install distruct

To run *diSTruct* it may be necessary to set your LD_LIBRARY_PATH to find the NetworKit extension.
*diSTruct* gives further instructions on importing.

## Publications

## Tests
To run the python tests run "make test".
