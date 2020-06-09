diSTruct
========

python package name: distruct
Implementation of Maxent-stress Optimization of Biomolecular systems.

*diSTruct* is in essence an implementation of the MaxEnt-Stress graph drawing algorithm (Gansner, E.; Hu, Y. and North, S. C.: "A Maxent-Stress Model for Graph Layout" in *IEE Trans. Vis. Comput. Graph.* 2013) for generating biomolecular structures from distance constraints.

## Requirements
The actual structure generation in diSTruct is build on the MaxEnt-stress graph drawing implementation in [NetworKit](https://networkit.iti.kit.edu).

It aims to conveniently combine this with the PDB module in [Biopython](https://biopython.org).

It also depends on *lxml* and requires *Cython*.

## Installation Instructions

First install NetworKit following the instructions in their [homepage](https://networkit.iti.kit.edu).
Make sure you can
import networkit
without error. Sometimes there are issues with automatically installing all dependenicies.
Simply install missing packages manually.

Install Cython, Biopython and lxml.
Install *diSTruct* with

pip install distruct

## Publications
Please cite 

    @article{10.1093/bioinformatics/btz578,
        author = {Taubert, Oskar and Reinartz, Ines and Meyerhenke, Henning and Schug, Alexander},
        title = "{diSTruct v1.0: generating biomolecular structures from distance constraints}",
        journal = {Bioinformatics},
        year = {2019},
        month = {07},
        abstract = "{The distance geometry problem is often encountered in molecular biology and the life sciences at large, as a host of experimental methods produce ambiguous and noisy distance data. In this note, we present diSTruct; an adaptation of the generic MaxEnt-Stress graph drawing algorithm to the domain of biological macromolecules. diSTruct is fast, provides reliable structural models even from incomplete or noisy distance data and integrates access to graph analysis tools.diSTruct is written in C++, Cython and Python 3. It is available from https://github.com/KIT-MBS/distruct.git or in the Python package index under the MIT license.Supplementary data are available at Bioinformatics online.}",
        issn = {1367-4803},
        doi = {10.1093/bioinformatics/btz578},
        url = {https://doi.org/10.1093/bioinformatics/btz578},
        note = {btz578},
        eprint = {http://oup.prod.sis.lan/bioinformatics/advance-article-pdf/doi/10.1093/bioinformatics/btz578/29161648/btz578.pdf},
    }

## Tests
To run the python tests install pytest an run python -m pytest in /path/to/distruct/tests
