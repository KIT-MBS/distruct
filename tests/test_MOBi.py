#!/usr/bin/env python3
#####################################
#
# Filename : MOBi_tests.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Tue 09 May 2017 13:53:41 CEST
#
# Last Modified : Fri 01 Jun 2018 02:46:54 PM CEST
#
#####################################


def test_NWK():
    import networkit as nwk
    g = nwk.Graph(5)
    g.addEdge(0, 1)
    g.toString()
    return

# TODO test other dependencies

# def test_GROMACS():
#     # TODO
#     return


def test_Distructure():
    from MOBi import Distructure
    from MOBi import data

    code = "1ptq"
    fileName = code + '.pdb'

    # generate distruct
    # create structure
    ds = Distructure(code, sequences)
    # create edges
    # create primary edges
    ds.generate_primary_edges()
    # create tertiary edges
    # get contact map
    contactMap = MOBi.tools.generate_atomic_contacts(refStructure)
    ds.generate_tertiary_edges(contactMap)
    # generate coordinates
    ds.run()

    # TODO structural alignment
    RMSD = MOBi.tools.align(ds, refStructure)

    return
