#!/usr/bin/env python3
#####################################
#
# Filename : test_distruct.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Tue 09 May 2017 13:53:41 CEST
#
# Last Modified : Thu 02 Aug 2018 01:39:09 PM CEST
#
#####################################


def test_NWK():
    import networkit as nwk
    g = nwk.Graph(5)
    g.addEdge(0, 1)
    g.toString()
    return


def test_protein_xray():
    assert False
    return


def test_protein_nmr():
    assert False
    return


def test_rna_xray():
    assert False
    return


def test_rna_nmr():
    assert False
    return
