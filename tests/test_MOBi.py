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
# Last Modified : Wed 20 Jun 2018 01:50:21 PM CEST
#
#####################################


def test_NWK():
    import networkit as nwk
    g = nwk.Graph(5)
    g.addEdge(0, 1)
    g.toString()
    return
