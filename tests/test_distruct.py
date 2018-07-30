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
# Last Modified : Mon 30 Jul 2018 04:23:09 PM CEST
#
#####################################


def test_NWK():
    import networkit as nwk
    g = nwk.Graph(5)
    g.addEdge(0, 1)
    g.toString()
    return
