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
# Last Modified : Thu 12 Mar 2020 06:19:48 PM CET
#
#####################################


def test_NWK():
    import networkit as nwk
    g = nwk.Graph(5)
    g.addEdge(0, 1)
    return
