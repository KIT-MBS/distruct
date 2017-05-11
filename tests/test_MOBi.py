#!/usr/bin/env python3
#####################################
#
# Filename : MOBi_tests.py
#
# Projectname :
#
# Author : Oskar Taubert
#
# Creation Date : Tue 09 May 2017 13:53:41 CEST
#
# Last Modified : Tue 09 May 2017 17:07:19 CEST
#
#####################################

# import MOBi


def test_NWK():
    try:
        import networkit as nwk
        g = nwk.Graph(5)
        g.addEdge(0, 1)
        g.toString()
    except:
        assert 0
