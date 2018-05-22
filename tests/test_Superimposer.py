#!/usr/bin/env python3
#####################################
#
# Filename : test_Superimposer.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Fri 18 May 2018 06:28:53 PM CEST
#
# Last Modified : Fri 18 May 2018 08:10:15 PM CEST
#
#####################################


import numpy as np


import MOBi


def test_superimposer():
    x = np.array([[51.65, -1.90, 50.07],
         [50.40, -1.23, 50.65],
         [50.68, -0.04, 51.54],
         [50.22, -0.02, 52.85]], 'f')

    y = np.array([[51.30, -2.99, 46.54],
         [51.09, -1.88, 47.58],
         [52.36, -1.20, 48.03],
         [52.71, -1.18, 49.38]], 'f')
    # x = np.array([[-1., 0., 5.],
    #              [1., 0., 5.]], 'f')
    # y = np.array([[0., 1., 0.],
    #              [0., -1., 0.]], 'f')

    sup = MOBi.Superimposer()
    sup.set_coords(x, y)
    from pytest import approx
    # TODO is this really that bad??
    assert sup.rms == approx(0., abs=1e-2)

