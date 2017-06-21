#!/usr/bin/env python3
#####################################
#
# Filename : math.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Tue 23 May 2017 17:01:55 CEST
#
# Last Modified : Wed 21 Jun 2017 14:45:36 CEST
#
#####################################

# TODO rename this
import math


def dist_from_angle(angle, d_ij, d_jk):  # AKA law of cosines
    print(angle, d_ij, d_jk)
    angle = angle * math.pi / 180.
    d_ik = math.sqrt(d_ij * d_ij + d_jk * d_jk - 2. * d_ij * d_jk * math.cos(angle))
    print(d_ik)
    return d_ik


# TODO this is probably not optimal. smaller numerical error?
def dist_from_improper(dihedral, d_ij, d_jk, d_kl, d_ik, d_jl):
    if dihedral not in (0., 180.):
        raise
    # TODO

    pass


# TODO this is needed for secondary structure stuff
def dist_from_dihedral():  # i-l edge for rings and proper-like improper dihedrals
    pass
