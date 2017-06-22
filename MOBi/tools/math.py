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
# Last Modified : Thu 22 Jun 2017 13:06:28 CEST
#
#####################################

# TODO rename this module

import math as m


def dist_from_angle(angle, d_ij, d_jk):  # AKA law of cosines
    print(angle, d_ij, d_jk)
    angle = angle * m.pi / 180.
    d_ik = m.sqrt(d_ij * d_ij + d_jk * d_jk - 2. * d_ij * d_jk * m.cos(angle))
    print(d_ik)
    return d_ik


# TODO this is probably not optimal. smaller numerical error?
def dist_from_improper(dihedral, d_ij, d_jk, d_kl, d_ik, d_jl):
    r_ijk = m.acos((d_ij * d_ij + d_jk * d_jk - d_ik * d_ik) / (2 * d_ij * d_jk))
    r_ljk = m.acos((d_jl * d_jl + d_jk * d_jk - d_kl * d_kl) / (2 * d_jl * d_jk))
    r_ijl = None
    if m.isclose(dihedral, 0.):
        r_ijl = r_ijk - r_ljk
    elif m.isclose(dihedral, 180.):
        r_ijl = r_ijk + r_ljk
    else:
        raise

    d_il = m.sqrt(d_ij * d_ij + d_jl * d_jl - 2. * d_ij * d_jl * m.cos(r_ijl))
    return d_il


# TODO this is needed for secondary structure stuff
def dist_from_dihedral(dihedral, d_ij, d_jk, d_kl, d_ik, d_jl):  # i-l edge for rings and proper-like improper dihedrals
    dihedral = dihedral * m.pi / 180.
    pass
