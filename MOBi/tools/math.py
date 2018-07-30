#!/usr/bin/env python3
#####################################
#
# Filename : math.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Tue 23 May 2017 17:01:55 CEST
#
# Last Modified : Mon 30 Jul 2018 03:02:16 PM CEST
#
#####################################

# TODO rename this module

import math as m


def dist_from_angle(angle, d_ij, d_jk):  # AKA law of cosines
    d_ik = m.sqrt(d_ij * d_ij + d_jk * d_jk - 2. * d_ij * d_jk * m.cos(angle))
    return d_ik


def angle_from_dist(a, b, c):
    gamma = m.acos((a * a + b * b - c * c) / (2 * a * b))
    return gamma


# TODO this is probably not optimal. smaller numerical error?
# TODO remove this?
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


def dist_from_dihedral(dihedral, d_ij, d_jk, d_kl, d_ik, d_jl):  # i-l edge for rings and proper-like improper dihedrals

    r_ijk = angle_from_dist(d_ij, d_jk, d_ik)
    r_ljk = angle_from_dist(d_jl, d_jk, d_kl)

    r_ijl = m.acos(m.cos(r_ijk) * m.cos(r_ljk) + m.sin(r_ijk) * m.sin(r_ljk) * m.cos(dihedral))
    d_il = dist_from_angle(r_ijl, d_ij, d_jl)
    return d_il
