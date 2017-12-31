#!/usr/bin/env python3
#####################################
#
# Filename : data.py
#
# Projectname :
#
# Author : Oskar Taubert
#
# Creation Date : Tue 15 Aug 2017 02:15:31 PM CEST
#
# Last Modified : Sun 31 Dec 2017 11:21:03 PM CET
#
#####################################

from Bio import Alphabet
from math import pi
# TODO handle alphabets here
# use generic base stuff and IUPAC
# {'ALA', 'Ala', 'A', 'a'}
# extended

AAs = ['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN', 'ARG', 'HIS', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 'PRO', 'CYS', 'MET']
# RNAs = ['RA', 'RG', 'RC', 'RU']
RNAs = ['A', 'G', 'C', 'U']
DNAs = ['DA', 'DG', 'DC', 'DT']


class ReducedPDBProtein(Alphabet.ThreeLetterProtein):
    letters = AAs
    pass


class PDBRNA(Alphabet.RNAAlphabet):
    letters = RNAs
    pass


PDBReducedProtein = ReducedPDBProtein()
PDBRNAalphabet = PDBRNA()


# TODO put in translation for ff names and other conventions
# NOTE for now the maximally protonated residues are used
# NOTE amber knows HISD, HISE and HISH but not HIS (and calls them HID, HIE and HIP)
amberAAs = ['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN', 'ARG', 'HIP', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 'PRO', 'CYS', 'MET']
amberAAts = ['C' + AA for AA in amberAAs] + ['N' + AA for AA in amberAAs]
# TODO amberNAs

# TODO pdb alphabet
# TODO pdb extended alphabet

alpha_phi = - 60 * pi / 180.
alpha_psi = - 45 * pi / 180.
beta_phi = -135 * pi / 180.
beta_psi = 135 * pi / 180.
hbond = 1.9  # hydrogen - acceptor
hangle = 180. * pi / 180.  # angle donor-hydrogen-acceptor
