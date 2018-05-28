#!/usr/bin/env python3
#####################################
#
# Filename : data.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Tue 15 Aug 2017 02:15:31 PM CEST
#
# Last Modified : Mon 28 May 2018 10:14:55 AM CEST
#
#####################################

from Bio import Alphabet
from math import pi
# TODO handle alphabets here
# use generic base stuff and IUPAC
# {'ALA', 'Ala', 'A', 'a'}
# extended

AAs = ['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN', 'ARG', 'HIS', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 'PRO', 'CYS', 'MET']
# RNAs = ['RA', 'RG', 'RC', 'RU']  # names as they appear in gromacs' amber99sb-ildn
RNAs = ['A', 'G', 'C', 'U']  # names as they appear in the PDB
DNAs = ['DA', 'DG', 'DC', 'DT']


class ReducedPDBProtein(Alphabet.ThreeLetterProtein):
    letters = AAs
    pass


class PDBRNA(Alphabet.RNAAlphabet):
    letters = RNAs
    pass


PDBReducedProtein = ReducedPDBProtein()
PDBRNAalphabet = PDBRNA()


# TODO handle terminals (C,N,3', and 5')
# TODO put in translation for ff names and other conventions
# NOTE for now the maximally protonated residues are used
# NOTE amber knows HISD, HISE and HISH but not HIS (and calls them HID, HIE and HIP)
# amberAAs = ['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN', 'ARG', 'HIP', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 'PRO', 'CYS', 'MET']
# amberAAts = ['C' + AA for AA in amberAAs] + ['N' + AA for AA in amberAAs]
# TODO amberNAs

# TODO pdb alphabet
# TODO pdb extended alphabet

# NOTE default values to determine edges for secondary structure
# TODO improve and cite these
alpha_phi = - 60 * pi / 180.
alpha_psi = - 45 * pi / 180.
beta_phi = -135 * pi / 180.
beta_psi = 135 * pi / 180.
hbond = 1.9  # hydrogen - acceptor
hangle = 180. * pi / 180.  # angle donor-hydrogen-acceptor

# TODO default omega value

# TODO maybe need nucleic acid specific stuff

# TODO put in alphabet with AAs and NAs
from .fileio import read_topology_file
defaultTopologyDB = read_topology_file("amber99sb-ildn_protein")
