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
# Last Modified : Thu 17 Aug 2017 02:24:50 PM CEST
#
#####################################

from Bio import Alphabet
# TODO handle alphabets here
# use generic base stuff and IUPAC
# {'ALA', 'Ala', 'A', 'a'}
# extended

AAs = ['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN', 'ARG', 'HIS', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 'PRO', 'CYS', 'MET']
RNAs = ['RA', 'RG', 'RC', 'RU']
DNAs = ['DA', 'DG', 'DC', 'DT']


class ReducedPDBProtein(Alphabet.ThreeLetterProtein):
    letters = AAs
    pass

PDBReducedProtein = ReducedPDBProtein()

# TODO put in translation for ff names and other conventions
# NOTE for now the maximally protonated residues are used
# NOTE amber knows HISD, HISE and HISH but not HIS (and calls them HID, HIE and HIP)
amberAAs = ['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN', 'ARG', 'HIP', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 'PRO', 'CYS', 'MET']
amberAAts = ['C' + AA for AA in amberAAs] + ['N' + AA for AA in amberAAs]
# TODO amberNAs

# TODO pdb alphabet
# TODO pdb extended alphabet
