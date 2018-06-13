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
# Last Modified : Wed 13 Jun 2018 06:07:22 PM CEST
#
#####################################

from Bio.Alphabet import IUPAC
from Bio import Alphabet
from Bio.Data import IUPACData
from math import pi
# TODO handle alphabets here
# use generic base stuff and IUPAC
# {'ALA', 'Ala', 'A', 'a'}
# extended

# PDBAAs = ['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN', 'ARG', 'HIS', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 'PRO', 'CYS', 'MET']
# # RNAs = ['RA', 'RG', 'RC', 'RU']  # names as they appear in gromacs' amber99sb-ildn
# PDBRNAs = ['A', 'G', 'C', 'U']  # names as they appear in the PDB
# PDBDNAs = ['DA', 'DG', 'DC', 'DT']

# NOTE Alphabets as they occur in PDB entries. Things like extended alphabets are not considered for now.
# NOTE Ambiguities are not useful in this context. Extensions might be.
AALetters = [IUPACData.protein_letters_1to3[letter].upper() for letter in IUPACData.protein_letters]
DNALetters = ['D' + letter for letter in IUPACData.unambiguous_dna_letters]
RNALetters = [letter for letter in IUPACData.unambiguous_rna_letters]

class AAAlphabet(Alphabet.ThreeLetterProtein):
    letters = AALetters
    pass

class DNAAlphabet(Alphabet.DNAAlphabet):
    letters = DNALetters
    size = 2
    pass

class RNAAlphabet(Alphabet.RNAAlphabet):
    letters = RNALetters
    pass

class MoleculeAlphabet(Alphabet.Alphabet()):
    letters = AALetters + DNALetters + RNALetters
    pass

AAalphabet = AAAlphabet()
DNAalphabet = DNAAlphabet()
RNAalphabet = RNAAlphabet()
alphabet = MoleculeAlphabet()

# TODO secondary structure

# TODO handle terminals (C,N,3', and 5')
# TODO put in translation for ff names and other conventions
# NOTE for now the maximally protonated residues are used
# NOTE amber knows HISD, HISE and HISH but not HIS (and calls them HID, HIE and HIP)
# amberAAs = ['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN', 'ARG', 'HIP', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 'PRO', 'CYS', 'MET']
# amberAAts = ['C' + AA for AA in amberAAs] + ['N' + AA for AA in amberAAs]

# NOTE default values to determine edges for secondary structure
# TODO improve and cite these
alpha_phi = - 60 * pi / 180.
alpha_psi = - 45 * pi / 180.
beta_phi = -135 * pi / 180.
beta_psi = 135 * pi / 180.
hbond = 1.9  # hydrogen - acceptor
hangle = 180. * pi / 180.  # angle donor-hydrogen-acceptor

# TODO default omega value

from .fileio import read_topology_file
defaultTopologyDB = read_topology_file("amber99sb-ildn")
