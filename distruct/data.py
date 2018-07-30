#!/usr/bin/env python3
#####################################
#
# Filename : data.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Tue 15 Aug 2017 02:15:31 PM CEST
#
# Last Modified : Mon 30 Jul 2018 03:03:05 PM CEST
#
#####################################

import os


from Bio.Alphabet import IUPAC
from Bio import Alphabet
from Bio.Data import IUPACData
from math import pi

# TODO remove these. have to work with just the default bio alphabets
# NOTE Alphabets as they occur in PDB entries.
# NOTE Things like extended alphabets are not considered for now.
# NOTE Ambiguities are not useful in this context. Extensions might be.
# AALetters = [IUPACData.protein_letters_1to3[letter].upper() for letter in IUPACData.protein_letters]
# DNALetters = ['D' + letter for letter in IUPACData.unambiguous_dna_letters]
# RNALetters = [letter for letter in IUPACData.unambiguous_rna_letters]
#
# class AAAlphabet(IUPAC.IUPACProtein):
#     letters = AALetters
#     size = 3
#     pass
#
# class DNAAlphabet(IUPAC.IUPACUnambiguousDNA):
#     letters = DNALetters
#     size = 2
#     pass
#
# class RNAAlphabet(IUPAC.IUPACUnambiguousRNA):
#     letters = RNALetters
#     size = 1
#     pass
#
# class MoleculeAlphabet(Alphabet.Alphabet):
#     letters = AALetters + DNALetters + RNALetters
#     pass
#
# AAalphabet = AAAlphabet()
# DNAalphabet = DNAAlphabet()
# RNAalphabet = RNAAlphabet()
# alphabet = MoleculeAlphabet()

# TODO should this be in the topDB?
def polymer_type(alphabet):
    if isinstance(alphabet, Alphabet.ProteinAlphabet):
        return "AA"
    elif isinstance(alphabet, Alphabet.DNAAlphabet):
        return "DNA"
    elif isinstance(alphabet, Alphabet.RNAAlphabet):
        return "RNA"
    else:
        print(alphabet)
        assert False # TODO helpful message
    return

# NOTE default values to determine edges for secondary structure
# TODO improve and cite these
alpha_phi = - 60 * pi / 180.
alpha_psi = - 45 * pi / 180.
beta_phi = -135 * pi / 180.
beta_psi = 135 * pi / 180.
hbond = 1.9  # hydrogen - acceptor
hangle = 180. * pi / 180.  # angle donor-hydrogen-acceptor

# TODO default omega value

defaultTopologyDB = None
from .config import data_path
if os.path.isfile(data_path + "amber99sb-ildn.xml"):
    from .fileio import read_topology_database
    defaultTopologyDB = read_topology_database("amber99sb-ildn")
