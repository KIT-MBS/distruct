#!/usr/bin/env python3
#####################################
#
# Filename : config.py
#
# Projectname :
#
# Author : Oskar Taubert
#
# Creation Date : Mon 22 May 2017 17:23:27 CEST
#
# Last Modified : Mon 22 May 2017 17:25:33 CEST
#
#####################################

import os

# TODO fix this in setup? don't use environment variable instead fix this at installation
data_path = os.environ['MOBi_DATAPATH']

# TODO do a setup function where it reads this from gromacs output
# (or tries to  and if it does not work, warns and tells user to set it themselves)
# gromacs_data_prefix = '/usr/'
gromacs_topology_path = '/usr/share/gromacs/top/'

AAs = ['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN', 'ARG', 'HIS', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 'PRO', 'CYS', 'MET']
NAs = ['A', 'G', 'C', 'U']

# TODO put in translation for ff names and other conventions
# NOTE amber knows HISD, HISE and HISH but not HIS (and calls them HID, HIE and HIP)
amberAAs = ['ALA', 'GLY', 'SER', 'THR', 'LEU', 'ILE', 'VAL', 'ASN', 'GLN', 'ARG', 'HIP', 'TRP', 'PHE', 'TYR', 'GLU', 'ASP', 'LYS', 'PRO', 'CYS', 'MET']
amberAAts = ['C' + AA for AA in amberAAs] + ['N' + AA for AA in amberAAs]
# TODO amberNAs

# TODO gromos
# TODO charmm
