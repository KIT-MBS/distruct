#!/usr/bin/env python3
#####################################
#
# Filename : __init__.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Tue 09 May 2017 13:33:38 CEST
#
# Last Modified : Mon 22 May 2017 14:28:17 CEST
#
#####################################

import os
from . import tools

# TODO fix this in setup
data_path = os.environ['MOBi_DATAPATH']
# gromacs_data_prefix = '/usr/'
# TODO do a setup function where it reads this from gromacs output
# (or tries to  and if it does not work, warns and tells user to set it themselves)
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

