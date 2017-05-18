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
# Last Modified : Thu 18 May 2017 18:06:27 CEST
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

AAs = []
AAts = ['C' + AA for AA in AAs] + ['N' + AA for AA in AAs]
# TODO
RAs = []
