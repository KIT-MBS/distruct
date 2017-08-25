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
# Last Modified : Tue 15 Aug 2017 02:27:13 PM CEST
#
#####################################

import os

# TODO fix this in setup? don't use environment variable instead fix this at installation
data_path = os.environ['MOBi_DATAPATH']

# TODO do a setup function where it reads this from gromacs output
# (or tries to  and if it does not work, warns and tells user to set it themselves)
# gromacs_data_prefix = '/usr/'
gromacs_topology_path = '/usr/share/gromacs/top/'
