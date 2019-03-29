#!/usr/bin/env python3
#####################################
#
# Filename : config.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Mon 22 May 2017 17:23:27 CEST
#
# Last Modified : Fri 29 Mar 2019 06:33:10 PM CET
#
#####################################

import os

data_path = os.path.split(os.path.abspath(__file__))[0] + '/share/'

# TODO do a setup function where it reads this from gromacs output
# (or tries to  and if it does not work, warns and tells user to set it themselves)
# gromacs_data_prefix = '/usr/'
gromacs_topology_path = '/usr/share/gromacs/top/'
