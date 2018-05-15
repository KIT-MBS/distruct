#!/usr/bin/env python3
#####################################
#
# Filename : config.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Mon 22 May 2017 17:23:27 CEST
#
# Last Modified : Wed 02 May 2018 04:09:30 PM CEST
#
#####################################

import os

data_path = os.path.split(os.path.abspath(__file__))[0] + '/../share/'

# TODO do a setup function where it reads this from gromacs output
# (or tries to  and if it does not work, warns and tells user to set it themselves)
# gromacs_data_prefix = '/usr/'
gromacs_topology_path = '/usr/share/gromacs/top/'
