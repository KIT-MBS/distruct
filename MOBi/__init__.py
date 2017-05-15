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
# Last Modified : Fri 12 May 2017 15:01:04 CEST
#
#####################################

import os
from . import tools

# TODO fix this in setup
data_path = os.environ['MOBi_DATAPATH']
gromacs_data_prefix = '/usr/'
