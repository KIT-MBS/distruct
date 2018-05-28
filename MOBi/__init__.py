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
# Last Modified : Mon 28 May 2018 11:10:18 AM CEST
#
#####################################

from . import config
from . import data
from . import seq
from . import fileio
from . import tools

from .superimposer import Superimposer

import sys
import os
# TODO is this too hacky? there probably is a proper way
sys.path.append(os.path.dirname(__file__) + "/../lib/")

# TODO convenience imports here
# from _MOBi import BioMaxentStress

# TODO remove this
from _MOBi import doublyWrappedMaxent

from _MOBi import Distruct
