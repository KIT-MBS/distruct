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
# Last Modified : Thu 28 Jun 2018 01:17:20 PM CEST
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
# sys.path.append(os.path.dirname(__file__) + "/../lib/")

from _MOBi import Distructure
