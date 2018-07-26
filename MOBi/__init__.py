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
# Last Modified : Thu 26 Jul 2018 07:04:43 PM CEST
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


try:
    from _MOBi import Distructure
except ImportError:
    # TODO set the LD_LIBRARY_PATH automatically
    import _NetworKit
    print("export LD_LIBRARY_PATH=LD_LIBRARY_PATH:" + os.path.split(_NetworKit.__file__)[0])
    print("or put it in the .bashrc")
    print("i know this is bananas")
    print("i will fix it eventually")
    raise ImportError
