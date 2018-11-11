#!/usr/bin/env python3
#####################################
#
# Filename : __init__.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Tue 09 May 2017 13:33:38 CEST
#
# Last Modified : Sun 11 Nov 2018 07:45:47 PM CET
#
#####################################

from . import config
from . import data
from . import fileio
from . import tools

from .superimposer import Superimposer

import sys
import os

name = "distruct"


try:
    from _diSTruct import Distructure
except ImportError:
    # TODO set the LD_LIBRARY_PATH automatically
    import _NetworKit
    print("####################################################################################")
    print("export LD_LIBRARY_PATH=LD_LIBRARY_PATH:" + os.path.split(_NetworKit.__file__)[0])
    print("or put it in the .bashrc")
    print("i know this is bananas")
    print("i will fix it eventually")
    print("also make sure all the dependencies are installed correctly")
    print("####################################################################################")
    raise ImportError
