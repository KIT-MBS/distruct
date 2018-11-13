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
# Last Modified : Tue 13 Nov 2018 11:57:36 PM CET
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
    print("An error occured when trying to import the diSTruct extension.")
    print("This may happen when the NetworKit extension could not be found.")
    print("export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:" + os.path.split(_NetworKit.__file__)[0])
    print("or put it in your .bashrc .")
    print("and make sure all the dependencies are installed correctly.")
    print("####################################################################################")
    raise ImportError
