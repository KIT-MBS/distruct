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
# Last Modified : Wed 11 Mar 2020 05:35:37 PM CET
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
    import networkit
    import _NetworKit

    path_to_networkit = os.path.split(_NetworKit.__file__)[0]
    os.environ['PATH'] = path_to_networkit + os.pathsep + os.environ['PATH']
    from _diSTruct import Distructure
