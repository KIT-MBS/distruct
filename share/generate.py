#!/usr/bin/env python3
#####################################
#
# Filename : generate.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Tue 05 Dec 2017 08:08:18 PM CET
#
# Last Modified : Fri 22 Jun 2018 05:27:22 PM CEST
#
#####################################

from MOBi.tools.ffparsergmx import generate
from MOBi.fileio import write_topology_database

from Bio.Alphabet.IUPAC import IUPACProtein, IUPACUnambiguousDNA, IUPACUnambiguousRNA

ff = "amber99sb-ildn"
topPath = "/usr/share/gromacs/top/"

alphabets = [IUPACProtein, IUPACUnambiguousDNA, IUPACUnambiguousRNA]

topDB = generate(ff, alphabets, topPath)

write_topology_database(topDB, ff)
