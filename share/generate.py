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
# Last Modified : Mon 18 Jun 2018 04:19:37 PM CEST
#
#####################################

from MOBi.tools.ffparsergmx import generate
from MOBi.data import AAalphabet, DNAalphabet, RNAalphabet, alphabet
from MOBi.fileio import write_topology_database

ff = "amber99sb-ildn"
topPath = "/usr/share/gromacs/top/"

AAtopDB = generate(ff, AAalphabet, topPath)
RNAtopDB = generate(ff, RNAalphabet, topPath)
DNAtopDB = generate(ff, DNAalphabet, topPath)

topDB = {}
topDB.update(AAtopDB)
topDB.update(RNAtopDB)
topDB.update(DNAtopDB)

write_topology_database(topDB, ff)
