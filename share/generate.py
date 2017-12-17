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
# Last Modified : Sun 17 Dec 2017 11:42:38 PM CET
#
#####################################

import MOBi

forcefields = ['amber99sb-ildn']

# TODO add RNA
for ff in forcefields:
    topDB = MOBi.tools.ffparsergmx.generate_chemical_primary_edge_database(ff, MOBi.data.PDBReducedProtein)

    databaseName = ff + "_protein"
    MOBi.fileio.write_primary_edge_database(topDB, databaseName)

    topDB2 = MOBi.tools.ffparsergmx.generate_chemical_primary_edge_database(ff, MOBi.data.PDBRNAalphabet)

    databaseName = ff + "_rna"
    MOBi.fileio.write_primary_edge_database(topDB, databaseName)

    topDB.update(topDB2)

    databaseName = ff
    MOBi.fileio.write_primary_edge_database(topDB, databaseName)

    pass
