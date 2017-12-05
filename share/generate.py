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
# Last Modified : Tue 05 Dec 2017 08:46:20 PM CET
#
#####################################

import MOBi

forcefields = ['amber99sb-ildn']

# TODO add RNA
for ff in forcefields:
    topDB = MOBi.tools.ffparsergmx.generate_chemical_primary_edge_database(ff, MOBi.data.PDBReducedProtein)

    databaseName = ff + "_protein"
    MOBi.fileio.write_primary_edge_database(topDB, databaseName)
    # test

    topDB2 = MOBi.fileio.parse_primary_edge_database(databaseName, './')

    # TODO add a real check / test case
    # assert topDB == topDB2
    # for k in topDB.keys():
    #     print(k)
    #     for k2 in topDB[k].keys():
    #         print(k2)
    #         print('topDB')
    #         print(topDB[k][k2])
    #         print('topDB2')
    #         print(topDB2[k][k2])
    #         print(topDB[k][k2] == topDB2[k][k2])
    #         pass
    #     pass

    pass
