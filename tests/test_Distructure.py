#!/usr/bin/env python3
#####################################
#
# Filename : test_Distructure.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Mon 11 Jun 2018 04:04:59 PM CEST
#
# Last Modified : Wed 20 Jun 2018 01:50:36 PM CEST
#
#####################################

from MOBi import Distructure
from MOBi import data


def test_Distructure():

    code = "1ptq"
    fileName = code + '.pdb'

    # generate distruct
    # create structure
    ds = Distructure(code, sequences)
    # create edges
    # create primary edges
    ds.generate_primary_edges()
    # create tertiary edges
    # get contact map
    contactMap = MOBi.tools.generate_atomic_contacts(refStructure)
    ds.generate_tertiary_edges(contactMap)
    # generate coordinates
    ds.run()

    # TODO structural alignment
    RMSD = MOBi.tools.align(ds, refStructure)

    return
