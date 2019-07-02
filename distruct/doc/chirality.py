#!/usr/bin/env python3
#####################################
#
# Filename : chirality.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Thu 25 Apr 2019 03:25:54 PM CEST
#
# Last Modified : Tue 02 Jul 2019 05:58:03 PM CEST
#
#####################################

import distruct as ds

# NOTE this is a test whether diSTruct generates consistent chirality when given no tertiary contacts whatsoever

# NOTE read a sequence
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein
from Bio.PDB import calc_dihedral

seqs = [Seq('HRFKVYNYMSPTFCDHCGSLLWGLVKQGLKCEDCGMNVHHKCREKVANLC', generic_protein)]

s = ds.Distructure('1ptq', seqs)
s.generate_primary_contacts()
s.run()

for r in s.get_residues():
    if r.get_resname() != 'GLY':
        v1 = r['CA'].get_vector()
        v2 = r['N'].get_vector()
        v3 = r['C'].get_vector()
        v4 = r['CB'].get_vector()
        theta = calc_dihedral(v1, v2, v3, v4)
        if theta < 0.:
            print("dihedral error at", r)
        else:
            print(r)
        pass
    pass

# NOTE write result to file
from Bio.PDB.mmcifio import MMCIFIO
io = MMCIFIO()
io.set_structure(s)
io.save('chirality_test.cif')
