#!/usr/bin/env python3
#####################################
#
# Filename : example.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Mon 01 Apr 2019 04:08:39 PM CEST
#
# Last Modified : Tue 02 Apr 2019 07:39:25 PM CEST
#
#####################################

import distruct as ds

# NOTE read a sequence
from Bio.Seq import Seq
from Bio.Alphabet import generic_protein

seqs = [Seq('HRFKVYNYMSPTFCDHCGSLLWGLVKQGLKCEDCGMNVHHKCREKVANLC', generic_protein)]

# NOTE read a list of contacts
import os
distFile = os.path.dirname(__file__) + '/../share/1ptq.dist'
contacts = []
with open(distFile, 'r') as f:
    for line in f:
        (c1, r1, a1, c2, r2, a2, distance, weight) = line.split()
        c = (((c1, int(r1), a1), (c2, int(r2), a2)), float(distance), float(weight))
        contacts.append(c)

# NOTE generate the structure
s = ds.Distructure('1ptq', seqs)
s.generate_primary_contacts()
s.set_tertiary_contacts(contacts)
s.run()

# NOTE write result to file
from Bio.PDB.mmcifio import MMCIFIO
io = MMCIFIO()
io.set_structure(s)
io.save('1ptq_out.cif')
