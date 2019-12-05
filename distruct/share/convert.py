#!/usr/bin/env python3
#####################################
#
# Filename : convert.py
#
# Projectname :
#
# Author : Oskar Taubert
#
# Creation Date : Mon 02 Dec 2019 07:41:07 PM CET
#
# Last Modified : Thu 05 Dec 2019 05:03:18 PM CET
#
#####################################

import distruct
from Bio import Alphabet
# from Bio.Alphabet.IUPAC import extended_protein

# alphabets = [extended_protein]

# NOTE convert default topology database to AA reduced and AA backbone only versions
tdb = distruct.fileio.read_topology_database("amber99sb-ildn")

# NOTE for this only include amino acids
for a in['DNA', 'RNA']:
    for letter in tdb['alphabets'][a]:
        r = tdb['alphabets'][a][letter]
        del tdb[r]
    del tdb['alphabets'][a]

backbone_atoms = {'N', 'CA', 'C'}
reduced_atoms = {'C', 'O'} | backbone_atoms

# NOTE do heavy atoms only version (i.e. no hydrogens)
for letter in tdb['alphabets']['AA']:
    r = tdb['alphabets']['AA'][letter]
    tdb[r]['vertices'] = [v for v in tdb[r]['vertices'] if v[1] is not 'H']

    # TODO delete all edges that don't have all the required atoms
tdb['alphabets']['AA']['X'] = 'XAA'
tdb['XAA'] = tdb['GLY']

distruct.fileio.write_topology_database(tdb, "amber99sb-ildn_heavy")


# NOTE include only backbone, CB and carboxyl oxygen
for letter in tdb['alphabets']['AA']:
    r = tdb['alphabets']['AA'][letter]

    tdb[r]['vertices'] = [v for v in tdb[r]['vertices'] if v[0] in reduced_atoms]

    # TODO delete all edges that don't have all the required atoms
distruct.fileio.write_topology_database(tdb, "amber99sb-ildn_reduced")

# NOTE include only backbone
for letter in tdb['alphabets']['AA']:
    r = tdb['alphabets']['AA'][letter]

    tdb[r]['vertices'] = [v for v in tdb[r]['vertices'] if v[0] in backbone_atoms]
    # TODO delete all edges that don't have all the required atoms
distruct.fileio.write_topology_database(tdb, "amber99sb-ildn_backbone")
