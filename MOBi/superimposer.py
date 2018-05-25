#!/usr/bin/env python3
#####################################
#
# Filename : superimposer.py
#
# Projectname : MOBi
#
# Author : Oskar Taubert
#
# Creation Date : Thu 17 May 2018 07:39:39 PM CEST
#
# Last Modified : Wed 23 May 2018 03:47:12 PM CEST
#
#####################################


import numpy as np

from Bio.PDB.PDBExceptions import PDBException


# TODO put this in C? might be necessary for large structures.
class Superimposer(object):
    """
    Rotate/translate/mirror one set of atoms op top of another to minimize RMSD.

    The same as in the Biopython package, but including mirror images (i.e.
    rotations with det -1).

    Reference:
    Gene H. Golub, Charles F. Van Loan: Matrix computations, 2nd edition (1989).

    """


    def __init__(self):
        self.rotran = None
        self.rms = None
        return


    def set_coords(self, fixedCoords, movingCoords):
        """Compute transformation minimizing RMSD."""

        if not len(fixedCoords) == len(movingCoords):
            raise PDBException("Lists of fixed and moving coordinates differ in size")
        l = len(fixedCoords)

        av1 = sum(movingCoords) / l
        av2 = sum(fixedCoords) / l

        movingCoords = movingCoords - av1
        fixedCoords = fixedCoords - av2

        A = np.dot(np.transpose(movingCoords), fixedCoords)
        U, S, Vt = np.linalg.svd(A)
        rot = np.transpose(np.dot(np.transpose(Vt), np.transpose(U)))
        tran = av2 - np.dot(av1, rot)

        # movingCoords = np.dot(movingCoords, rot) + tran
        movingCoords = np.dot(movingCoords, rot)

        diff = fixedCoords - movingCoords
        self.rms = np.sqrt(sum(sum(diff * diff)) / l)
        self.rotran = rot, tran
        return


    def set_atoms(self, fixed, moving):
        """Extract coordinates from atoms and continue with coordinates"""

        if not len(fixed) == len(moving):
            raise PDBException("Lists of fixed and moving atoms differ in size")
        l = len(fixed)
        fixedCoords = np.zeros((l, 3))
        movingCoords = np.zeros((l, 3))
        for i in range(l):
            fixedCoords[i] = fixed[i].get_coord()
            movingCoords[i] = moving[i].get_coord()
            pass

        self.set_coords(fixedCoords, movingCoords)
        return

    def apply(self, atom_list):
        """Apply transformation to atom_list."""
        if self.rotran is None:
            raise PDBException("No transformation has been calculated yet.")
        rot, tran = self.rotran
        rot = rot.astype('f')
        tran = tran.astype('f')

        for atom in atom_list:
            atom.transform(rot, tran)
            pass

        return