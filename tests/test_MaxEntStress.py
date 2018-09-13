#!/usr/bin/env python3
#####################################
#
# Filename : test_MaxEntStress.py
#
# Projectname : diSTruct
#
# Author : Oskar Taubert
#
# Creation Date : Tue 17 Apr 2018 04:25:49 PM CEST
#
# Last Modified : Thu 13 Sep 2018 01:05:46 PM CEST
#
#####################################

from pytest import approx

from Bio.PDB import PDBParser
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet.IUPAC import unambiguous_rna

from distruct import Superimposer
from distruct import Distructure
from distruct.tools.pdb import get_contacts


from distruct import config
testFilePath = config.data_path + 'tests/'


def test_maxent_from_contacts():
    # TODO use different molecule here, that has all AAs
    code = '1ptq'

    fileName = testFilePath + code + '.pdb'

    refStructure = PDBParser().get_structure(code, fileName)

    # NOTE remove all unneeded residues to make sure, atom naming is consistent
    rToRemove = [r for r in refStructure.get_residues() if r.get_id()[0] != ' ']
    for r in rToRemove:
        refStructure[0]['A'].detach_child(r.get_id())
        pass
    contacts = get_contacts(refStructure[0], cutOff=5., minSeqDist=0)

    sequences = []
    with open(fileName, 'rU') as f:
        sequences = [r.seq for r in SeqIO.parse(f, "pdb-seqres")]
        pass

    ds = Distructure(
            'test',
            sequences,
            [
                [r.get_id() for r in c if r.get_id()[0] == ' ']
                for c in refStructure[0]
            ]
    )

    ds.generate_primary_contacts()
    ds.set_tertiary_contacts(contacts)
    ds.run()

    sup = Superimposer()
    sup.set_structures(refStructure, ds)

    RMSD = sup.rms
    # TODO
    print("this test sometimes fails depending on initialization.")
    print("test with deterministic initialization is coming.")
    assert RMSD < 0.15
    return


def test_RNA():

    code = "2gis"
    fileName = testFilePath + code + '.pdb'

    refStructure = PDBParser().get_structure(code, fileName)

    rToRemove = [r for r in refStructure.get_residues() if r.get_id()[0] != ' ']
    for r in rToRemove:
        refStructure[0]['A'].detach_child(r.get_id())
        pass
    contacts = get_contacts(refStructure[0], cutOff=6., minSeqDist=0)

    sequences = []
    for chain in refStructure[0]:
        s = ''.join([r.get_resname().strip() for r in chain if r.get_id()[0]==' '])
        s = Seq(s, unambiguous_rna)
        sequences.append(s)
        pass

    ds = Distructure('test', sequences)

    ds.generate_primary_contacts()
    ds.set_tertiary_contacts(contacts)
    ds.run()

    sup = Superimposer()
    sup.set_structures(refStructure, ds)

    # from Bio.PDB.PDBIO import PDBIO
    # io = PDBIO()
    # io.set_structure(ds)
    # io.save("test.pdb")

    RMSD = sup.rms

    # TODO
    print("this test sometimes fails depending on initialization.")
    print("test with deterministic initialization is coming.")
    assert RMSD < 6.
    return
